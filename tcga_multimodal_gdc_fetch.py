#!/usr/bin/env python3
"""
TCGA multi‑modal fetcher (Scheme A, GDC API)
-------------------------------------------
Purpose
  • From TCGA-LIHC, find patients who concurrently have RNA (STAR - Counts), CNV (Masked Copy Number Segment),
    and pathology whole-slide images (WSI; Diagnostic Slide, SVS), then export per‑modality manifests for gdc-client.

Deliverables (written into --outdir)
  • rna.tsv, cnv.tsv, wsi_raw.tsv  – raw file inventories
  • patients_complete.txt          – 12-char patient IDs with all three modalities
  • wsi.tsv                        – deduped representative WSI per patient (prefers -DX1, else largest file)
  • patient_file_map.csv           – patient → {rna_file_id, cnv_file_id, wsi_file_id} mapping
  • manifest_rna.txt / manifest_cnv.txt / manifest_wsi.txt – manifests ready for `gdc-client download -m ...`

Usage
  python tcga_multimodal_gdc_fetch.py \
      --project TCGA-LIHC \
      --outdir out \
      [--download]    # if you already have gdc-client in PATH and want to download immediately

Notes
  • Only open-access data categories are used; no dbGaP token required.
  • Primary Tumor (01) only, to keep clinical intent tight for prognosis modeling.
  • Harmonized GDC data; RNA filtered to STAR - Counts to avoid legacy HTSeq.
"""

import argparse
import json
import os
import sys
import time
from typing import Dict, List

import pandas as pd
import requests

GDC_API = "https://api.gdc.cancer.gov"
HEADERS = {"Content-Type": "application/json"}


def _post_gdc(path: str, payload: Dict, as_json: bool = True, max_retries: int = 6, timeout: int = 120):
    url = f"{GDC_API}{path}"
    backoff = 1.5
    for attempt in range(max_retries):
        try:
            r = requests.post(url, headers=HEADERS, data=json.dumps(payload), timeout=timeout)
        except requests.RequestException as e:
            if attempt == max_retries - 1:
                raise
            time.sleep(backoff ** attempt)
            continue
        if r.status_code == 200:
            return r.json() if as_json else r.text
        if r.status_code in (429, 500, 502, 503, 504):
            # Respect server load; exponential backoff
            time.sleep(backoff ** attempt)
            continue
        # Hard error
        r.raise_for_status()
    raise RuntimeError(f"GDC API retry limit reached for {path}")


def _query_files(filters: Dict, fields: List[str], size: int = 10000) -> pd.DataFrame:
    payload = {
        "filters": filters,
        "fields": ",".join(sorted(set(fields))),
        "size": size,
    }
    data = _post_gdc("/files", payload, as_json=True)
    hits = (data or {}).get("data", {}).get("hits", [])
    # Flatten nested structures of interest
    rows = []
    for h in hits:
        base = {
            "file_id": h.get("file_id"),
            "file_name": h.get("file_name"),
            "file_size": h.get("file_size"),
            "data_format": h.get("data_format"),
        }
        if isinstance(h.get("analysis"), dict):
            base["analysis.workflow_type"] = h["analysis"].get("workflow_type")
        submitters = {c.get("submitter_id") for c in h.get("cases", []) if c.get("submitter_id")}
        if not submitters:
            submitters = {None}
        for sid in submitters:
            row = dict(base)
            row["cases.submitter_id"] = sid
            rows.append(row)
    return pd.DataFrame(rows)


def _manifest_from_file_ids(file_ids: List[str]) -> str:
    payload = {"ids": list(map(str, file_ids))}
    return _post_gdc("/manifest", payload, as_json=False)


# ----- Filters (Scheme A) -----

def _filters_rna(project: str) -> Dict:
    return {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}},
            {"op": "in", "content": {"field": "files.data_category", "value": ["Transcriptome Profiling"]}},
            {"op": "in", "content": {"field": "files.data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
            {"op": "in", "content": {"field": "cases.samples.sample_type", "value": ["Primary Tumor"]}},
        ],
    }


def _filters_cnv(project: str) -> Dict:
    return {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}},
            {"op": "in", "content": {"field": "files.data_category", "value": ["Copy Number Variation"]}},
            {"op": "in", "content": {"field": "files.data_type", "value": ["Masked Copy Number Segment"]}},
            {"op": "in", "content": {"field": "cases.samples.sample_type", "value": ["Primary Tumor"]}},
        ],
    }


def _filters_wsi(project: str) -> Dict:
    return {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}},
            {"op": "in", "content": {"field": "files.data_category", "value": ["Biospecimen"]}},
            {"op": "in", "content": {"field": "files.data_type", "value": ["Slide Image"]}},
            {"op": "in", "content": {"field": "files.experimental_strategy", "value": ["Diagnostic Slide"]}},
            {"op": "in", "content": {"field": "files.data_format", "value": ["SVS"]}},
        ],
    }


# ----- Utilities -----

def _patient12(series: pd.Series) -> pd.Series:
    return series.astype(str).str.slice(0, 12)


def _pick_one_per_patient(df: pd.DataFrame, prefer_col: str = "file_size") -> pd.DataFrame:
    # Keep exactly one file per patient; tie-breaker by descending prefer_col
    out = df.copy()
    out["patient"] = _patient12(out["cases.submitter_id"])
    out = out.sort_values(["patient", prefer_col], ascending=[True, False])
    out = out.groupby("patient", as_index=False).first()
    return out


def _pick_wsi(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["patient"] = _patient12(out["cases.submitter_id"])
    out["is_dx1"] = out["file_name"].fillna("").str.contains("-DX1", case=False, regex=False)
    out = out.sort_values(["patient", "is_dx1", "file_size"], ascending=[True, False, False])
    out = out.groupby("patient", as_index=False).first()
    out = out.drop(columns=["is_dx1"])  # clean
    return out


# ----- Main -----

def main():
    parser = argparse.ArgumentParser(description="Fetch GDC files and manifests for TCGA multi-modal (RNA/CNV/WSI)")
    parser.add_argument("--project", default="TCGA-LIHC", help="GDC project ID (default: TCGA-LIHC)")
    parser.add_argument("--outdir", default="out", help="Output directory")
    parser.add_argument("--download", action="store_true", help="If set, will run gdc-client download for each manifest if available in PATH")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 1) Query three modalities
    print("[1/6] Querying RNA (STAR - Counts)…", flush=True)
    rna_df = _query_files(_filters_rna(args.project), [
        "file_id",
        "file_name",
        "file_size",
        "analysis.workflow_type",
        "cases.submitter_id",
    ])
    rna_path = os.path.join(args.outdir, "rna.tsv")
    rna_df.to_csv(rna_path, sep="\t", index=False)
    print(f"  → RNA files: {len(rna_df)}  saved: {rna_path}")

    print("[2/6] Querying CNV (Masked Copy Number Segment)…", flush=True)
    cnv_df = _query_files(_filters_cnv(args.project), [
        "file_id",
        "file_name",
        "file_size",
        "cases.submitter_id",
    ])
    cnv_path = os.path.join(args.outdir, "cnv.tsv")
    cnv_df.to_csv(cnv_path, sep="\t", index=False)
    print(f"  → CNV files: {len(cnv_df)}  saved: {cnv_path}")

    print("[3/6] Querying WSI (Diagnostic Slide, SVS)…", flush=True)
    wsi_df = _query_files(_filters_wsi(args.project), [
        "file_id",
        "file_name",
        "file_size",
        "data_format",
        "cases.submitter_id",
    ])
    wsi_raw_path = os.path.join(args.outdir, "wsi_raw.tsv")
    wsi_df.to_csv(wsi_raw_path, sep="\t", index=False)
    print(f"  → WSI files: {len(wsi_df)}  saved: {wsi_raw_path}")

    # 2) Intersection at patient level (12-char)
    p_rna = set(_patient12(rna_df["cases.submitter_id"].dropna()).unique())
    p_cnv = set(_patient12(cnv_df["cases.submitter_id"].dropna()).unique())
    p_wsi = set(_patient12(wsi_df["cases.submitter_id"].dropna()).unique())
    patients = sorted(list(p_rna & p_cnv & p_wsi))

    pts_path = os.path.join(args.outdir, "patients_complete.txt")
    with open(pts_path, "w") as f:
        for p in patients:
            f.write(p + "\n")
    print(f"[4/6] Patients with all three modalities: {len(patients)}  saved: {pts_path}")

    # 3) Subset to intersected patients & pick representatives
    rna_sub = rna_df[_patient12(rna_df["cases.submitter_id"]).isin(patients)]
    cnv_sub = cnv_df[_patient12(cnv_df["cases.submitter_id"]).isin(patients)]
    wsi_sub = wsi_df[_patient12(wsi_df["cases.submitter_id"]).isin(patients)]

    rna_rep = _pick_one_per_patient(rna_sub, "file_size")
    cnv_rep = _pick_one_per_patient(cnv_sub, "file_size")
    wsi_rep = _pick_wsi(wsi_sub)

    wsi_path = os.path.join(args.outdir, "wsi.tsv")
    wsi_rep.to_csv(wsi_path, sep="\t", index=False)
    print(f"[5/6] Representative WSI per patient saved: {wsi_path}")

    # 4) Mapping table (patient → per-modality file_ids)
    mapping = pd.DataFrame({"patient": patients})
    mapping = mapping.merge(
        rna_rep[["patient", "file_id", "file_name"]].rename(columns={
            "file_id": "rna_file_id",
            "file_name": "rna_file_name",
        }),
        on="patient",
        how="left",
    )
    mapping = mapping.merge(
        cnv_rep[["patient", "file_id", "file_name"]].rename(columns={
            "file_id": "cnv_file_id",
            "file_name": "cnv_file_name",
        }),
        on="patient",
        how="left",
    )
    mapping = mapping.merge(
        wsi_rep[["patient", "file_id", "file_name"]].rename(columns={
            "file_id": "wsi_file_id",
            "file_name": "wsi_file_name",
        }),
        on="patient",
        how="left",
    )
    map_path = os.path.join(args.outdir, "patient_file_map.csv")
    mapping.to_csv(map_path, index=False)
    print(f"[5.5/6] Mapping CSV saved: {map_path}  (rows={len(mapping)})")

    # 5) Manifests for gdc-client
    mani_rna = _manifest_from_file_ids(list(rna_rep["file_id"].dropna().unique()))
    mani_cnv = _manifest_from_file_ids(list(cnv_rep["file_id"].dropna().unique()))
    mani_wsi = _manifest_from_file_ids(list(wsi_rep["file_id"].dropna().unique()))

    m_rna_path = os.path.join(args.outdir, "manifest_rna.txt")
    m_cnv_path = os.path.join(args.outdir, "manifest_cnv.txt")
    m_wsi_path = os.path.join(args.outdir, "manifest_wsi.txt")

    with open(m_rna_path, "w") as f:
        f.write(mani_rna)
    with open(m_cnv_path, "w") as f:
        f.write(mani_cnv)
    with open(m_wsi_path, "w") as f:
        f.write(mani_wsi)

    print(f"[6/6] Manifests written →\n  {m_rna_path}\n  {m_cnv_path}\n  {m_wsi_path}")

    # Optional: trigger downloads via gdc-client
    if args.download:
        import shutil
        import subprocess

        gdc = shutil.which("gdc-client")
        if not gdc:
            print("[WARN] gdc-client not found in PATH. Skipping downloads.", file=sys.stderr)
        else:
            for mp in (m_rna_path, m_cnv_path, m_wsi_path):
                print(f"[DL] gdc-client download -m {mp}")
                subprocess.run([gdc, "download", "-m", mp], check=True)


if __name__ == "__main__":
    main()
