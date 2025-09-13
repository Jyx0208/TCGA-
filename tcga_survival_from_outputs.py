#!/usr/bin/env python3
"""
Fetch TCGA survival for patients listed by tcga_multimodal_gdc_fetch.py.

Inputs (from the existing pipeline's outdir):
  - patients_complete.txt  (one 12-char patient ID per line), or
  - patient_file_map.csv   (will use the 'patient' column if the txt is missing)

Outputs:
  - survival.csv with columns:
      patient, os_days, os_months, os_event, gender, age_at_diagnosis_years, tumor_stage, case_id

Notes:
  - Uses open-access GDC /cases endpoint (no token required).
  - If multiple diagnoses per case, picks the one with the longest available follow-up
    (max of coalesce(days_to_death, days_to_last_follow_up)).
"""

import argparse
import json
import os
import time
from typing import Dict, List, Iterable, Any

import pandas as pd
import requests

GDC_API = "https://api.gdc.cancer.gov"
HEADERS = {"Content-Type": "application/json"}

def _post_gdc(path: str, payload: Dict, as_json: bool = True,
              max_retries: int = 6, timeout: int = 120):
    url = f"{GDC_API}{path}"
    backoff = 1.5
    for attempt in range(max_retries):
        try:
            r = requests.post(url, headers=HEADERS, data=json.dumps(payload), timeout=timeout)
        except requests.RequestException:
            if attempt == max_retries - 1:
                raise
            time.sleep(backoff ** attempt)
            continue
        if r.status_code == 200:
            return r.json() if as_json else r.text
        if r.status_code in (429, 500, 502, 503, 504):
            time.sleep(backoff ** attempt)
            continue
        r.raise_for_status()
    raise RuntimeError(f"GDC API retry limit reached for {path}")

def _patient12(series: pd.Series) -> pd.Series:
    return series.astype(str).str.slice(0, 12)

def _chunked(seq: List[Any], n: int) -> Iterable[List[Any]]:
    for i in range(0, len(seq), n):
        yield seq[i:i+n]

def _coalesce(*vals):
    for v in vals:
        if v is not None:
            return v
    return None

def _read_patient_list(outdir: str, patients_file: str | None) -> List[str]:
    # Priority: explicit file > patients_complete.txt > patient_file_map.csv
    if patients_file and os.path.isfile(patients_file):
        s = pd.read_csv(patients_file, header=None, names=["patient"])
        return sorted(_patient12(s["patient"]).dropna().unique())
    txt_path = os.path.join(outdir, "patients_complete.txt")
    if os.path.isfile(txt_path):
        s = pd.read_csv(txt_path, header=None, names=["patient"])
        return sorted(_patient12(s["patient"]).dropna().unique())
    map_path = os.path.join(outdir, "patient_file_map.csv")
    if os.path.isfile(map_path):
        s = pd.read_csv(map_path)
        return sorted(_patient12(s["patient"]).dropna().unique())
    raise FileNotFoundError(
        "Could not find patients_complete.txt or patient_file_map.csv. "
        "Please provide --patients or run the upstream script first."
    )

def _query_cases_by_patients(project: str, patients: List[str]) -> pd.DataFrame:
    """
    Query GDC /cases for the given 12-char submitter IDs, constrained to a project.
    Flattens demographics and best-available diagnosis per case.
    """
    rows = []
    fields = ",".join(sorted(set([
        "case_id",
        "submitter_id",
        "demographic.gender",
        # 年龄应来自 diagnoses，而不是 demographic
        "diagnoses.age_at_diagnosis",
        "diagnoses.vital_status",
        "diagnoses.days_to_death",
        "diagnoses.days_to_last_follow_up",
        "diagnoses.tumor_stage",
    ])))
    for chunk in _chunked(patients, 200):
        payload = {
            "filters": {
                "op": "and",
                "content": [
                    {"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}},
                    {"op": "in", "content": {"field": "cases.submitter_id", "value": chunk}},
                ],
            },
            "fields": fields,
            "size": max(1000, len(chunk)),
        }
        data = _post_gdc("/cases", payload, as_json=True)
        hits = (data or {}).get("data", {}).get("hits", [])
        for h in hits:
            base = {
                "case_id": h.get("case_id"),
                "patient": h.get("submitter_id"),
                "gender": (h.get("demographic") or {}).get("gender"),
            }
            diags = h.get("diagnoses") or []
            # 选“随访/死亡时间最长”的诊断条目
            def diag_score(d: Dict) -> float:
                dt = d.get("days_to_death")
                lf = d.get("days_to_last_follow_up")
                v = dt if dt is not None else lf
                return float(v) if v is not None else -1.0

            d_best = max(diags, key=diag_score) if diags else {}

            row = dict(base)
            row["vital_status"] = d_best.get("vital_status")
            row["days_to_death"] = d_best.get("days_to_death")
            row["days_to_last_follow_up"] = d_best.get("days_to_last_follow_up")
            row["tumor_stage"] = d_best.get("tumor_stage")
            # 诊断年龄（天）
            row["age_at_diagnosis_days"] = d_best.get("age_at_diagnosis")

            # 计算 OS
            vs = (row["vital_status"] or "").strip().lower()
            if vs in {"dead", "deceased"}:
                row["os_days"] = d_best.get("days_to_death")
                row["os_event"] = 1
            else:
                row["os_days"] = d_best.get("days_to_last_follow_up")
                row["os_event"] = 0

            rows.append(row)

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    # ---- 关键：在 round 之前把列强制为数值 ----
    for col in ["os_days", "age_at_diagnosis_days"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # 去重：同一 patient 保留 os_days 最长的一行
    df["os_days_rank"] = df["os_days"].fillna(-1).astype(float)
    df = df.sort_values(["patient", "os_days_rank"], ascending=[True, False]).drop_duplicates("patient")
    df = df.drop(columns=["os_days_rank"])

    # 衍生列
    df["os_months"] = (df["os_days"] / 30.44).round(2)
    df["age_at_diagnosis_years"] = (df["age_at_diagnosis_days"] / 365.25).round(2)

    return df[
        ["patient", "os_days", "os_months", "os_event",
         "gender", "age_at_diagnosis_years", "tumor_stage", "case_id"]
    ]


def main():
    ap = argparse.ArgumentParser(description="Fetch survival for patients produced by tcga_multimodal_gdc_fetch.py")
    ap.add_argument("--project", default="TCGA-LIHC", help="GDC project ID (default: TCGA-LIHC)")
    ap.add_argument("--outdir", default="out", help="Directory where upstream outputs live and where survival.csv will be written")
    ap.add_argument("--patients", default=None, help="Optional path to a text/CSV file of patient IDs (12-char). If omitted, will use patients_complete.txt in --outdir")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    patients = _read_patient_list(args.outdir, args.patients)
    if not patients:
        raise SystemExit("No patients found to query.")

    print(f"[1/3] Patients to query: {len(patients)}")
    print(f"[2/3] Querying GDC /cases for survival fields in project {args.project} …", flush=True)
    df = _query_cases_by_patients(args.project, patients)

    out_path = os.path.join(args.outdir, "survival.csv")
    df.to_csv(out_path, index=False)
    print(f"[3/3] Saved survival table → {out_path}  (rows={len(df)})")

if __name__ == "__main__":
    main()
