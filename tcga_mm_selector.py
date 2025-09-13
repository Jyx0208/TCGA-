#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TCGA 多模态数据筛选器：
- 按癌种筛选同时具备 RNA、CNV、WSI 和生存信息的病人
- 支持样本级（root sample）跨模态匹配
- 生成 GDC manifest 以下载各模态文件
"""

import argparse
import os
from typing import List, Dict, Any, Optional

import requests
import pandas as pd

GDC_API = "https://api.gdc.cancer.gov"


def gdc_post(path: str, payload: Dict[str, Any], stream: bool = False) -> requests.Response:
    url = f"{GDC_API}/{path.lstrip('/')}"
    r = requests.post(url, json=payload, stream=stream)
    r.raise_for_status()
    return r


def gdc_get(path: str, params: Dict[str, Any] = None, stream: bool = False) -> requests.Response:
    url = f"{GDC_API}/{path.lstrip('/')}"
    r = requests.get(url, params=params or {}, stream=stream)
    r.raise_for_status()
    return r


def list_tcga_projects() -> pd.DataFrame:
    filters = {
        "op": "in",
        "content": {
            "field": "programs.name",
            "value": ["TCGA"],
        },
    }
    fields = ",".join(
        [
            "project_id",
            "name",
            "primary_site",
            "disease_type",
            "summary.case_count",
        ]
    )
    payload = {"filters": filters, "fields": fields, "format": "JSON", "size": 1000}
    res = gdc_post("/projects", payload).json()
    hits = res.get("data", {}).get("hits", [])
    if not hits:
        return pd.DataFrame()
    df = pd.json_normalize(hits)
    df = df.rename(columns={"summary.case_count": "case_count"})
    return df[
        ["project_id", "name", "primary_site", "disease_type", "case_count"]
    ].sort_values("project_id")


def build_in_clause(field: str, values: List[Any]) -> Dict[str, Any]:
    return {"op": "in", "content": {"field": field, "value": values}}


def get_files_for_modality(
    project_id: str,
    data_category: str,
    data_type: Optional[str] = None,
    workflow_type: Optional[str] = None,
    sample_type: str = "Primary Tumor",
    size: int = 10000,
) -> pd.DataFrame:
    clauses = [
        build_in_clause("cases.project.project_id", [project_id]),
        build_in_clause("files.data_category", [data_category]),
        build_in_clause("cases.samples.sample_type", [sample_type]),
    ]
    if data_type:
        clauses.append(build_in_clause("files.data_type", [data_type]))
    if workflow_type:
        clauses.append(build_in_clause("files.analysis.workflow_type", [workflow_type]))
        clauses_alt = clauses[:-1] + [build_in_clause("files.workflow_type", [workflow_type])]
    else:
        clauses_alt = None

    fields = ",".join(
        [
            "file_id",
            "file_name",
            "data_type",
            "data_category",
            "experimental_strategy",
            "platform",
            "analysis.workflow_type",
            "workflow_type",
            "cases.case_id",
            "cases.submitter_id",
            "cases.samples.sample_id",
            "cases.samples.submitter_id",
            "cases.samples.sample_type",
        ]
    )

    payload = {
        "filters": {"op": "and", "content": clauses},
        "fields": fields,
        "format": "JSON",
        "size": size,
    }

    def post_and_collect(pld):
        res = gdc_post("/files", pld).json()
        hits = res.get("data", {}).get("hits", [])
        rows = []
        for h in hits:
            cases = h.get("cases") or []
            if not cases:
                continue
            c = cases[0]
            case_sid = c.get("submitter_id")
            case_id = c.get("case_id")
            samples = c.get("samples") or []
            samp = None
            for s in samples:
                if s.get("sample_type") == sample_type:
                    samp = s
                    break
            if not samp and samples:
                samp = samples[0]
            rows.append(
                {
                    "file_id": h.get("file_id"),
                    "file_name": h.get("file_name"),
                    "data_type": h.get("data_type"),
                    "data_category": h.get("data_category"),
                    "workflow_type": h.get("analysis", {}).get("workflow_type")
                    or h.get("workflow_type"),
                    "experimental_strategy": h.get("experimental_strategy"),
                    "platform": h.get("platform"),
                    "case_id": case_id,
                    "case_submitter_id": case_sid,
                    "sample_id": (samp or {}).get("sample_id"),
                    "sample_submitter_id": (samp or {}).get("submitter_id"),
                    "sample_type": (samp or {}).get("sample_type"),
                }
            )
        return pd.DataFrame(rows)

    df = post_and_collect(payload)
    if workflow_type and df.empty and clauses_alt is not None:
        payload_alt = {
            "filters": {"op": "and", "content": clauses_alt},
            "fields": fields,
            "format": "JSON",
            "size": size,
        }
        df = post_and_collect(payload_alt)

    return df


def to_root_sample(sample_submitter_id: Optional[str]) -> Optional[str]:
    if not sample_submitter_id:
        return None
    return sample_submitter_id[:15]


def get_survival_for_project(project_id: str) -> pd.DataFrame:
    filters = {"op": "in", "content": {"field": "projects.project_id", "value": [project_id]}}
    fields = ",".join(
        [
            "case_id",
            "submitter_id",
            "diagnoses.vital_status",
            "diagnoses.days_to_death",
            "diagnoses.days_to_last_follow_up",
        ]
    )
    payload = {"filters": filters, "fields": fields, "format": "JSON", "size": 10000}
    res = gdc_post("/cases", payload).json()
    hits = res.get("data", {}).get("hits", [])
    rows = []
    for h in hits:
        case_sid = h.get("submitter_id")
        case_id = h.get("case_id")
        diags = h.get("diagnoses") or []
        max_dod = None
        max_dlfu = None
        vital = None
        for d in diags:
            vital = d.get("vital_status") or vital
            dd = d.get("days_to_death")
            lf = d.get("days_to_last_follow_up")

            def to_num(x):
                try:
                    return int(float(x))
                except Exception:
                    return None

            dd = to_num(dd)
            lf = to_num(lf)
            if dd is not None:
                max_dod = max(max_dod or dd, dd)
            if lf is not None:
                max_dlfu = max(max_dlfu or lf, lf)
        os_event = None
        os_time = None
        if (vital or "").lower() == "dead" and max_dod is not None:
            os_event = 1
            os_time = max_dod
        else:
            if max_dlfu is not None:
                os_event = 0
                os_time = max_dlfu
        rows.append(
            {
                "case_id": case_id,
                "case_submitter_id": case_sid,
                "vital_status": vital,
                "os_time": os_time,
                "os_event": os_event,
            }
        )
    df = pd.DataFrame(rows)
    df = df.dropna(subset=["os_time", "os_event"])
    return df


def pick_one_per_root(df: pd.DataFrame, by_fields: List[str]) -> pd.DataFrame:
    if df.empty:
        return df
    df = df.copy().sort_values(by=["case_submitter_id", "root_sample", "file_name"])
    return df.drop_duplicates(subset=by_fields, keep="first")


def select_complete_cases(
    project_id: str,
    rna_workflow: str = "HTSeq - Counts",
    cnv_data_type: str = "Copy Number Segment",
    sample_type: str = "Primary Tumor",
    require_sample_match: bool = True,
) -> Dict[str, pd.DataFrame]:
    rna_df = get_files_for_modality(
        project_id,
        "Transcriptome Profiling",
        data_type="Gene Expression Quantification",
        workflow_type=rna_workflow,
        sample_type=sample_type,
    )
    cnv_df = get_files_for_modality(
        project_id, "Copy Number Variation", data_type=cnv_data_type, sample_type=sample_type
    )
    wsi_df = get_files_for_modality(
        project_id, "Biospecimen", data_type="Slide Image", sample_type=sample_type
    )
    surv_df = get_survival_for_project(project_id)

    for df in (rna_df, cnv_df, wsi_df):
        if not df.empty:
            df["root_sample"] = df["sample_submitter_id"].apply(to_root_sample)

    set_rna = set(rna_df.case_submitter_id.dropna()) if not rna_df.empty else set()
    set_cnv = set(cnv_df.case_submitter_id.dropna()) if not cnv_df.empty else set()
    set_wsi = set(wsi_df.case_submitter_id.dropna()) if not wsi_df.empty else set()
    set_sur = set(surv_df.case_submitter_id.dropna()) if not surv_df.empty else set()

    case_intersection = set_rna & set_cnv & set_wsi & set_sur

    if not case_intersection:
        return {"rna": rna_df, "cnv": cnv_df, "wsi": wsi_df, "surv": surv_df, "complete": pd.DataFrame()}

    def case_to_roots(df: pd.DataFrame) -> Dict[str, set]:
        m = {}
        if df.empty:
            return m
        for _, row in df.iterrows():
            c = row["case_submitter_id"]
            r = row.get("root_sample")
            if not r:
                continue
            m.setdefault(c, set()).add(r)
        return m

    rna_roots = case_to_roots(rna_df)
    cnv_roots = case_to_roots(cnv_df)
    wsi_roots = case_to_roots(wsi_df)

    selected_rows = []
    for case in sorted(case_intersection):
        roots_rna = rna_roots.get(case, set())
        roots_cnv = cnv_roots.get(case, set())
        roots_wsi = wsi_roots.get(case, set())
        roots_rc = roots_rna & roots_cnv
        if require_sample_match:
            roots_rcw = roots_rc & roots_wsi
            if not roots_rcw:
                continue
            chosen_root = sorted(roots_rcw)[0]
        else:
            chosen_root = sorted(roots_rc)[0] if roots_rc else None
            if not chosen_root:
                continue

        rna_pick = rna_df[(rna_df.case_submitter_id == case) & (rna_df.root_sample == chosen_root)]
        cnv_pick = cnv_df[(cnv_df.case_submitter_id == case) & (cnv_df.root_sample == chosen_root)]
        wsi_pick = wsi_df[(wsi_df.case_submitter_id == case) & (wsi_df.root_sample == chosen_root)]
        if rna_pick.empty or cnv_pick.empty or (require_sample_match and wsi_pick.empty):
            continue
        rna_pick = pick_one_per_root(rna_pick, ["case_submitter_id", "root_sample"])
        cnv_pick = pick_one_per_root(cnv_pick, ["case_submitter_id", "root_sample"])
        wsi_pick = pick_one_per_root(wsi_pick, ["case_submitter_id", "root_sample"])

        surv_pick = surv_df[surv_df.case_submitter_id == case]
        if surv_pick.empty:
            continue
        s = surv_pick.iloc[0]

        selected_rows.append(
            {
                "case_submitter_id": case,
                "root_sample": chosen_root,
                "rna_file_id": rna_pick.iloc[0]["file_id"],
                "rna_file_name": rna_pick.iloc[0]["file_name"],
                "cnv_file_id": cnv_pick.iloc[0]["file_id"],
                "cnv_file_name": cnv_pick.iloc[0]["file_name"],
                "wsi_file_id": wsi_pick.iloc[0]["file_id"] if not wsi_pick.empty else None,
                "wsi_file_name": wsi_pick.iloc[0]["file_name"] if not wsi_pick.empty else None,
                "os_time": s["os_time"],
                "os_event": s["os_event"],
            }
        )

    complete_df = pd.DataFrame(selected_rows)
    return {"rna": rna_df, "cnv": cnv_df, "wsi": wsi_df, "surv": surv_df, "complete": complete_df}


def write_manifest(file_ids: List[str], out_path: str):
    if not file_ids:
        return
    url = f"{GDC_API}/manifest"
    r = requests.post(url, data={"ids": ",".join(file_ids)})
    r.raise_for_status()
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(r.text)


def main():
    ap = argparse.ArgumentParser(description="TCGA 多模态（RNA/CNV/WSI/生存）筛选器")
    sub = ap.add_subparsers(dest="cmd", required=True)

    sp1 = sub.add_parser("list-projects", help="列出 TCGA 可选项目")
    sp1.add_argument("--filter", help="按关键字过滤 project_id/name", default=None)

    sp2 = sub.add_parser("select", help="按项目筛选完整多模态病人并导出下载清单")
    sp2.add_argument("--project", required=True, help="TCGA 项目，如 TCGA-BRCA")
    sp2.add_argument("--out", required=True, help="输出目录")
    sp2.add_argument("--sample-type", default="Primary Tumor", help="样本类型，默认 Primary Tumor")
    sp2.add_argument("--rna-workflow", default="HTSeq - Counts", help="RNA 工作流程，默认 HTSeq - Counts")
    sp2.add_argument(
        "--cnv-type", default="Copy Number Segment", help="CNV 数据类型，默认 Copy Number Segment"
    )
    sp2.add_argument(
        "--allow-case-level",
        action="store_true",
        help="允许仅病人级匹配（WSI 不强制匹配同一 root 样本）",
    )

    args = ap.parse_args()

    if args.cmd == "list-projects":
        df = list_tcga_projects()
        if args.filter:
            mask = df.apply(lambda r: args.filter.lower() in str(r.values).lower(), axis=1)
            df = df[mask]
        if df.empty:
            print("未找到项目。")
            return
        print(df.to_string(index=False))
        return

    if args.cmd == "select":
        os.makedirs(args.out, exist_ok=True)
        print(f"[Info] 开始筛选 {args.project} ...")
        results = select_complete_cases(
            args.project,
            rna_workflow=args.rna_workflow,
            cnv_data_type=args.cnv_type,
            sample_type=args.sample_type,
            require_sample_match=(not args.allow_case_level),
        )

        complete_df = results["complete"]
        if complete_df.empty:
            print("[Warn] 未找到满足条件的完整多模态病人。可尝试加 --allow-case-level 放宽匹配。")
            results["rna"].to_csv(os.path.join(args.out, "rna_all.csv"), index=False)
            results["cnv"].to_csv(os.path.join(args.out, "cnv_all.csv"), index=False)
            results["wsi"].to_csv(os.path.join(args.out, "wsi_all.csv"), index=False)
            results["surv"].to_csv(os.path.join(args.out, "survival_all.csv"), index=False)
            return

        complete_csv = os.path.join(args.out, f"{args.project}_complete_cases.csv")
        complete_df.to_csv(complete_csv, index=False)
        print(f"[Info] 完整多模态病人数：{len(complete_df)}")
        print(f"[Info] 对齐表已保存：{complete_csv}")

        rna_ids = sorted(complete_df["rna_file_id"].dropna().unique().tolist())
        cnv_ids = sorted(complete_df["cnv_file_id"].dropna().unique().tolist())
        wsi_ids = sorted(complete_df["wsi_file_id"].dropna().unique().tolist())

        if rna_ids:
            p = os.path.join(args.out, f"{args.project}_manifest_rna.tsv")
            write_manifest(rna_ids, p)
            print(f"[Info] RNA manifest: {p} ({len(rna_ids)} files)")
        if cnv_ids:
            p = os.path.join(args.out, f"{args.project}_manifest_cnv.tsv")
            write_manifest(cnv_ids, p)
            print(f"[Info] CNV manifest: {p} ({len(cnv_ids)} files)")
        if wsi_ids:
            p = os.path.join(args.out, f"{args.project}_manifest_wsi.tsv")
            write_manifest(wsi_ids, p)
            print(f"[Info] WSI manifest: {p} ({len(wsi_ids)} files)")

        results["rna"].to_csv(os.path.join(args.out, "rna_all.csv"), index=False)
        results["cnv"].to_csv(os.path.join(args.out, "cnv_all.csv"), index=False)
        results["wsi"].to_csv(os.path.join(args.out, "wsi_all.csv"), index=False)
        results["surv"].to_csv(os.path.join(args.out, "survival_all.csv"), index=False)
        print("[Done] 完成。")


if __name__ == "__main__":
    main()

