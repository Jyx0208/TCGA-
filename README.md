TCGA 多模态数据筛选器（RNA/CNV/WSI/生存）

概览
- 这是一个用于从 GDC/TCGA 中筛选“多模态完整病人”的轻量工具，支持跨模态对齐 RNA、CNV、病理全视野图像（WSI）及生存信息，并生成可用于 gdc-client 的下载清单（manifest）。
- 按癌种（如 `TCGA-LIHC`）筛选，默认只保留“Primary Tumor”样本，并在样本级（root sample，TCGA 条码前 15 位）跨模态严格匹配。
- 输出包含：对齐表（病例、root 样本、各模态 file_id、OS 时间/事件）以及各模态 manifest 下载清单。

特性
- 多癌种选择：列出并筛选任意 TCGA 项目（如 BRCA、LIHC 等）。
- 多模态对齐：RNA（HTSeq - Counts，可切换）、CNV（Copy Number Segment）、WSI（Slide Image）。
- 生存计算：基于 GDC `diagnoses` 的 `days_to_death` 与 `days_to_last_follow_up` 计算 OS 时间/事件。
- 一键导出：生成 `*_manifest_*.tsv` 供 `gdc-client` 批量下载。

安装
- 先安装 Python（>=3.8）。
- 安装依赖：
  - Windows: `py -m pip install -r requirements.txt`
  - macOS/Linux: `python3 -m pip install -r requirements.txt`

快速开始
- 列出 TCGA 项目（支持关键字过滤）
  - `python tcga_mm_selector.py list-projects`
  - `python tcga_mm_selector.py list-projects --filter LIHC`

- 按癌种筛选完整多模态病人（严格样本级匹配）
  - `python tcga_mm_selector.py select --project TCGA-LIHC --out ./out_lihc`

- 若严格样本级匹配过于严格，可放宽到病人级（不强制 WSI 与 RNA/CNV 在相同 root 样本）
  - `python tcga_mm_selector.py select --project TCGA-LIHC --out ./out_lihc --allow-case-level`

输出说明
- `out_lihc/TCGA-LIHC_complete_cases.csv`：多模态对齐表（病例、root 样本、各模态 file_id、OS 时间/事件）。
- `out_lihc/TCGA-LIHC_manifest_rna.tsv`：RNA 下载清单（manifest）。
- `out_lihc/TCGA-LIHC_manifest_cnv.tsv`：CNV 下载清单（manifest）。
- `out_lihc/TCGA-LIHC_manifest_wsi.tsv`：WSI 下载清单（manifest）。
- `out_lihc/{rna_all,cnv_all,wsi_all,survival_all}.csv`：原始汇总，便于复核与排查。

使用 gdc-client 下载
1) 到 GDC 官方页面下载 `gdc-client` 二进制并将其加入 PATH。
   - https://gdc.cancer.gov/access-data/gdc-data-transfer-tool
2) 使用 manifest 批量下载：
   - RNA：`gdc-client download -m out_lihc/TCGA-LIHC_manifest_rna.tsv -d data_lihc/rna`
   - CNV：`gdc-client download -m out_lihc/TCGA-LIHC_manifest_cnv.tsv -d data_lihc/cnv`
   - WSI：`gdc-client download -m out_lihc/TCGA-LIHC_manifest_wsi.tsv -d data_lihc/wsi`

命令行参数
- `list-projects`
  - `--filter`：按关键字过滤项目 ID/名称。
- `select`
  - `--project`：TCGA 项目（如 `TCGA-LIHC`）。
  - `--out`：输出目录。
  - `--sample-type`：样本类型（默认 `Primary Tumor`）。
  - `--rna-workflow`：RNA 工作流（默认 `HTSeq - Counts`，可设为 `HTSeq - FPKM` 等）。
  - `--cnv-type`：CNV 数据类型（默认 `Copy Number Segment`）。
  - `--allow-case-level`：放宽为病人级匹配（WSI 不强制同一 root 样本）。

生存定义与匹配逻辑
- Root 样本定义：TCGA 条形码前 15 位（例如 `TCGA-XX-YYYY-01`）。
- 匹配策略：默认严格在 root 样本级对齐 RNA、CNV、WSI；如需放宽，使用 `--allow-case-level`。
- OS 计算：
  - 若 `vital_status == Dead` 且有 `days_to_death`，则 `os_event=1`，`os_time=days_to_death`；
  - 否则使用 `days_to_last_follow_up`，`os_event=0`；
  - 若两者皆缺失，则该病例不计入完整生存集。

实践建议
- WSI 体积巨大（单张百 MB 至数 GB），建议先按 `*_complete_cases.csv` 抽样验证流程。
- 国内网络访问 GDC 可能较慢，建议配置稳定网络或重试策略。
- 若输出为空：
  - 放宽到病人级匹配：添加 `--allow-case-level`；
  - 检查 `--sample-type` 是否与该癌种的可用样本一致；
  - 检查 `--rna-workflow` 与 `--cnv-type` 是否存在对应数据。

项目结构
- `tcga_mm_selector.py`：主脚本与 CLI。
- `requirements.txt`：Python 依赖。
- `README.md`：使用说明。

开源许可
- 本仓库未附带许可证文件。若你需要将其正式开源发布，请添加 `LICENSE` 文件（例如 MIT 或 Apache-2.0），并在 `README.md` 中注明。

致谢与引用
- 数据来源：GDC（Genomic Data Commons）/TCGA。使用数据请遵守 GDC/TCGA 的数据使用政策和条款。

