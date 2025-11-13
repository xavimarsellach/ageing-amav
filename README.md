# Ageing is not just ageing — Data Analysis Pipeline

This repository contains a Quarto document that reconstructs `AMAV_DATA.xlsx` directly from the consolidated workbook `data/Supplemental_Table_1.xlsx`.  
It implements the AMAV/AMAV-POS computation and exports the aggregated sheets **`AMAV-P`**, **`AMAV-F-Y`**, **`AMAV-F-R`**, and a **`LOG_used_column`**.

---

## Repository structure

```
.
├── README.md
├── LICENSE                  # Code licence (e.g. MIT)
├── LICENSE-Data             # Data licence (e.g. CC BY 4.0)
├── LICENSE-Docs             # Docs/figures licence (e.g. CC BY 4.0)
├── CITATION.cff
├── environment.yml          # Conda environment (Python)
├── analysis/
│   └── Ageing_pipeline.qmd  # Quarto document (R/knitr orchestrating Python)
├── scripts/
│   └── build_amav_from_supplemental_Table_1.py
├── data/
│   └── Supplemental_Table_1.xlsx
├── output/                  # Rendered outputs (gitignored)
└── .github/workflows/
    └── render.yml           # (Optional) CI to render Quarto
```

---

## What the pipeline does (summary)

- Reads the single-table input `Supplemental_Table_1.xlsx` (one row = one prevalence trend; a `Citation` column, if present, is ignored).
- For each phenotype (disease), builds piecewise linear slopes per study between observed points.
- Computes yearly **MAV** (mean slope across studies) and integrates within the final contiguous block to obtain **AMAV**.
- If AMAV dips below zero in that block, builds **AMAV-POS** (accumulating from the year after the AMAV valley; negatives allowed).
- Aggregates outputs across all phenotypes into `AMAV_DATA.xlsx` with sheets:
  - `AMAV-P` (per-disease AMAV/AMAV-POS by calendar year)
  - `AMAV-F-Y` (normalised by first positive AMAV/AMAV-POS value, by year)
  - `AMAV-F-R` (same normalisation, re-indexed from 0 per disease)
  - `LOG_used_column` (phenotype name and number of trends parsed)

---

## Requirements

- **Conda/Miniforge** (or Mamba)
- **Python 3.12+** with:
  - `pandas`, `numpy`, `openpyxl`, `XlsxWriter`, `matplotlib`
- **Quarto** (≥ 1.4) — only if you want to render the `.qmd` to HTML/PDF
- For **PDF** rendering: a TeX distribution with **LuaLaTeX** (e.g. MacTeX/TeX Live)

All Python packages are pinned in `environment.yml`.

---

## Setup

### 1) Create and activate the environment

```bash
conda env create -f environment.yml
conda activate ageing-amav
```

### 2) (Optional) Render the Quarto document

```bash
quarto render analysis/Ageing_is_not_Just-ageing_AMAV_Data_Process_Pipeline.qmd
```

- Rendering writes the Excel output to `output/AMAV_DATA.xlsx` (created if missing).
- The document also renders `analysis/Ageing_is_not_Just-ageing_AMAV_Data_Process_Pipeline.html` (and, if enabled, PDF).

---

## Script entry points

You can run the Python script directly (bypassing Quarto):

```bash
conda activate ageing-amav
python scripts/build_amav_from_supplemental_Table_1.py        data/Supplemental_Table_1.xlsx        AMAV_DATA.xlsx
```

- First argument: input workbook (default search: `data/Supplemental_Table_1.xlsx`).
- Second argument: output file (default: `AMAV_DATA.xlsx` at the repo root).

---

## Reproducibility notes

- Results depend only on the content of `data/Supplemental_Table_1.xlsx`.
- Year columns are auto-detected (e.g. `1940 … 2023`); a `Citation` column is ignored if present.
- Numeric parsing is robust to percent signs, comma decimals, and thousand separators.

For very large `.xlsx`, consider enabling **Git LFS**:

```bash
git lfs install
git lfs track "*.xlsx"
git add .gitattributes
```

---

## Licence

- **Code:** see `LICENSE` (e.g. MIT).
- **Data & Docs:** see `LICENSE-Data` / `LICENSE-Docs` (e.g. CC BY 4.0).

State any third-party data restrictions here if applicable.

---

## Citation

[![DOI](https://zenodo.org/badge/1095869874.svg)](https://doi.org/10.5281/zenodo.17602427)
Please cite this repository as:

> Marsellach, X. (2025). *Ageing is not just ageing — Data Analysis Pipeline* (v1.0.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.XXXXXX

Machine-readable citation metadata are provided in `CITATION.cff`.
---

## Contact

For questions or issues, please open a GitHub issue or contact **Xavi Marsellach**.
