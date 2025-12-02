# Ageing is not just ageing & Rising disease prevalence — Data Analysis Pipeline

## DOI

### 🔵 Software pipeline (Concept DOI)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17602427.svg)](https://doi.org/10.5281/zenodo.17602427)

### 🟣 Associated manuscripts

**Ageing is not just ageing**  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17662596.svg)](https://doi.org/10.5281/zenodo.17662596)

**Rising disease prevalence signals epigenetic degeneration in humans**  
[![DOI](https://img.shields.io/badge/Preprints-10.20944/preprints202508.2157.v2-blue.svg)](https://doi.org/10.20944/preprints202508.2157.v2)

---

This repository contains the full data-analysis pipeline used in:

1. **Ageing is not just ageing**  
2. **Rising disease prevalence signals epigenetic degeneration in humans**

Both manuscripts rely on the same underlying epidemiological reconstruction of  
**AMAV (Accumulated Mean Annual Variation)** and the associated  
**fold-increase metrics**.

The pipeline reconstructs `AMAV_DATA.xlsx` directly from the consolidated  
prevalence workbooks:

- `data/Supplemental_Table_1.xlsx`
- `data/Supplementary_Table_1.xlsx`

Either filename is accepted automatically.

---

## Repository structure

```text
.
├── README.md
├── LICENSE
├── LICENSE-Data
├── LICENSE-Docs
├── CITATION.cff
├── environment.yml
├── analysis/
│   └── Ageing_is_not_Just-ageing_AMAV_Data_Process_Pipeline.qmd
├── scripts/
│   ├── build_amav_from_supplemental_Table_1.py
│   ├── Figure_1_RAW.R
│   ├── Figure_1_LOSS.R
│   ├── Figure_2_3_from_AMAV_DATA.R
│   ├── FIGURE_4_PANELS_A_and_B.R
│   ├── FIGURE_4_PANEL_C.R
│   ├── FIGURE_S1B.R
│   ├── make_supplemental_figures.py
│   └── remove_monotonic_prefix_FOLD_RELATIVE.R
├── data/
│   ├── Supplemental_Table_1.xlsx
│   ├── Supplementary_Table_1.xlsx
│   ├── Spore_Survival_Loss_F1.csv
│   ├── Spore_Survival_Loss.csv
│   ├── Spore_Survival_RAW_Values_F1.csv
│   └── Spore_Survival_RAW_Values.csv
├── output/
└── .github/workflows/
    └── render.yml
```

> Note: some filenames may differ slightly as the repository evolves,  
> but the overall structure and logic remain the same.

---

## What the pipeline does (summary)

- Reads a single-table prevalence workbook (`Supplemental_Table_1.xlsx` or `Supplementary_Table_1.xlsx`).  
- Each row is one prevalence trend (study). A `Citation` column, if present, is ignored.  
- For each phenotype (disease):
  - Builds piecewise linear slopes between observed points.
  - Computes yearly **MAV** (mean slope across studies).
  - Integrates MAV across the final contiguous block to obtain **AMAV**.
  - If AMAV dips below zero in that block, builds **AMAV-POS** (positive-limb accumulation with configurable handling of negative MAV).

- Aggregates outputs across all phenotypes into `AMAV_DATA.xlsx`, containing:

  - **`AMAV`** – per-disease AMAV/AMAV-POS by calendar year.  
  - **`FOLD_YEARLY`** – yearly fold-increase (normalised by the first positive AMAV/AMAV-POS value).  
  - **`FOLD_RELATIVE`** – the same fold-increase series re-indexed 0…n per disease.  
  - **`LOG_used_column`** – number of trends parsed per phenotype (disease name and number of contributing studies).  

These sheets are then used by the R scripts to generate all **main figures (1–4)** and **supplemental figures** for the two manuscripts.

---

## Requirements

- **Conda/Miniforge** (or Mamba)  
- **Python ≥ 3.12** with:
  - `pandas`, `numpy`, `openpyxl`, `XlsxWriter`, `matplotlib`, `statsmodels`  
- **Quarto ≥ 1.4**  
- **R** with packages:
  - `readxl`, `ggplot2`, `reshape2`, `scales`, `writexl`, `dplyr`, `tidyr`  
- For **PDF** rendering: a TeX distribution with **LuaLaTeX** (e.g. MacTeX/TeX Live)  

All core Python packages are listed in `environment.yml`.

---

## Setup

Create and activate the Conda environment:

```bash
conda env create -f environment.yml
conda activate ageing-amav
```

Install the required R packages in your R installation (outside Conda, if you prefer):

```r
install.packages(c(
  "readxl", "ggplot2", "reshape2", "scales",
  "writexl", "dplyr", "tidyr", "forcats", "patchwork"
))
```

---

## Rendering the Quarto document

From the repository root:

```bash
conda activate ageing-amav
quarto render analysis/Ageing_is_not_Just-ageing_AMAV_Data_Process_Pipeline.qmd
```

This will:

1. Rebuild `output/AMAV_DATA.xlsx` from the appropriate `Supplemental(ary)_Table_1.xlsx`.  
2. Run all R scripts in `scripts/` to generate:
   - Figure 1 (raw prevalence and loss of studies)  
   - Figures 2–3 (AMAV-derived disease trajectories)  
   - Figure 4 (mental vs physical fold-increase dynamics)  
   - Supplemental figures (e.g. S1B and others)  
3. Render HTML and, if enabled, PDF for the Quarto document.

---

## Running the Python script directly

You can also run the Python AMAV builder without Quarto:

```bash
conda activate ageing-amav
python scripts/build_amav_from_supplemental_Table_1.py
```

or explicitly specifying input and output:

```bash
python scripts/build_amav_from_supplemental_Table_1.py   data/Supplementary_Table_1.xlsx   output/AMAV_DATA.xlsx
```

---

## Reproducibility notes

- Year columns are auto-detected (e.g. 1940…2023).  
- Numeric parsing is robust to percentages, comma decimals, and thousand separators.  
- Missing or irregular values are handled gracefully.  
- Results depend only on the content of the `data/` Excel workbooks and CSV files.  

For large `.xlsx` files, Git LFS is recommended:

```bash
git lfs install
git lfs track "*.xlsx"
git add .gitattributes
```

---

## Licence

- **Code:** see `LICENSE`.  
- **Data:** see `LICENSE-Data`.  
- **Docs/Figures:** see `LICENSE-Docs`.  

Please check individual files or subdirectories for any additional restrictions.

---

## Citation

### 1. Software pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17602427.svg)](https://doi.org/10.5281/zenodo.17602427)

> Marsellach, X. (2025).  
> *Ageing is not just ageing & Rising disease prevalence — Data Analysis Pipeline* (v1.1.0) [Software].  
> Zenodo. https://doi.org/10.5281/zenodo.17602427

### 2. Associated manuscripts

> Marsellach, X. (2025). *Ageing is not just ageing*.  
> Zenodo. https://doi.org/10.5281/zenodo.17662596 

> Marsellach, X. (2025). *Rising disease prevalence signals epigenetic degeneration in humans*.  
> Preprints. https://doi.org/10.20944/preprints202508.2157.v2  

Machine-readable citation metadata are provided in `CITATION.cff`.

---

## Contact

For questions or issues, please open a GitHub issue or contact **Xavi Marsellach**.
