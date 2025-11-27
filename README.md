# Ageing is not just ageing & Rising disease prevalence — Data Analysis Pipeline

This repository contains the full data-analysis pipeline used in:

1. **Ageing is not just ageing**  
2. **Rising disease prevalence signals epigenetic degeneration in humans**

Both manuscripts rely on the same underlying epidemiological reconstruction of  
**AMAV (Accumulated Mean Annual Variation)** and the associated  
**fold-increase metrics**.

The pipeline reconstructs `AMAV_DATA.xlsx` directly from the consolidated  
prevalence workbooks:

- `data/Supplemental_Table_1.xlsx` (used in *Ageing is not just ageing*)
- `data/Supplementary_Table_1.xlsx` (used in *Rising disease prevalence…*)

Either filename is accepted automatically.

---

## Repository structure

```
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
│   ├── Figure_1_LOSS.R
│   ├── Figure_1_RAW.R
│   ├── Figure_2_3_from_AMAV_DATA.R
│   ├── FIGURE_4_PANEL_C.R
│   ├── FIGURE_4_PANELS_A_and_B.R
│   ├── FIGURE_S1B.R
│   ├── make_supplemental_figres.py
│   └── remove_monotonic_prefix_FOLD_RELATIVE.R
├── data/
│   ├── Supplemental_Table_1.xlsx
│   ├── Supplementary_Table_1.xlsx
│   ├── Spore_Survival_Loss_F1.csv
│   ├── Spore_Survival_Loss.csv
│   ├── Spore_Survival_RAW_Values_F1.csv
│   └── Spore_Survival_RAW_Values.csv
├── output/
│   ├── AMAV_DATA.xlsx
│   ├── FIGURE_4A.png
│   ├── FIGURE_4A_B.png
│   ├── FIGURE_4B.png
│   ├── FIGURE_4C.png
│   └── GEOMETRIC_MEANS_MENTAL_vs_PHYSICAL_ANCOVA_thresholds_k0_to_k10.xlsx
└── .github/workflows/
    └── render.yml
```

---

## What the pipeline does (summary)

- Reads a single-table prevalence workbook (`Supplemental_Table_1.xlsx` or `Supplementary_Table_1.xlsx`).
- Each row is one prevalence trend (study). A `Citation` column is ignored.
- For each phenotype:
  - Builds piecewise linear slopes between observed points.
  - Computes yearly **MAV** (mean slope across studies).
  - Integrates MAV across the final contiguous block to obtain **AMAV**.
  - If AMAV dips below zero, builds **AMAV-POS**.

- Aggregates outputs across all phenotypes into `AMAV_DATA.xlsx`, containing:

  - **`AMAV`**  
    Per-disease AMAV/AMAV-POS by calendar year.

  - **`FOLD_YEARLY`**  
    Yearly fold-increase (normalised by the first positive AMAV value).

  - **`FOLD_RELATIVE`**  
    Same fold-increase series re-indexed 0…n per disease.

  - **`LOG_used_column`**  
    Number of trends parsed per phenotype.

In addition, the `data/Spore_Survival_*.csv` files contain the raw spore-survival
data used to generate **Figure 1** (raw prevalence and loss of studies).

---

## Requirements

- Conda/Miniforge
- Python ≥ 3.12  
  (`pandas`, `numpy`, `openpyxl`, `XlsxWriter`, `matplotlib`)
- Quarto ≥ 1.4
- LuaLaTeX (if rendering PDF)

---

## Setup

```bash
conda env create -f environment.yml
conda activate ageing-amav
```

---

## Render the Quarto document

```bash
quarto render analysis/Ageing_is_not_Just-ageing_AMAV_Data_Process_Pipeline.qmd
```

This:

- Rebuilds `output/AMAV_DATA.xlsx`
- Generates all main and supplemental figures
- Renders HTML/PDF depending on options

---

## Running the script directly

```bash
python scripts/build_amav_from_supplemental_Table_1.py
```

or explicitly:

```bash
python scripts/build_amav_from_supplemental_Table_1.py   data/Supplementary_Table_1.xlsx   output/AMAV_DATA.xlsx
```

---

## Reproducibility notes

- Year columns auto-detected.
- Numeric parsing robust to %, comma decimals, thousand separators.
- Missing or irregular values handled gracefully.
- The Quarto document orchestrates the full pipeline:
  1. Rebuilds `AMAV_DATA.xlsx`.
  2. Runs all figure scripts in `scripts/`.
  3. Writes all derived files into `output/`.

Git LFS recommended for large `.xlsx` files:

```bash
git lfs install
git lfs track "*.xlsx"
git add .gitattributes
```

---

## Licence

Code → `LICENSE`  
Data & Docs → `LICENSE-Data`, `LICENSE-Docs`

---

## Citation

### Software pipeline

> Marsellach, X. (2025).  
> *Ageing is not just ageing & Rising disease prevalence — Data Analysis Pipeline* (v1.1.0).  
> Zenodo. https://doi.org/10.5281/zenodo.17602428

### Manuscripts

> Marsellach, X. (2025). *Ageing is not just ageing*.  
> https://doi.org/10.5281/zenodo.15596247

> Marsellach, X. (2025). *Rising disease prevalence signals epigenetic degeneration in humans*.  
> https://doi.org/10.20944/preprints202508.2157.v2

---

## Contact

For questions or issues, please open a GitHub issue or contact **Xavi Marsellach**.
