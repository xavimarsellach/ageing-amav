# Ageing-amav – Technical README

This repository contains the full AMAV-based epidemiological pipeline used in:

1. **Ageing is not just ageing**  
2. **Rising disease prevalence signals epigenetic degeneration in humans**

Both manuscripts rely on the same reconstruction of:

- **AMAV (Accumulated Mean Annual Variation)**
- Derived **fold-increase metrics** (`FOLD_YEARLY`, `FOLD_RELATIVE`)

The pipeline rebuilds `output/AMAV_DATA.xlsx` from a single prevalence workbook:

- `data/Supplemental_Table_1.xlsx` **or**
- `data/Supplementary_Table_1.xlsx`

The Python script auto-detects which one is present.

---

## Quick start

```bash
# 1. Clone and enter the repo
git clone https://github.com/xavimarsellach/ageing-amav.git
cd ageing-amav

# 2. Create and activate Conda environment
conda env create -f environment.yml
conda activate ageing-amav

# 3. Install required R packages in your R installation
R -q -e 'install.packages(c(
  "readxl","ggplot2","reshape2","scales",
  "writexl","dplyr","tidyr"
), repos = "https://cloud.r-project.org")'
```

---

## Main commands (via Makefile)

After activating the Conda environment:

```bash
# Rebuild AMAV_DATA.xlsx, figures and Quarto document
make all

# Only rebuild AMAV_DATA.xlsx
make amav

# Only regenerate figures
make figures

# Only render the Quarto document
make render
```

---

## How the pipeline works

### 1. Build AMAV and fold metrics (Python)

The script:

```text
scripts/build_amav_from_supplemental_Table_1.py
```

- Parses all studies for each disease  
- Calculates piecewise slopes  
- Computes MAV per year  
- Integrates MAV → **AMAV**  
- Normalises AMAV → **FOLD_YEARLY**  
- Re-indexes the series → **FOLD_RELATIVE**  
- Writes all final sheets into `output/AMAV_DATA.xlsx`

### 2. Generate figures (R)

Scripts in `scripts/` produce:

- Figure 1 RAW trends  
- Figure 1 LOSS-of-studies  
- Figures 2–3 (trajectories from AMAV_DATA.xlsx)  
- Figure 4 (mental vs physical categories, fold dynamics)  
- Supplemental figures (e.g., S1B)

### 3. Render Quarto document

The full pipeline (AMAV + figures) is integrated in:

```text
analysis/Ageing_is_not_Just-ageing_AMAV_Data_Process_Pipeline.qmd
```

```bash
quarto render analysis/Ageing_is_not_Just-ageing_AMAV_Data_Process_Pipeline.qmd
```

---

## Dependency summary

### Python (via environment.yml)

- pandas  
- numpy  
- openpyxl  
- XlsxWriter  
- matplotlib  
- statsmodels  

### R

- readxl  
- ggplot2  
- reshape2  
- scales  
- writexl  
- dplyr  
- tidyr  

### System

- Conda/Miniforge  
- Quarto ≥ 1.4  
- TeX Live / LuaLaTeX (if PDF rendering)

---

## Notes

- `Supplemental_Table_1.xlsx` and `Supplementary_Table_1.xlsx` are interchangeable.  
- Automatic year detection (e.g. 1940–2023).  
- Robust parsing of comma decimals, percentages and missing values.  
- Large `.xlsx` files should use Git LFS:

```bash
git lfs install
git lfs track "*.xlsx"
git add .gitattributes
```

---

## Contact

For issues or questions, please open a GitHub issue or contact **Xavi Marsellach**.
