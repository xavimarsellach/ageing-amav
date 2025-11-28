# ----------------------------
# Makefile for ageing-amav repo
# ----------------------------

# Activate environment automatically if possible
SHELL := /bin/bash
CONDA_ENV = ageing-amav

all: amav figures render

# 1) Rebuild AMAV_DATA.xlsx
amav:
	conda run -n $(CONDA_ENV) python scripts/build_amav_from_supplemental_Table_1.py

# 2) Regenerate all figures
figures:
	conda run -n $(CONDA_ENV) Rscript scripts/Figure_1_RAW.R
	conda run -n $(CONDA_ENV) Rscript scripts/Figure_1_LOSS.R
	conda run -n $(CONDA_ENV) Rscript scripts/Figure_2_3_from_AMAV_DATA.R
	conda run -n $(CONDA_ENV) Rscript scripts/FIGURE_4_PANELS_A_and_B.R
	conda run -n $(CONDA_ENV) Rscript scripts/FIGURE_4_PANEL_C.R
	conda run -n $(CONDA_ENV) python scripts/make_supplemental_figures.py

# 3) Render Quarto document
render:
	conda run -n $(CONDA_ENV) quarto render analysis/Ageing_is_not_Just-ageing_AMAV_Data_Process_Pipeline.qmd

# Utility: clean output
clean:
	rm -f output/*.xlsx
	rm -f output/*.png
	rm -f analysis/*.html
	rm -f analysis/*.pdf