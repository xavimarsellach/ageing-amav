#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Rebuild AMAV_DATA.xlsx from a single-table workbook (Supplemental(ary)_Table_1.xlsx).

Input expectation (neutral description):
- One Excel workbook containing a single wide table with:
    Phenotype | Citation | 1941 | 1942 | ... | 2023
- Each row is one prevalence trend (study) for a given phenotype/disease.
- A "Citation" column may be present and is ignored for calculations.

What the script computes:
- For each phenotype, piecewise-linear yearly slopes per study between observed points.
- Yearly MAV (mean of study slopes).
- AMAV: cumulative MAV restricted to the last contiguous block of defined MAV.
- AMAV-POS: if AMAV in that block dips below zero, a positive-limb accumulation
  starting after the AMAV valley (MODE='prev'), optionally including negative MAV.
- Aggregated outputs (written to AMAV_DATA.xlsx):
    • AMAV           : one column per phenotype (indexed by Year) using AMAV-POS
                       if it exists, otherwise AMAV.
    • FOLD_YEARLY    : AMAV normalised by the first positive value (by Year).
    • FOLD_RELATIVE  : Same normalisation but re-indexed 0..n per phenotype
                       (relative time).
    • LOG_used_column: Phenotype name and number of trends parsed.

Default paths (relative to current working directory):
- Input  : data/Supplemental_Table_1.xlsx or data/Supplementary_Table_1.xlsx
          (auto-discovery also checks analysis/data and ./)
- Output : output/AMAV_DATA.xlsx             (directory is created if missing)

CLI:
    python scripts/build_amav_from_supplemental_Table_1.py [input.xlsx] [output.xlsx]

Programmatic (e.g., from R with reticulate):
    from scripts.build_amav_from_supplemental_Table_1 import build_amav
    out_path = build_amav()  # or build_amav("path/to/input.xlsx", "output/AMAV_DATA.xlsx")
"""

from __future__ import annotations

from pathlib import Path
import math
import re
import sys
import warnings
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ---- Behaviour switches ----
MODE = "prev"        # accumulate MAV from the previous year
INCLUDE_NEG = True   # include negative MAV in AMAV-POS accumulation
EPS = 1e-12

# ---- Default locations (relative to working directory) ----
DEFAULT_INPUT  = Path("data/Supplemental_Table_1.xlsx")
DEFAULT_OUTPUT = Path("output/AMAV_DATA.xlsx")


# ---------------- Utilities ----------------
def to_number(series: pd.Series) -> pd.Series:
    """Robust numeric conversion: handles %, comma decimals, and thousand separators."""
    def _one(x):
        if pd.isna(x):
            return np.nan
        if isinstance(x, (int, float, np.number)):
            return float(x)
        s = str(x).strip().replace(" ", "")
        is_pct = s.endswith("%")
        s = s.replace("%", "").replace(",", ".")
        # remove thousands separators like 1.234 or 1'234
        s = re.sub(r"(?<=\d)[\.,'’](?=\d{3}\b)", "", s)
        try:
            v = float(s)
        except ValueError:
            return np.nan
        return v / 100.0 if is_pct else v
    return series.apply(_one)


def runs_true(mask: np.ndarray) -> List[Tuple[int, int]]:
    """Return [(start, end)] index intervals where mask is True and contiguous."""
    idx = np.where(mask)[0]
    if idx.size == 0:
        return []
    runs, start, prev = [], idx[0], idx[0]
    for k in idx[1:]:
        if k == prev + 1:
            prev = k
        else:
            runs.append((start, prev))
            start = prev = k
    runs.append((start, prev))
    return runs


# -------- AMAV / AMAV-POS core --------
def build_perstudy_linearslope(raw_df: pd.DataFrame) -> pd.DataFrame:
    """
    raw_df columns: Year | T1 | T2 | ... | Tk
    For each trend, build a piecewise-constant slope between successive observed points.
    """
    df = raw_df.sort_values("Year").reset_index(drop=True)
    years = df["Year"].to_numpy()
    out = pd.DataFrame({"Year": years})

    for col in df.columns:
        if col == "Year":
            continue
        s = pd.to_numeric(df[col], errors="coerce")
        mask = s.notna().to_numpy()
        obs_y = years[mask]
        obs_v = s[mask].to_numpy(float)

        arr = np.full(len(years), np.nan, float)
        for i in range(len(obs_y) - 1):
            y0, y1 = obs_y[i], obs_y[i + 1]
            v0, v1 = obs_v[i],  obs_v[i + 1]
            if y1 > y0:
                slope = (v1 - v0) / (y1 - y0)
                seg = (years > y0) & (years <= y1)
                arr[seg] = slope
        out[col] = arr

    return out


def build_summary_mav_amav(slope_df: pd.DataFrame) -> pd.DataFrame:
    """
    - MAV: mean across study slopes per year.
    - AMAV: cumulative MAV within the last contiguous block where MAV is defined.
    - AMAV-POS: if AMAV dips below zero in that block, accumulate MAV from the year
      after the valley (MODE='prev'), including negative MAV if INCLUDE_NEG=True.
    """
    df = slope_df.sort_values("Year").reset_index(drop=True)
    years = df["Year"].to_numpy()

    study_cols = [c for c in df.columns if c != "Year"]
    M = df[study_cols].to_numpy(float)

    datapoints = np.sum(~np.isnan(M), axis=1)
    with np.errstate(all="ignore"):
        mav = np.nanmean(M, axis=1)

    # AMAV within the last contiguous MAV block
    amav = np.full(len(years), np.nan, float)
    finite_mav = ~np.isnan(mav)
    mav_runs = runs_true(finite_mav)
    last_s = last_e = None
    if mav_runs:
        last_s, last_e = mav_runs[-1]
        acc = 0.0
        # i from last_s+1 .. last_e+1 inclusive (AMAV[i] integrates mav[i-1])
        for i in range(last_s + 1, last_e + 2):
            prev = i - 1
            if not math.isnan(mav[prev]):
                acc += mav[prev]
            amav[i] = acc

    summary = pd.DataFrame({
        "Year": years,
        "DataPoints": datapoints,
        "MAV": mav,
        "AMAV": amav,
    })

    # AMAV-POS only if AMAV dips below zero in the final block
    if last_s is None:
        return summary

    lo, hi = last_s + 1, last_e + 1  # indices where AMAV is defined
    window_am = amav[lo:hi + 1]
    has_negative = np.isfinite(window_am).any() and (window_am[np.isfinite(window_am)] < -EPS).any()
    if not has_negative:
        return summary

    valley_rel = int(np.nanargmin(window_am))
    valley = lo + valley_rel

    amav_pos = np.full(len(amav), np.nan, float)
    run = 0.0

    # Start at valley+1 for MODE='prev'
    start_i = max(valley + 1, lo) if MODE == "prev" else max(valley, lo)
    for i in range(start_i, hi + 1):
        idx_mav = (i - 1) if MODE == "prev" else i
        v = mav[idx_mav]
        if np.isfinite(v):
            contrib = v if INCLUDE_NEG else max(v, 0.0)
            run += contrib
            amav_pos[i] = run

    insert_at = list(summary.columns).index("AMAV") + 1
    summary.insert(insert_at, "AMAV-POS", amav_pos)
    return summary


def first_positive_baseline(s: pd.Series) -> float:
    """First non-zero, non-NA value."""
    for v in s:
        if pd.notna(v) and v != 0:
            return float(v)
    return np.nan


# -------- Read single table from Supplemental(ary)_Table_1 --------
def find_table_sheet(xlsx: Path) -> str:
    """Find the sheet whose (0,0) cell is 'Phenotype' (fallback: first sheet)."""
    xls = pd.ExcelFile(xlsx)
    for sh in xls.sheet_names:
        test = pd.read_excel(xlsx, sheet_name=sh, nrows=1, header=None)
        if not test.empty and str(test.iat[0, 0]).strip().lower() == "phenotype":
            return sh
    return xls.sheet_names[0]


def read_single_table(input_xlsx: Path) -> tuple[pd.DataFrame, List[int]]:
    """
    Returns:
    - df: columns = Phenotype | <year_int_1> | ... | <year_int_k>
    - year_order: sorted list of integer years

    Drops the 'Citation' column. Removes a header-like duplicate row if needed.
    """
    sheet = find_table_sheet(input_xlsx)
    df = pd.read_excel(input_xlsx, sheet_name=sheet, header=0)
    if df.empty:
        raise SystemExit("Input sheet appears to be empty.")

    if "Phenotype" not in df.columns:
        raise SystemExit("Cannot find 'Phenotype' column in header.")

    # Drop Citation column if present
    drop_cit = [c for c in df.columns if str(c).strip().lower() == "citation"]
    if drop_cit:
        df = df.drop(columns=drop_cit)

    # Identify year columns (numeric or string)
    year_cols_orig = []
    for c in df.columns:
        if c == "Phenotype":
            continue
        if isinstance(c, (int, float)) and 1800 <= int(c) <= 2100:
            year_cols_orig.append(c)
        elif isinstance(c, str) and re.fullmatch(r"\s*\d{4}\s*", c or ""):
            year_cols_orig.append(c)

    if not year_cols_orig:
        raise SystemExit("No year-like columns found in the table.")

    # Normalize year columns to integer names
    year_map = {c: int(str(c).strip()) for c in year_cols_orig}
    df = df.rename(columns=year_map)
    year_cols = sorted(year_map.values())

    # Keep only Phenotype + years
    keep = ["Phenotype"] + year_cols
    df = df[keep].copy()
    df["Phenotype"] = df["Phenotype"].astype(str).str.strip()

    # Remove bogus header-like row
    df = df[df["Phenotype"].str.strip().str.lower() != "phenotype"]

    # Coerce year values to numeric
    for c in year_cols:
        df[c] = to_number(df[c])

    return df, year_cols


def build_disease_series(df: pd.DataFrame, year_order: List[int]) -> Dict[str, pd.Series]:
    """
    For each phenotype, build a per-disease AMAV/AMAV-POS series,
    preferring AMAV-POS when present, else AMAV.
    """
    series_by: Dict[str, pd.Series] = {}
    years = list(year_order)

    for disease, grp in df.groupby("Phenotype", dropna=True):
        raw = pd.DataFrame({"Year": years})
        # Each row in grp is one trend (T1..Tk)
        for i, (_, row) in enumerate(grp.iterrows(), start=1):
            raw[f"T{i}"] = row[years].values

        # Skip if all trends are NA
        if raw.drop(columns=["Year"]).isna().all().all():
            continue

        slope = build_perstudy_linearslope(raw)
        summary = build_summary_mav_amav(slope)

        chosen_col = "AMAV-POS" if ("AMAV-POS" in summary.columns and summary["AMAV-POS"].notna().any()) else "AMAV"
        ser = summary.set_index("Year")[chosen_col].copy()
        ser.name = disease
        series_by[disease] = ser

    return series_by


# ---------------- I/O helpers ----------------
def resolve_input_output(
    input_xlsx: Path | None,
    output_xlsx: Path | None
) -> tuple[Path, Path]:
    """Resolve paths with sensible defaults and fallbacks."""
    cwd = Path.cwd()

    # Input resolution
    if input_xlsx is None:
        candidates = [
            cwd / "data" / "Supplemental_Table_1.xlsx",
            cwd / "data" / "Supplementary_Table_1.xlsx",
            cwd / "analysis" / "data" / "Supplemental_Table_1.xlsx",
            cwd / "analysis" / "data" / "Supplementary_Table_1.xlsx",
            cwd / "Supplemental_Table_1.xlsx",
            cwd / "Supplementary_Table_1.xlsx",
        ]
        inp = next((p for p in candidates if p.exists()), None)
        if inp is None:
            raise SystemExit(
                "Input file not found. Expected at one of: "
                "data/Supplemental_Table_1.xlsx, data/Supplementary_Table_1.xlsx, "
                "analysis/data/Supplemental_Table_1.xlsx, analysis/data/Supplementary_Table_1.xlsx, "
                "./Supplemental_Table_1.xlsx, ./Supplementary_Table_1.xlsx"
            )
    else:
        inp = Path(input_xlsx)
        if not inp.exists():
            raise SystemExit(f"Input file not found: {inp}")

    # Output resolution
    out = Path(output_xlsx) if output_xlsx is not None else (cwd / DEFAULT_OUTPUT)
    out.parent.mkdir(parents=True, exist_ok=True)

    return inp, out


def write_workbook(amav_df: pd.DataFrame, folds_y: pd.DataFrame,
                   folds_rel: pd.DataFrame, log_df: pd.DataFrame, out_path: Path) -> None:
    """Write the 4 output sheets to an Excel workbook."""
    with pd.ExcelWriter(out_path, engine="xlsxwriter") as xw:
        amav_df.reset_index().rename(columns={"index": "Year"}).to_excel(
            xw, index=False, sheet_name="AMAV"
        )
        folds_y.reset_index().rename(columns={"index": "Year"}).to_excel(
            xw, index=False, sheet_name="FOLD_YEARLY"
        )
        folds_rel.rename_axis("RelYear").reset_index().to_excel(
            xw, index=False, sheet_name="FOLD_RELATIVE"
        )
        log_df.to_excel(xw, index=False, sheet_name="LOG_used_column")


# ---------------- Orchestrator ----------------
def build_amav(input_xlsx: str | Path | None = None,
               output_xlsx: str | Path | None = None) -> Path:
    """
    High-level entry point. Returns the Path of the created Excel workbook.
    """
    inp, out = resolve_input_output(
        Path(input_xlsx) if input_xlsx is not None else None,
        Path(output_xlsx) if output_xlsx is not None else None,
    )

    df, year_order = read_single_table(inp)
    series_by = build_disease_series(df, year_order)
    if not series_by:
        raise SystemExit("No diseases could be parsed into AMAV/AMAV-POS series.")

    # AMAV: union of years across all diseases
    all_years = sorted(set().union(*[set(s.index) for s in series_by.values()]))
    amav_df = pd.DataFrame(index=all_years)
    for disease, s in series_by.items():
        amav_df[disease] = s.reindex(all_years)

    # FOLD_YEARLY: normalised by first positive baseline, indexed by Year
    folds_y = pd.DataFrame(index=amav_df.index)
    for col in amav_df.columns:
        v = amav_df[col]
        b = first_positive_baseline(v.dropna())
        folds_y[col] = v / b if pd.notna(b) else np.nan

    # FOLD_RELATIVE: same but re-indexed 0..n per disease
    rel_dict: Dict[str, pd.Series] = {}
    max_len = 0
    for col in amav_df.columns:
        v = amav_df[col].dropna()
        b = first_positive_baseline(v)
        rel = (v / b) if pd.notna(b) else pd.Series(dtype=float)
        rel.index = range(len(rel))
        rel_dict[col] = rel
        max_len = max(max_len, len(rel))

    folds_rel = pd.DataFrame(index=range(max_len))
    for col, s in rel_dict.items():
        folds_rel[col] = s.reindex(folds_rel.index)

    # LOG sheet: phenotype + number of trends parsed
    log_df = (
        df.groupby("Phenotype", dropna=True)["Phenotype"]
        .count()
        .rename("Trends")
        .reset_index()
        .rename(columns={"Phenotype": "Disease"})
        .sort_values("Disease")
    )

    write_workbook(amav_df, folds_y, folds_rel, log_df, out)
    print(f"✅ OK: created {out}")
    return out


# ---------------- CLI ----------------
def main() -> None:
    # CLI overrides: [input.xlsx] [output.xlsx]
    argv = sys.argv
    in_arg = Path(argv[1]) if len(argv) >= 2 else None
    out_arg = Path(argv[2]) if len(argv) >= 3 else None
    build_amav(in_arg, out_arg)


if __name__ == "__main__":
    main()
