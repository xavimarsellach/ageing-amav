#| label: make_supplemental_figures.py
import math
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt

# ===========================
# CONFIG
# ===========================

INPUT_XLSX = Path("data/Supplemental_Table_1.xlsx")   # canvia el path si cal
OUT_DIR    = Path("output/Supplemental_Figures_S2-S32")            # carpeta on es guardaran els PNG

# Diccionari per afinar malaltia per malaltia (pots deixar-lo buit)
# Exemple:
# PLOT_OVERRIDES = {
#     "Asthma":  {"xmin": 1980, "xmax": 2018},
#     "Anaemia": {"xmin": 2002},
# }
PLOT_OVERRIDES: Dict[str, Dict[str, float]] = {}


# ===========================
# Utilitats generals
# ===========================

def slugify(name: str) -> str:
    """Safe file name from disease name."""
    return (
        name.replace("/", "_")
            .replace(" ", "_")
            .replace("(", "")
            .replace(")", "")
            .replace("'", "")
    )


def to_number(series: pd.Series) -> pd.Series:
    """Robust numeric conversion (%, commas, thousands separators, etc.)."""
    import re as _re

    def _one(x):
        if pd.isna(x):
            return np.nan
        if isinstance(x, (int, float, np.number)):
            return float(x)
        s = str(x).strip().replace(" ", "")
        is_pct = s.endswith("%")
        s = s.replace("%", "").replace(",", ".")
        # thousands separators: 1.234 / 1'234 / 1 234
        s = _re.sub(r"(?<=\d)[\.,'’](?=\d{3}\b)", "", s)
        try:
            v = float(s)
        except ValueError:
            return np.nan
        return v / 100.0 if is_pct else v

    return series.apply(_one)


def runs_true(mask: np.ndarray):
    """Return contiguous runs where mask is True."""
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


# ===========================
# AMAV / AMAV-POS
# ===========================

def build_perstudy_linearslope(raw_df: pd.DataFrame) -> pd.DataFrame:
    """
    raw_df: Year | T1 | T2 | ... | Tk (prevalence per study).
    Builds piecewise-constant slopes between observed points for each study.
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
    - MAV: mean of slopes per year.
    - AMAV: cumulative MAV within the last contiguous block with MAV defined.
    - AMAV-POS: if AMAV goes below zero in that block, accumulate from the valley.
    """
    df = slope_df.sort_values("Year").reset_index(drop=True)
    years = df["Year"].to_numpy()

    study_cols = [c for c in df.columns if c != "Year"]
    M = df[study_cols].to_numpy(float)

    # Manual mean to avoid warnings on empty rows
    mav = np.full(len(years), np.nan, float)
    for i in range(len(years)):
        row = M[i, :]
        mask = np.isfinite(row)
        if mask.any():
            mav[i] = row[mask].mean()

    datapoints = np.sum(np.isfinite(M), axis=1)

    # AMAV within the last contiguous MAV block
    n = len(mav)
    amav = np.full(n, np.nan, float)
    finite_mav = np.isfinite(mav)
    mav_runs = runs_true(finite_mav)
    last_s = last_e = None
    if mav_runs:
        last_s, last_e = mav_runs[-1]
        acc = 0.0
        hi_idx = min(last_e + 1, n - 1)
        for i in range(last_s + 1, hi_idx + 1):
            prev = i - 1
            if np.isfinite(mav[prev]):
                acc += mav[prev]
            amav[i] = acc

    summary = pd.DataFrame({
        "Year": years,
        "DataPoints": datapoints,
        "MAV": mav,
        "AMAV": amav,
    })

    # AMAV-POS only if there are negative AMAV values in the final block
    if last_s is None:
        return summary

    lo, hi = last_s + 1, last_e + 1
    window_am = amav[lo:hi + 1]
    valid_window = np.isfinite(window_am)
    has_negative = valid_window.any() and (window_am[valid_window] < 0).any()
    if not has_negative:
        return summary

    valley_rel = int(np.nanargmin(window_am))
    valley = lo + valley_rel

    amav_pos = np.full(len(amav), np.nan, float)
    run = 0.0
    start_i = max(valley + 1, lo)
    for i in range(start_i, min(hi, n - 1) + 1):
        idx_mav = i - 1
        v = mav[idx_mav]
        if np.isfinite(v):
            run += v
            amav_pos[i] = run

    insert_at = list(summary.columns).index("AMAV") + 1
    summary.insert(insert_at, "AMAV-POS", amav_pos)
    return summary


# ===========================
# Reading Supplemental_Table_1
# ===========================

def read_single_table(xlsx: Path) -> (pd.DataFrame, List[int]):
    """Read big table and return df + list of year columns."""
    xls = pd.ExcelFile(xlsx)
    sheet = None
    for sh in xls.sheet_names:
        test = pd.read_excel(xlsx, sheet_name=sh, nrows=1, header=None)
        if not test.empty and str(test.iat[0, 0]).strip().lower() == "phenotype":
            sheet = sh
            break
    if sheet is None:
        sheet = xls.sheet_names[0]

    df = pd.read_excel(xlsx, sheet_name=sheet, header=0)
    if "Phenotype" not in df.columns:
        raise SystemExit("Cannot find 'Phenotype' column in input file.")

    # Drop "Citation" if present
    drop_cit = [c for c in df.columns if str(c).strip().lower() == "citation"]
    if drop_cit:
        df = df.drop(columns=drop_cit)

    # Year columns
    year_cols = []
    for c in df.columns:
        if c == "Phenotype":
            continue
        if isinstance(c, (int, float)) and 1800 <= int(c) <= 2100:
            year_cols.append(c)
        elif isinstance(c, str) and c.strip().isdigit() and 1800 <= int(c) <= 2100:
            year_cols.append(c)

    if not year_cols:
        raise SystemExit("No year-like columns found in the table.")

    year_map = {c: int(str(c).strip()) for c in year_cols}
    df = df.rename(columns=year_map)
    year_cols = sorted(year_map.values())

    keep = ["Phenotype"] + year_cols
    df = df[keep].copy()
    df["Phenotype"] = df["Phenotype"].astype(str).str.strip()
    df = df[df["Phenotype"].str.lower() != "phenotype"]

    for c in year_cols:
        df[c] = to_number(df[c])

    return df, year_cols


# ===========================
# Linear fit over AMAV
# ===========================

def fit_linear(x: np.ndarray, y: np.ndarray):
    """
    Fit y ~ a + b*x.
    Returns (x_sorted, y_pred, p_value). If not enough points, x_sorted is empty.
    """
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if len(y) < 3:
        return np.array([]), np.array([]), np.nan

    x_center = x - x.mean()
    X = sm.add_constant(x_center)
    model = sm.OLS(y, X)
    res = model.fit()

    order = np.argsort(x)
    x_sorted = x[order]
    X_sorted = sm.add_constant(x_sorted - x.mean())
    y_pred = res.predict(X_sorted)

    p_val = res.pvalues[1] if len(res.params) > 1 else np.nan
    return x_sorted, y_pred, p_val


# ===========================
# Plot per disease
# ===========================

def plot_disease(
    disease: str,
    grp: pd.DataFrame,
    years_all: List[int],
    out_dir: Path,
):
    # Prevalence per study
    trends = []
    for _, row in grp.iterrows():
        vals = row[years_all].to_numpy(dtype=float)
        trends.append(vals)

    if len(trends) == 0:
        return

    years_all_arr = np.array(years_all, dtype=float)

    # Years where at least one study has prevalence data
    vals_stack = np.vstack(trends)  # (n_trends, n_years)
    mask_any = np.isfinite(vals_stack).any(axis=0)
    if not mask_any.any():
        return

    years_valid = years_all_arr[mask_any]
    if len(years_valid) < 3:
        return  # too few points

    # Rebuild raw_df only with valid years
    raw = pd.DataFrame({"Year": years_valid})
    for i, vals in enumerate(trends, start=1):
        raw[f"T{i}"] = vals[mask_any]

    if raw.drop(columns=["Year"]).isna().all().all():
        return

    # Slopes + AMAV
    slope_df = build_perstudy_linearslope(raw)
    summary = build_summary_mav_amav(slope_df)

    years_vec = summary["Year"].to_numpy(dtype=float)
    amav_all = summary["AMAV"].to_numpy(dtype=float)
    has_pos = "AMAV-POS" in summary.columns and summary["AMAV-POS"].notna().any()
    amav_pos = summary["AMAV-POS"].to_numpy(dtype=float) if has_pos else None

    # Series used for the fit
    if has_pos:
        y_for_fit = amav_pos
        fit_label = "Linear fit (POS)"
    else:
        y_for_fit = amav_all
        fit_label = "Linear fit"

    # Linear fit
    x_fit, y_fit, p_val = fit_linear(years_vec, y_for_fit)

    # Optional overrides
    overrides = PLOT_OVERRIDES.get(disease, {})
    xmin = overrides.get("xmin", years_valid.min() - 1)
    xmax = overrides.get("xmax", years_valid.max() + 1)

    # ---------- Plot ----------
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, ax_left = plt.subplots(figsize=(8, 5))
    ax_right = ax_left.twinx()

    # Grey lines: prevalence per study (with points at observed years)
    for vals in trends:
        # Restrict to years where at least one study has data
        vals_valid = vals[mask_any]              # shape = len(years_valid)
        mask_line = np.isfinite(vals_valid)

        # Only draw lines with ≥ 2 points
        if mask_line.sum() >= 2:
            other_years = years_valid[mask_line]
            other_vals  = vals_valid[mask_line]

            ax_left.plot(
                other_years,
                other_vals,
                color="lightgray",
                linewidth=1,
                marker="o",
                markersize=2,
                alpha=0.5,
            )

    # AMAV (ALL) in red
    mask_amav = np.isfinite(amav_all)
    if mask_amav.any():
        ax_right.plot(
            years_vec[mask_amav],
            amav_all[mask_amav],
            color="red",
            marker="o",
            linewidth=2,
            label="AMAV (ALL)",
        )

    # AMAV-POS in blue if present
    if has_pos:
        mask_pos = np.isfinite(amav_pos)
        if mask_pos.any():
            ax_right.plot(
                years_vec[mask_pos],
                amav_pos[mask_pos],
                color="blue",
                marker="o",
                linewidth=2,
                label="AMAV-POS",
            )

    # Linear fit in black
    if len(x_fit) > 0:
        ax_right.plot(
            x_fit,
            y_fit,
            color="black",
            linestyle="-",
            linewidth=2,
            label=fit_label,
        )

    # Labels and limits
    ax_left.set_xlabel("Year")
    ax_left.set_ylabel("Prevalence (original trends)", color="0.3")
    ax_right.set_ylabel("AMAV / AMAV-POS (arbitrary units)", color="0.3")
    ax_left.set_title(disease, fontweight="bold")
    ax_left.set_xlim(xmin, xmax)

    # Optional custom y-limits
    if "ymin_left" in overrides or "ymax_left" in overrides:
        ymin_l, ymax_l = ax_left.get_ylim()
        ymin_l = overrides.get("ymin_left", ymin_l)
        ymax_l = overrides.get("ymax_left", ymax_l)
        ax_left.set_ylim(ymin_l, ymax_l)

     # P-value text in bold
    if not np.isnan(p_val):
        ax_right.text(
            0.02, 0.95,
            rf"$\mathbf{{p\ (linear)\ =\ {p_val:.2e}}}$",
            transform=ax_right.transAxes,
            ha="left", va="top",
            fontsize=9,
        )
    
    # Legend (only right-axis series)
    lines, labels = ax_right.get_legend_handles_labels()
    if lines:
        ax_right.legend(lines, labels, loc="lower right", fontsize=9)

    fig.tight_layout()
    fname = out_dir / f"{slugify(disease)}.png"
    fig.savefig(fname, dpi=300)
    plt.close(fig)
    print(f"Saved figure for {disease!r} -> {fname}")


# ===========================
# MAIN
# ===========================

def main():
    if not INPUT_XLSX.exists():
        raise SystemExit(f"Input file not found: {INPUT_XLSX}")

    df, year_cols = read_single_table(INPUT_XLSX)

    for disease, grp in df.groupby("Phenotype", dropna=True):
        plot_disease(disease, grp, year_cols, OUT_DIR)


if __name__ == "__main__":
    main()
