#| label: LINEAR_vs_EXPONENTIAL.py
import pandas as pd
import numpy as np
import statsmodels.api as sm
from pathlib import Path
import math

# ===========================
# CONFIG
# ===========================
# We read directly from AMAV_DATA.xlsx, sheet "AMAV-P"
# Sheet layout: Year | Disease1 | Disease2 | ... | DiseaseN
excel_path = Path("data/AMAV_DATA.xlsx")
sheet_name = "AMAV"

# Minimum number of points
MIN_POINTS_LINEAR = 2   # minimum to fit a line
MIN_POINTS_EXP    = 3   # minimum to fit an exponential
MIN_POINTS_CV     = 5   # minimum to run LOOCV reliably

# Decision rule for preferring exponential
AIC_MARGIN = 4.0        # exponential must improve AIC by at least this
RMSE_RATIO = 0.98       # exponential RMSE must be <= 0.98 * linear RMSE


# ===========================
# Helpers: fitting
# ===========================

def fit_linear(x, y):
    """Fit a simple linear model: y ~ a + b * x."""
    X = sm.add_constant(x)
    res = sm.OLS(y, X).fit()
    return res, res.aic, res.bic


def fit_exponential(x, y):
    """
    Fit an exponential model y = a * exp(b * x).
    Implemented as log(y) ~ A + B * x.
    """
    if np.any(y <= 0):
        raise ValueError("Exponential model requires y > 0.")
    logy = np.log(y)
    X = sm.add_constant(x)
    res = sm.OLS(logy, X).fit()
    return res, res.aic, res.bic


# ===========================
# Helpers: LOOCV
# ===========================

def loocv_rmse_linear(x, y):
    """
    Leave-One-Out Cross-Validation RMSE for the linear model.
    Uses manual prediction (beta0 + beta1 * x_test) to avoid
    shape issues with statsmodels.predict.
    """
    n = len(y)
    if n < MIN_POINTS_CV:
        return np.nan

    se_sum = 0.0
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        x_train, y_train = x[mask], y[mask]
        x_test,  y_test  = x[~mask], y[~mask]

        res, _, _ = fit_linear(x_train, y_train)
        beta0, beta1 = res.params
        y_pred = beta0 + beta1 * x_test
        se_sum += np.sum((y_test - y_pred) ** 2)

    rmse = math.sqrt(se_sum / n)
    return rmse


def loocv_rmse_exponential(x, y):
    """
    Leave-One-Out Cross-Validation RMSE for the exponential model.
    RMSE is computed on the original y scale.
    """
    n = len(y)
    if n < MIN_POINTS_CV or np.any(y <= 0):
        return np.nan

    se_sum = 0.0
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        x_train, y_train = x[mask], y[mask]
        x_test,  y_test  = x[~mask], y[~mask]

        res, _, _ = fit_exponential(x_train, y_train)
        beta0, beta1 = res.params
        log_pred = beta0 + beta1 * x_test
        y_pred = np.exp(log_pred)
        se_sum += np.sum((y_test - y_pred) ** 2)

    rmse = math.sqrt(se_sum / n)
    return rmse


# ===========================
# Model choice per disease
# ===========================

def choose_model_for_disease(df_one):
    """
    Decide whether AMAV (already AMAV-POS in AMAV_DATA.xlsx)
    for one disease is better described by a linear or exponential trend.
    df_one has columns: Year, Disease, AMAV.
    """
    n_all = len(df_one)
    n_pos = n_all  # all values are already AMAV-POS

    if n_all < MIN_POINTS_LINEAR:
        return {
            "n_all": n_all,
            "n_pos": n_pos,
            "model": "INSUFFICIENT_DATA",
            "AIC_linear": np.nan,
            "AIC_exponential": np.nan,
            "BIC_linear": np.nan,
            "BIC_exponential": np.nan,
            "RMSEcv_linear": np.nan,
            "RMSEcv_exponential": np.nan,
        }

    # Centre years for numerical stability
    x = df_one["Year"].to_numpy(dtype=float)
    x = x - np.mean(x)
    y = df_one["AMAV"].to_numpy(dtype=float)

    # If we have too few points for exponential, fit only linear
    if n_pos < MIN_POINTS_EXP:
        res_lin, aic_lin, bic_lin = fit_linear(x, y)
        return {
            "n_all": n_all,
            "n_pos": n_pos,
            "model": "LINEAR_ONLY",
            "AIC_linear": aic_lin,
            "AIC_exponential": np.nan,
            "BIC_linear": bic_lin,
            "BIC_exponential": np.nan,
            "RMSEcv_linear": np.nan,
            "RMSEcv_exponential": np.nan,
        }

    # Linear model
    res_lin, aic_lin, bic_lin = fit_linear(x, y)
    rmse_lin = loocv_rmse_linear(x, y)

    # Exponential model
    try:
        res_exp, aic_exp, bic_exp = fit_exponential(x, y)
        rmse_exp = loocv_rmse_exponential(x, y)
    except ValueError:
        # Some y <= 0 (should not happen with AMAV-POS, but keep safe)
        aic_exp = bic_exp = rmse_exp = np.nan
        chosen = "LINEAR"
    else:
        delta_aic = aic_exp - aic_lin  # negative favours exponential

        if np.isnan(rmse_lin) or np.isnan(rmse_exp):
            # Decide using only AIC if CV is not available
            chosen = "EXPONENTIAL" if delta_aic <= -AIC_MARGIN else "LINEAR"
        else:
            rmse_ratio = rmse_exp / rmse_lin if rmse_lin > 0 else math.inf
            if (delta_aic <= -AIC_MARGIN) and (rmse_ratio <= RMSE_RATIO):
                chosen = "EXPONENTIAL"
            else:
                chosen = "LINEAR"

    return {
        "n_all": n_all,
        "n_pos": n_pos,
        "model": chosen,
        "AIC_linear": aic_lin,
        "AIC_exponential": aic_exp,
        "BIC_linear": bic_lin,
        "BIC_exponential": bic_exp,
        "RMSEcv_linear": rmse_lin,
        "RMSEcv_exponential": rmse_exp,
    }


# ===========================
# MAIN
# ===========================

def main():
    # Read AMAV sheet (wide format: one column per disease)
    wide = pd.read_excel(excel_path, sheet_name=sheet_name)

    # Expect a 'Year' column + disease columns
    if "Year" not in wide.columns:
        raise ValueError("Column 'Year' not found in AMAV sheet.")

    # Convert to long format: Year | Disease | AMAV
    long_df = wide.melt(
        id_vars=["Year"],
        var_name="Disease",
        value_name="AMAV"
    )

    # Drop missing AMAV values
    long_df = long_df.dropna(subset=["AMAV"])

    results = []
    for disease, df_one in long_df.groupby("Disease"):
        metrics = choose_model_for_disease(df_one)
        metrics["Disease"] = disease
        results.append(metrics)

    res_df = pd.DataFrame(results)

    # Order columns nicely
    cols_order = [
        "Disease", "model", "n_all", "n_pos",
        "AIC_linear", "AIC_exponential",
        "BIC_linear", "BIC_exponential",
        "RMSEcv_linear", "RMSEcv_exponential",
    ]
    cols_order = [c for c in cols_order if c in res_df.columns]
    res_df = res_df[cols_order]

    out_path = Path("output/model_choice_AMAV.csv")
    res_df.to_csv(out_path, index=False)
    print(f"Saved model comparison to: {out_path.resolve()}")


if __name__ == "__main__":
    main()
