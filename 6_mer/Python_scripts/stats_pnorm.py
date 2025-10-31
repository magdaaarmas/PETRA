import os
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from scipy.stats import norm

########################################################################################################################

path=os.getcwd()
gene=sys.argv[1]
path_filtered=f"{path}/filtered_scores/"
output_path=f"{path_filtered}stats_pnorm/"
os.makedirs(output_path, exist_ok=True)


########################################################################################################################
def open_pickle(file_path: str) -> pd.DataFrame:
    with open(file_path, "rb") as fp:
        return pickle.load(fp)

def pnorm_2_sided(df: pd.DataFrame, score_col: str) -> pd.DataFrame:
    """Two-sided normal p-values for the given column vs its own mean/std."""
    x = df[score_col].astype(float)
    mu = x.mean()
    sd = x.std(ddof=1)
    if not np.isfinite(sd) or sd == 0:
        # avoid divide-by-zero; assign p=1 if no variance
        pvals = np.ones(len(x))
    else:
        z = (x - mu) / sd
        pvals = 2 * norm.sf(np.abs(z))
    out = df.copy()
    out["p_val"] = pvals
    out["sig_05"] = out["p_val"] < 0.05
    out["sig_01"] = out["p_val"] < 0.01
    return out

def qq_plot(df: pd.DataFrame, pval_col: str, out_png: str, title: str):
    """Simple QQ plot of p-values."""
    p = df[pval_col].dropna().values
    if p.size == 0:
        return
    p.sort()
    n = p.size
    expected = -np.log10(np.linspace(1/(n+1), 1, n))
    observed = -np.log10(p)

    plt.figure(figsize=(6, 6))
    plt.plot(expected, observed, "o", markersize=3, label="Observed")
    m = max(expected.max(), observed.max())
    plt.plot([0, m], [0, m], "--", color="lightgrey", label="Expected")
    plt.xlabel("Expected -log10(p)")
    plt.ylabel("Observed -log10(p)")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
########################################################################################################################

scores=open_pickle(f"{path_filtered}{gene}_filtered_scores.pkl")
statistics=pnorm_2_sided(scores, 'filtered_score')
qq_plot(statistics, 'p_val', 2,f"QQ plot for {gene} pnorm")
statistics.to_csv(f"{path_filtered}{gene}_filtered_scores_statistics.csv")




