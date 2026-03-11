"""Differential Expression Gene (DEG) Analysis.

Implements statistical tests for identifying differentially expressed genes
between two or more groups.
"""

import numpy as np
from scipy import stats


def compute_deg(
    expression: np.ndarray,
    gene_names: list,
    sample_groups: dict,
    method: str = "wilcoxon",
    log2fc_threshold: float = 1.0,
    pvalue_threshold: float = 0.05,
) -> dict:
    """Run DEG analysis between groups.

    Parameters
    ----------
    expression : np.ndarray
        Normalized expression matrix (genes × samples).
    gene_names : list
        Gene symbol list matching rows of expression.
    sample_groups : dict
        {group_name: [sample_index, ...]}. Must have exactly 2 groups.
    method : str
        Statistical test: 'wilcoxon' (Mann-Whitney U) or 'ttest'.
    log2fc_threshold : float
        |log2FC| threshold for significance (default 1.0).
    pvalue_threshold : float
        Adjusted p-value threshold (default 0.05).

    Returns
    -------
    dict with keys:
        - results: list of per-gene dicts sorted by adj_pvalue
        - summary: dict with counts of up/down/not significant
        - thresholds: dict with fc/pvalue cutoffs used
    """
    group_names = list(sample_groups.keys())
    if len(group_names) != 2:
        raise ValueError(f"DEG requires exactly 2 groups, got {len(group_names)}")

    idx_g1 = np.array(sample_groups[group_names[0]])
    idx_g2 = np.array(sample_groups[group_names[1]])

    n_genes = expression.shape[0]
    results = []

    for i in range(n_genes):
        expr_g1 = expression[i, idx_g1]
        expr_g2 = expression[i, idx_g2]

        # Mean expression per group
        mean_g1 = float(np.mean(expr_g1))
        mean_g2 = float(np.mean(expr_g2))

        # log2 Fold Change: group1 vs group2
        # If data is already log-transformed, difference = log2FC
        # If not, compute from raw means with pseudocount
        log2fc = mean_g1 - mean_g2

        # Statistical test
        pval = 1.0
        if method == "wilcoxon":
            try:
                _, pval = stats.mannwhitneyu(
                    expr_g1, expr_g2, alternative="two-sided"
                )
            except ValueError:
                pval = 1.0
        elif method == "ttest":
            try:
                _, pval = stats.ttest_ind(expr_g1, expr_g2, equal_var=False)
                if np.isnan(pval):
                    pval = 1.0
            except Exception:
                pval = 1.0

        results.append(
            {
                "gene": gene_names[i],
                "gene_idx": i,
                "log2fc": float(log2fc),
                "pvalue": float(pval),
                "mean_g1": mean_g1,
                "mean_g2": mean_g2,
            }
        )

    # Multiple testing correction (Benjamini-Hochberg)
    pvals = np.array([r["pvalue"] for r in results])
    adj_pvals = _benjamini_hochberg(pvals)

    n_up = 0
    n_down = 0
    n_ns = 0

    for i, r in enumerate(results):
        r["adj_pvalue"] = float(adj_pvals[i])
        r["neg_log10_p"] = float(-np.log10(max(r["adj_pvalue"], 1e-300)))

        # Classification
        is_sig = r["adj_pvalue"] < pvalue_threshold
        if is_sig and r["log2fc"] > log2fc_threshold:
            r["direction"] = "up"
            n_up += 1
        elif is_sig and r["log2fc"] < -log2fc_threshold:
            r["direction"] = "down"
            n_down += 1
        else:
            r["direction"] = "ns"
            n_ns += 1

    # Sort by adjusted p-value
    results.sort(key=lambda r: r["adj_pvalue"])

    return {
        "results": results,
        "summary": {
            "n_up": n_up,
            "n_down": n_down,
            "n_not_significant": n_ns,
            "n_total": n_genes,
        },
        "group_names": group_names,
        "thresholds": {
            "log2fc": log2fc_threshold,
            "pvalue": pvalue_threshold,
        },
        "method": method,
    }


def _benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    sorted_idx = np.argsort(pvals)
    sorted_pvals = pvals[sorted_idx]

    # BH adjustment: p_adj[i] = p[i] * n / rank
    adj_pvals = np.zeros(n)
    for i in range(n):
        adj_pvals[sorted_idx[i]] = sorted_pvals[i] * n / (i + 1)

    # Enforce monotonicity (working backwards through sorted order)
    adj_sorted = adj_pvals[sorted_idx].copy()
    for i in range(n - 2, -1, -1):
        adj_sorted[i] = min(adj_sorted[i], adj_sorted[i + 1])
    adj_pvals[sorted_idx] = adj_sorted

    return np.minimum(adj_pvals, 1.0)
