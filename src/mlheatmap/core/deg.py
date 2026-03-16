"""Differential Expression Gene (DEG) Analysis.

Implements statistical tests for identifying differentially expressed genes
between two or more groups.
"""

import numpy as np
from scipy import stats

from mlheatmap.core.cancellation import raise_if_cancelled

MIN_POSITIVE_PVALUE = float(np.nextafter(0.0, 1.0))


def compute_deg(
    expression: np.ndarray,
    gene_names: list,
    sample_groups: dict,
    method: str = "wilcoxon",
    log2fc_threshold: float = 1.0,
    pvalue_threshold: float = 0.05,
    use_raw_pvalue: bool = False,
    effect_size_data: np.ndarray = None,
    effect_size_basis: str = "",
    cancel_check=None,
) -> dict:
    """Run DEG analysis between groups.

    Parameters
    ----------
    expression : np.ndarray
        Normalized expression matrix (genes × samples).
        Used for the statistical test.
    gene_names : list
        Gene symbol list matching rows of expression.
    sample_groups : dict
        {group_name: [sample_index, ...]}. Must have exactly 2 groups.
    method : str
        Statistical test: 'wilcoxon' (Mann-Whitney U) or 'ttest'.
    log2fc_threshold : float
        |log2FC| threshold for significance (default 1.0).
    pvalue_threshold : float
        P-value threshold for significance (default 0.05).
    use_raw_pvalue : bool
        If True, use raw p-value instead of FDR-adjusted for significance.
    effect_size_data : np.ndarray, optional
        Linear-scale abundance matrix (genes × samples) used to compute
        effect sizes as log2(mean+1) differences. For example:
        raw counts, TPM, or size-factor-normalized counts.
    effect_size_basis : str
        Human-readable label for `effect_size_data`.

    Returns
    -------
    dict with keys:
        - results: list of per-gene dicts sorted by the chosen p-value
        - summary: dict with counts of up/down/not significant
        - thresholds: dict with fc/pvalue cutoffs used
        - comparison_group: str, group_names[0] — numerator in log2FC
        - reference_group: str, group_names[1] — denominator in log2FC
    """
    group_names = list(sample_groups.keys())
    if len(group_names) != 2:
        raise ValueError(f"DEG requires exactly 2 groups, got {len(group_names)}")

    idx_g1 = np.array(sample_groups[group_names[0]])
    idx_g2 = np.array(sample_groups[group_names[1]])

    n_genes = expression.shape[0]
    results = []

    for i in range(n_genes):
        if i % 128 == 0:
            raise_if_cancelled(cancel_check)
        expr_g1 = expression[i, idx_g1]
        expr_g2 = expression[i, idx_g2]

        # Replace non-finite values (from log2(0) = -inf, etc.)
        expr_g1 = np.nan_to_num(expr_g1, nan=0.0, posinf=0.0, neginf=0.0)
        expr_g2 = np.nan_to_num(expr_g2, nan=0.0, posinf=0.0, neginf=0.0)

        # Mean expression per group (normalized, for reporting)
        mean_g1 = float(np.mean(expr_g1))
        mean_g2 = float(np.mean(expr_g2))

        # log2 Fold Change always uses a linear-scale abundance basis when available.
        if effect_size_data is not None:
            basis_g1 = np.nan_to_num(effect_size_data[i, idx_g1].astype(np.float64),
                                     nan=0.0, posinf=0.0, neginf=0.0)
            basis_g2 = np.nan_to_num(effect_size_data[i, idx_g2].astype(np.float64),
                                     nan=0.0, posinf=0.0, neginf=0.0)
            mean_basis_g1 = float(np.mean(basis_g1))
            mean_basis_g2 = float(np.mean(basis_g2))
            log2fc = float(np.log2(mean_basis_g1 + 1) - np.log2(mean_basis_g2 + 1))
        else:
            # Fallback only for legacy callers that do not provide a linear basis.
            log2fc = mean_g1 - mean_g2
        if not np.isfinite(log2fc):
            log2fc = 0.0

        # Statistical test
        pval = 1.0
        if method == "wilcoxon":
            try:
                _, pval = stats.mannwhitneyu(
                    expr_g1, expr_g2, alternative="two-sided"
                )
                if not np.isfinite(pval):
                    pval = 1.0
            except ValueError:
                pval = 1.0
        elif method == "ttest":
            try:
                _, pval = stats.ttest_ind(expr_g1, expr_g2, equal_var=False)
                if not np.isfinite(pval):
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
    raise_if_cancelled(cancel_check)
    pvals = np.array([r["pvalue"] for r in results])
    adj_pvals = _benjamini_hochberg(pvals)

    n_up = 0
    n_down = 0
    n_ns = 0

    for i, r in enumerate(results):
        r["adj_pvalue"] = float(adj_pvals[i])

        # Use raw or adjusted p-value for significance and display
        p_for_sig = r["pvalue"] if use_raw_pvalue else r["adj_pvalue"]
        clipped_p = float(np.clip(p_for_sig, MIN_POSITIVE_PVALUE, 1.0))
        r["neg_log10_p"] = float(-np.log10(clipped_p))

        # Classification
        is_sig = p_for_sig < pvalue_threshold
        if is_sig and r["log2fc"] > log2fc_threshold:
            r["direction"] = "up"
            n_up += 1
        elif is_sig and r["log2fc"] < -log2fc_threshold:
            r["direction"] = "down"
            n_down += 1
        else:
            r["direction"] = "ns"
            n_ns += 1

    # Sort by the chosen p-value type
    sort_key = "pvalue" if use_raw_pvalue else "adj_pvalue"
    results.sort(key=lambda r: r[sort_key])

    return {
        "results": results,
        "summary": {
            "n_up": n_up,
            "n_down": n_down,
            "n_not_significant": n_ns,
            "n_total": n_genes,
        },
        "group_names": group_names,
        "comparison_group": group_names[0],
        "reference_group": group_names[1],
        "thresholds": {
            "log2fc": log2fc_threshold,
            "pvalue": pvalue_threshold,
        },
        "method": method,
        "pvalue_type": "raw" if use_raw_pvalue else "fdr",
        "effect_size_basis": effect_size_basis or "normalized_expression",
    }


def _benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    if n == 0:
        return np.zeros(0, dtype=np.float64)
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
