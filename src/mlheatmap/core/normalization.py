"""Normalization methods for RNA-seq count data."""

import numpy as np


def deseq2_normalize(counts: np.ndarray, return_size_factors: bool = False):
    """DESeq2-like median-of-ratios normalization + VST.

    Args:
        counts: Raw count matrix (genes x samples)
        return_size_factors: If True, return (vst, size_factors) tuple.

    Returns:
        VST-transformed matrix (genes x samples), or
        (vst, size_factors) if return_size_factors=True.
    """
    counts = counts.astype(np.float64)

    # Step 1: Geometric mean per gene (log-space to avoid overflow)
    log_counts = np.full(counts.shape, np.nan, dtype=np.float64)
    positive_mask = counts > 0
    log_counts[positive_mask] = np.log(counts[positive_mask])
    geo_mean_log = np.nanmean(log_counts, axis=1)
    valid = np.isfinite(geo_mean_log)
    geo_means = np.exp(geo_mean_log)

    # Step 2: Ratios to reference
    ratios = np.where(
        (counts[valid] > 0) & (geo_means[valid, np.newaxis] > 0),
        counts[valid] / geo_means[valid, np.newaxis],
        np.nan,
    )

    # Step 3: Size factors = median of ratios per sample
    size_factors = np.nanmedian(ratios, axis=0)
    size_factors = np.where(size_factors > 0, size_factors, 1.0)

    # Step 4: Normalized counts
    normalized = counts / size_factors[np.newaxis, :]

    # Step 5: VST approximation
    vst = _vst_transform(normalized)

    if return_size_factors:
        return vst, size_factors
    return vst


def _vst_transform(normalized: np.ndarray) -> np.ndarray:
    """Variance-stabilizing transform (simplified).

    Uses dispersion estimation + arcsinh transform.
    """
    n_samples = normalized.shape[1]
    gene_means = np.mean(normalized, axis=1)
    gene_vars = np.var(normalized, axis=1, ddof=1) if n_samples > 1 else np.zeros(normalized.shape[0])

    # Fit mean-variance: var ≈ mean + alpha * mean^2
    valid = gene_means > 0.5
    if np.sum(valid) > 100:
        x = gene_means[valid]
        y = gene_vars[valid]
        alpha_est = (y - x) / (x ** 2 + 1e-6)
        pos = alpha_est > 0
        if np.sum(pos) > 10:
            alpha = float(np.median(alpha_est[pos]))
        else:
            alpha = 0.05
        alpha = max(min(alpha, 10.0), 0.01)
    else:
        alpha = 0.05

    # VST: 2/sqrt(alpha) * arcsinh(sqrt(alpha * x))
    result = (2.0 / np.sqrt(alpha)) * np.arcsinh(np.sqrt(alpha * np.maximum(normalized, 0)))
    return result


def tpm_abundance(counts: np.ndarray, gene_lengths: np.ndarray = None) -> np.ndarray:
    """Return linear TPM abundance or CPM fallback when gene lengths are absent."""
    counts = counts.astype(np.float64)

    if gene_lengths is not None and len(gene_lengths) == counts.shape[0]:
        rpk = counts / (gene_lengths[:, np.newaxis] / 1000.0 + 1e-6)
        scale = np.sum(rpk, axis=0) / 1e6
        scale = np.maximum(scale, 1e-6)
        return rpk / scale[np.newaxis, :]

    lib_sizes = np.sum(counts, axis=0)
    lib_sizes = np.maximum(lib_sizes, 1.0)
    return counts / lib_sizes[np.newaxis, :] * 1e6


def tpm_normalize(counts: np.ndarray, gene_lengths: np.ndarray = None) -> np.ndarray:
    """TPM normalization with log2 transform.

    If gene_lengths not provided, uses CPM + log2 instead.
    """
    tpm = tpm_abundance(counts, gene_lengths=gene_lengths)
    return np.log2(tpm + 1)


def log2_normalize(counts: np.ndarray) -> np.ndarray:
    """Simple log2(count + 1) normalization."""
    return np.log2(counts.astype(np.float64) + 1)
