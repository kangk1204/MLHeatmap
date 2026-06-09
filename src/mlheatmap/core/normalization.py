"""Normalization methods for RNA-seq count data.

Two usage modes are supported:

* Global normalization (the historical default): ``deseq2_normalize`` /
  ``tpm_normalize`` / ``log2_normalize`` compute every statistic from the matrix
  they are given. This is what the application calls on the full uploaded matrix.
* Fold-isolated normalization (``FoldNormalizer``): the cross-sample statistics
  (DESeq2-like size-factor reference and the VST dispersion) are *fit* on a
  training matrix and then *applied* unchanged to held-out samples. This lets the
  biomarker module normalize inside each cross-validation fold so that test-fold
  samples never influence the transform applied to training samples.

The global helpers are implemented on top of the same primitives as
``FoldNormalizer`` so that ``deseq2_normalize(X)`` is numerically identical to
``FoldNormalizer("deseq2").fit_transform(X)``.
"""

from __future__ import annotations

import warnings

import numpy as np


# ---------------------------------------------------------------------------
# DESeq2-like primitives (median-of-ratios size factors + simplified VST)
# ---------------------------------------------------------------------------

def _size_factor_reference(counts: np.ndarray) -> np.ndarray:
    """Per-gene geometric-mean reference (log-space mean of positive counts).

    Returns an array of length ``n_genes``; entries are ``nan`` for genes that
    are zero across every sample in ``counts`` (i.e. carry no reference).
    """
    log_counts = np.full(counts.shape, np.nan, dtype=np.float64)
    positive_mask = counts > 0
    log_counts[positive_mask] = np.log(counts[positive_mask])
    with warnings.catch_warnings():
        # All-zero genes produce an empty slice; the resulting nan reference is
        # intended (those genes simply carry no median-of-ratios contribution).
        warnings.simplefilter("ignore", category=RuntimeWarning)
        geo_mean_log = np.nanmean(log_counts, axis=1)
    return np.exp(geo_mean_log)


def _size_factors_from_reference(counts: np.ndarray, geo_means: np.ndarray) -> np.ndarray:
    """Median-of-ratios size factor per sample, relative to ``geo_means``.

    Only genes with a positive reference and a positive count contribute, which
    matches DESeq2's behaviour and lets the same reference be reused for held-out
    samples without leaking their distribution into the reference itself.
    """
    gm = geo_means[:, np.newaxis]
    ratios = np.where((counts > 0) & (gm > 0), counts / gm, np.nan)
    size_factors = np.nanmedian(ratios, axis=0)
    return np.where(size_factors > 0, size_factors, 1.0)


def _estimate_vst_alpha(normalized: np.ndarray) -> float:
    """Estimate the single global dispersion used by the simplified VST."""
    n_samples = normalized.shape[1]
    gene_means = np.mean(normalized, axis=1)
    gene_vars = np.var(normalized, axis=1, ddof=1) if n_samples > 1 else np.zeros(normalized.shape[0])

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
    return alpha


def _apply_vst(normalized: np.ndarray, alpha: float) -> np.ndarray:
    """Variance-stabilizing transform: 2/sqrt(alpha) * arcsinh(sqrt(alpha * x))."""
    return (2.0 / np.sqrt(alpha)) * np.arcsinh(np.sqrt(alpha * np.maximum(normalized, 0)))


def deseq2_normalize(counts: np.ndarray, return_size_factors: bool = False):
    """DESeq2-like median-of-ratios normalization + VST (global).

    Args:
        counts: Raw count matrix (genes x samples)
        return_size_factors: If True, return (vst, size_factors) tuple.

    Returns:
        VST-transformed matrix (genes x samples), or
        (vst, size_factors) if return_size_factors=True.
    """
    counts = counts.astype(np.float64)
    geo_means = _size_factor_reference(counts)
    size_factors = _size_factors_from_reference(counts, geo_means)
    normalized = counts / size_factors[np.newaxis, :]
    vst = _apply_vst(normalized, _estimate_vst_alpha(normalized))

    if return_size_factors:
        return vst, size_factors
    return vst


def _vst_transform(normalized: np.ndarray) -> np.ndarray:
    """Variance-stabilizing transform (simplified): dispersion fit + arcsinh."""
    return _apply_vst(normalized, _estimate_vst_alpha(normalized))


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


# ---------------------------------------------------------------------------
# Fold-isolated normalization (fit on training samples, apply to held-out ones)
# ---------------------------------------------------------------------------

NORMALIZATION_METHODS = ("deseq2", "tpm", "log2")


class FoldNormalizer:
    """Fit cross-sample normalization parameters on training data, then apply
    them unchanged to any matrix (e.g. a held-out CV fold).

    For ``deseq2`` the fitted parameters are the per-gene geometric-mean
    reference and the single VST dispersion ``alpha``. For ``tpm`` and ``log2``
    the transform is fully per-sample, so ``fit`` stores nothing and
    ``transform`` is identical to the global function (kept here so callers can
    treat all methods uniformly).

    Matrices are ``genes x samples`` in every method, matching the global
    helpers above.
    """

    def __init__(self, method: str = "deseq2"):
        if method not in NORMALIZATION_METHODS:
            raise ValueError(f"Unknown method: {method}. Choose from {NORMALIZATION_METHODS}.")
        self.method = method
        self._geo_means: np.ndarray | None = None
        self._alpha: float | None = None
        self._fitted = False

    def fit(self, train_counts: np.ndarray) -> "FoldNormalizer":
        train = np.asarray(train_counts, dtype=np.float64)
        if self.method == "deseq2":
            self._geo_means = _size_factor_reference(train)
            size_factors = _size_factors_from_reference(train, self._geo_means)
            normalized = train / size_factors[np.newaxis, :]
            self._alpha = _estimate_vst_alpha(normalized)
        self._fitted = True
        return self

    def transform(self, counts: np.ndarray) -> np.ndarray:
        if not self._fitted:
            raise RuntimeError("FoldNormalizer.transform called before fit().")
        c = np.asarray(counts, dtype=np.float64)
        if self.method == "deseq2":
            size_factors = _size_factors_from_reference(c, self._geo_means)
            normalized = c / size_factors[np.newaxis, :]
            return _apply_vst(normalized, self._alpha)
        if self.method == "tpm":
            return tpm_normalize(c)
        return log2_normalize(c)

    def fit_transform(self, counts: np.ndarray) -> np.ndarray:
        return self.fit(counts).transform(counts)
