"""Biomarker discovery using ML + SHAP with nested outer CV."""

from __future__ import annotations

import numpy as np
from sklearn.metrics import auc, roc_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder, label_binarize, StandardScaler

from mlheatmap.core.capabilities import normalize_model_name
from mlheatmap.core.cancellation import raise_if_cancelled

PREFILTER_MAX_GENES = 2000
MAX_PANEL_GENES = 15


def _build_model(model_name: str, n_estimators: int = 500, n_samples: int = 100):
    """Return (clf, use_shap_tree, needs_scaling) for the requested model."""
    model_name = normalize_model_name(model_name)

    if model_name in ("rf", "random_forest"):
        from sklearn.ensemble import RandomForestClassifier

        max_d = 3 if n_samples < 20 else (6 if n_samples < 50 else None)
        clf = RandomForestClassifier(
            n_estimators=n_estimators,
            max_depth=max_d,
            min_samples_leaf=max(2, n_samples // 10) if n_samples < 30 else max(1, n_samples // 20),
            class_weight="balanced",
            random_state=42,
            n_jobs=-1,
        )
        return clf, True, False

    if model_name == "xgboost":
        try:
            from xgboost import XGBClassifier
        except ImportError as exc:
            raise ImportError(
                "XGBoost is not available in this installation. "
                "Install mlheatmap[full] or use Docker for the full model set."
            ) from exc
        if n_samples < 20:
            n_est = min(n_estimators, 50)
            max_d, lr = 2, 0.1
            sub, col_sub = 1.0, 0.6
            reg_a, reg_l = 1.0, 5.0
        elif n_samples < 50:
            n_est = min(n_estimators, 200)
            max_d, lr = 3, 0.05
            sub, col_sub = 0.9, 0.7
            reg_a, reg_l = 0.5, 2.0
        else:
            n_est = n_estimators
            max_d, lr = 6, 0.1
            sub, col_sub = 0.8, 0.8
            reg_a, reg_l = 0.1, 1.0
        clf = XGBClassifier(
            n_estimators=n_est,
            max_depth=max_d,
            learning_rate=lr,
            subsample=sub,
            colsample_bytree=col_sub,
            reg_alpha=reg_a,
            reg_lambda=reg_l,
            random_state=42,
            n_jobs=-1,
            verbosity=0,
        )
        return clf, True, False

    if model_name == "lightgbm":
        try:
            from lightgbm import LGBMClassifier
        except ImportError as exc:
            raise ImportError(
                "LightGBM is not available in this installation. "
                "Install mlheatmap[full] or use Docker for the full model set."
            ) from exc
        if n_samples < 20:
            n_est = min(n_estimators, 50)
            n_leaves = 4
            min_child = max(2, n_samples // 5)
            lr = 0.1
            reg_a, reg_l = 1.0, 5.0
            sub, col_sub = 1.0, 0.6
        elif n_samples < 50:
            n_est = min(n_estimators, 150)
            n_leaves = 8
            min_child = max(2, n_samples // 10)
            lr = 0.05
            reg_a, reg_l = 0.5, 2.0
            sub, col_sub = 0.9, 0.7
        else:
            n_est = n_estimators
            n_leaves = 31
            min_child = max(1, n_samples // 20)
            lr = 0.05
            reg_a, reg_l = 0.1, 1.0
            sub, col_sub = 0.8, 0.8
        clf = LGBMClassifier(
            n_estimators=n_est,
            num_leaves=n_leaves,
            max_depth=3 if n_samples < 20 else -1,
            learning_rate=lr,
            subsample=sub,
            subsample_freq=1 if sub < 1.0 else 0,
            colsample_bytree=col_sub,
            min_child_samples=min_child,
            reg_alpha=reg_a,
            reg_lambda=reg_l,
            is_unbalance=True,
            random_state=42,
            n_jobs=-1,
            verbose=-1,
        )
        return clf, True, False

    if model_name == "logistic_regression":
        from sklearn.linear_model import LogisticRegression

        clf = LogisticRegression(
            penalty="l1",
            C=1.0,
            solver="saga",
            max_iter=2000,
            class_weight="balanced",
            random_state=42,
            n_jobs=-1,
        )
        return clf, False, True

    if model_name == "svm_linear":
        from sklearn.svm import SVC

        clf = SVC(
            kernel="linear",
            C=1.0,
            class_weight="balanced",
            probability=True,
            random_state=42,
        )
        return clf, False, True

    raise ValueError(
        f"Unknown model: {model_name}. "
        "Choose from: rf, xgboost, lightgbm, logistic_regression, svm_linear"
    )


def _model_display_name(model_name: str) -> str:
    """Human-readable model name."""
    model_name = normalize_model_name(model_name)
    names = {
        "rf": "Random Forest",
        "xgboost": "XGBoost",
        "lightgbm": "LightGBM",
        "logistic_regression": "Logistic Regression (L1)",
        "svm_linear": "SVM (Linear)",
    }
    return names.get(model_name, model_name)


def _feature_importance_vector(clf) -> np.ndarray:
    if hasattr(clf, "feature_importances_"):
        return np.asarray(clf.feature_importances_, dtype=np.float64)
    if hasattr(clf, "coef_"):
        coef = np.abs(np.asarray(clf.coef_, dtype=np.float64))
        return np.mean(coef, axis=0) if coef.ndim > 1 else coef.ravel()
    return np.ones(getattr(clf, "n_features_in_", 1), dtype=np.float64)


def _predict_proba_full(clf, X: np.ndarray, n_classes: int) -> np.ndarray:
    raw = np.asarray(clf.predict_proba(X), dtype=np.float64)
    if raw.ndim == 1:
        raw = raw.reshape(-1, 1)
    full = np.zeros((raw.shape[0], n_classes), dtype=np.float64)
    clf_classes = np.asarray(getattr(clf, "classes_", np.arange(raw.shape[1])), dtype=int)
    for col_idx, class_idx in enumerate(clf_classes):
        if 0 <= int(class_idx) < n_classes:
            full[:, int(class_idx)] = raw[:, col_idx]
    return full


def _mean_multiclass_auc(y_true: np.ndarray, y_prob: np.ndarray, n_classes: int) -> float:
    y_bin = label_binarize(y_true, classes=range(n_classes))
    if n_classes == 2:
        y_bin = np.column_stack([1 - y_bin.ravel(), y_bin.ravel()])

    aucs = []
    for class_idx in range(n_classes):
        target = y_bin[:, class_idx]
        if target.min() == target.max():
            continue
        try:
            fpr, tpr, _ = roc_curve(target, y_prob[:, class_idx])
        except ValueError:
            continue
        aucs.append(auc(fpr, tpr))
    return float(np.mean(aucs)) if aucs else 0.5


def _compute_test_auc(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    y_test: np.ndarray,
    *,
    model: str,
    n_estimators: int,
    needs_scaling: bool,
    n_classes: int,
) -> float:
    clf, _, _ = _build_model(model, n_estimators, len(y_train))
    if needs_scaling:
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)
    clf.fit(X_train, y_train)
    y_prob = _predict_proba_full(clf, X_test, n_classes)
    return _mean_multiclass_auc(y_test, y_prob, n_classes)


def _evaluate_prefix_auc_curve(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    y_test: np.ndarray,
    ordered_idx: list[int],
    gene_names: list[str],
    *,
    model: str,
    n_estimators: int,
    needs_scaling: bool,
    n_classes: int,
    cancel_check=None,
) -> list[dict[str, float | int | str]]:
    curve = []
    for k in range(1, len(ordered_idx) + 1):
        raise_if_cancelled(cancel_check)
        subset = ordered_idx[:k]
        try:
            score = _compute_test_auc(
                X_train[:, subset],
                y_train,
                X_test[:, subset],
                y_test,
                model=model,
                n_estimators=n_estimators,
                needs_scaling=needs_scaling,
                n_classes=n_classes,
            )
        except Exception:
            score = 0.5
        curve.append(
            {
                "n_genes": k,
                "auc": float(score),
                "gene_added": gene_names[ordered_idx[k - 1]],
            }
        )
    return curve


def _quick_cv_auc(
    X: np.ndarray,
    y: np.ndarray,
    n_classes: int,
    *,
    model: str,
    n_estimators: int,
    needs_scaling: bool,
    cv,
) -> float:
    if cv is None:
        return 0.5
    aucs = []
    quick_n = min(100, n_estimators)
    for train_idx, test_idx in cv.split(X, y):
        try:
            score = _compute_test_auc(
                X[train_idx],
                y[train_idx],
                X[test_idx],
                y[test_idx],
                model=model,
                n_estimators=quick_n,
                needs_scaling=needs_scaling,
                n_classes=n_classes,
            )
        except Exception:
            continue
        aucs.append(score)
    return float(np.mean(aucs)) if aucs else 0.5


def _make_inner_cv(y: np.ndarray, n_splits: int):
    min_class = int(min(np.bincount(y)))
    if min_class < 2:
        return None
    inner_splits = min(n_splits, min_class)
    if inner_splits < 2:
        return None
    return StratifiedKFold(n_splits=inner_splits, shuffle=True, random_state=123)


def _build_inner_auc_curve(
    X_all: np.ndarray,
    y: np.ndarray,
    ordered_idx: list[int],
    gene_names: list[str],
    *,
    model: str,
    n_estimators: int,
    needs_scaling: bool,
    n_splits: int,
    cancel_check=None,
) -> tuple[list[dict[str, float | int | str]], int]:
    n_classes = len(np.unique(y))
    inner_cv = _make_inner_cv(y, n_splits)

    auc_curve = []
    best_auc = -1.0
    best_n = 0
    for k in range(1, len(ordered_idx) + 1):
        raise_if_cancelled(cancel_check)
        subset = ordered_idx[:k]
        try:
            mean_auc = _quick_cv_auc(
                X_all[:, subset],
                y,
                n_classes,
                model=model,
                n_estimators=n_estimators,
                needs_scaling=needs_scaling,
                cv=inner_cv,
            )
        except Exception:
            mean_auc = 0.5
        auc_curve.append(
            {
                "n_genes": k,
                "auc": float(mean_auc),
                "gene_added": gene_names[ordered_idx[k - 1]],
            }
        )
        if mean_auc > best_auc:
            best_auc = mean_auc
            best_n = k
    return auc_curve, best_n


def _forward_order(
    X_all: np.ndarray,
    y: np.ndarray,
    candidate_idx: list[int],
    *,
    model: str,
    n_estimators: int,
    needs_scaling: bool,
    n_splits: int,
    max_genes: int,
    cancel_check=None,
) -> list[int]:
    n_classes = len(np.unique(y))
    inner_cv = _make_inner_cv(y, n_splits)
    if inner_cv is None:
        return list(candidate_idx[:max_genes])

    ordered = []
    pool = list(candidate_idx[:max_genes])
    for _ in range(min(max_genes, len(pool))):
        raise_if_cancelled(cancel_check)
        best_gene = None
        best_auc = -1.0
        for gene_idx in pool:
            if gene_idx in ordered:
                continue
            trial = ordered + [gene_idx]
            try:
                score = _quick_cv_auc(
                    X_all[:, trial],
                    y,
                    n_classes,
                    model=model,
                    n_estimators=n_estimators,
                    needs_scaling=needs_scaling,
                    cv=inner_cv,
                )
            except Exception:
                continue
            if score > best_auc:
                best_auc = score
                best_gene = gene_idx
        if best_gene is None:
            break
        ordered.append(best_gene)
    return ordered


def _lasso_order(
    X_all: np.ndarray,
    y: np.ndarray,
    candidate_idx: list[int],
    *,
    max_genes: int,
    cancel_check=None,
) -> list[int]:
    raise_if_cancelled(cancel_check)
    from sklearn.linear_model import LogisticRegression

    X_cand = X_all[:, candidate_idx]
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_cand)
    n_classes = len(np.unique(y))
    multi = "ovr" if n_classes > 2 else "auto"

    best_sel: list[int] = []
    best_coef = np.zeros(len(candidate_idx), dtype=np.float64)
    for c_val in np.logspace(-3, 1, 20):
        try:
            lr = LogisticRegression(
                penalty="l1",
                C=c_val,
                solver="saga",
                max_iter=5000,
                multi_class=multi,
                class_weight="balanced",
                random_state=42,
            )
            lr.fit(X_scaled, y)
        except Exception:
            continue
        coef = np.abs(np.asarray(lr.coef_, dtype=np.float64))
        coef = np.mean(coef, axis=0) if coef.ndim > 1 else coef.ravel()
        nonzero = np.where(coef > 1e-10)[0]
        if 1 <= len(nonzero) <= max_genes and len(nonzero) > len(best_sel):
            best_sel = nonzero.tolist()
            best_coef = coef

    if best_sel:
        ranked = sorted(best_sel, key=lambda idx: best_coef[idx], reverse=True)
        return [candidate_idx[idx] for idx in ranked]

    try:
        lr = LogisticRegression(
            penalty="l1",
            C=1.0,
            solver="saga",
            max_iter=5000,
            multi_class=multi,
            class_weight="balanced",
            random_state=42,
        )
        lr.fit(X_scaled, y)
        coef = np.abs(np.asarray(lr.coef_, dtype=np.float64))
        coef = np.mean(coef, axis=0) if coef.ndim > 1 else coef.ravel()
        ranked = np.argsort(coef)[::-1][:max_genes]
        return [candidate_idx[idx] for idx in ranked]
    except Exception:
        return list(candidate_idx[:max_genes])


def _stability_order(
    X_all: np.ndarray,
    y: np.ndarray,
    candidate_idx: list[int],
    *,
    max_genes: int,
    cancel_check=None,
) -> list[int]:
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedShuffleSplit

    X_cand = X_all[:, candidate_idx]
    freq = np.zeros(len(candidate_idx), dtype=np.float64)
    n_classes = len(np.unique(y))
    multi = "ovr" if n_classes > 2 else "auto"
    splitter = StratifiedShuffleSplit(n_splits=100, test_size=0.2, random_state=42)

    for seed, (train_idx, _) in enumerate(splitter.split(X_cand, y)):
        raise_if_cancelled(cancel_check)
        X_sub = X_cand[train_idx]
        y_sub = y[train_idx]
        scaler = StandardScaler()
        X_sub = scaler.fit_transform(X_sub)
        try:
            lr = LogisticRegression(
                penalty="l1",
                C=0.5,
                solver="saga",
                max_iter=2000,
                multi_class=multi,
                class_weight="balanced",
                random_state=seed,
            )
            lr.fit(X_sub, y_sub)
        except Exception:
            continue
        coef = np.abs(np.asarray(lr.coef_, dtype=np.float64))
        coef = np.mean(coef, axis=0) if coef.ndim > 1 else coef.ravel()
        freq[coef > 1e-10] += 1

    freq /= 100.0
    threshold = 0.7
    while threshold >= 0.3:
        if np.sum(freq >= threshold) >= 2:
            break
        threshold -= 0.1
    if np.sum(freq >= threshold) < 2:
        ranked = np.argsort(freq)[::-1][:max_genes]
    else:
        stable_idx = np.where(freq >= threshold)[0]
        ranked = sorted(stable_idx, key=lambda idx: freq[idx], reverse=True)[:max_genes]
    return [candidate_idx[idx] for idx in ranked]


def _mrmr_order(
    X_all: np.ndarray,
    y: np.ndarray,
    candidate_idx: list[int],
    *,
    max_genes: int,
    cancel_check=None,
) -> list[int]:
    from sklearn.feature_selection import mutual_info_classif, mutual_info_regression

    X_cand = X_all[:, candidate_idx]
    n_cand = X_cand.shape[1]
    relevance = mutual_info_classif(X_cand, y, random_state=42)
    selected_local: list[int] = []
    remaining = set(range(n_cand))

    for _ in range(min(max_genes, n_cand)):
        raise_if_cancelled(cancel_check)
        if not remaining:
            break
        if not selected_local:
            best_idx = max(remaining, key=lambda idx: relevance[idx])
        else:
            best_idx = None
            best_score = -np.inf
            for cand in remaining:
                redundancy = 0.0
                for sel in selected_local:
                    redundancy += mutual_info_regression(
                        X_cand[:, sel].reshape(-1, 1),
                        X_cand[:, cand],
                        random_state=42,
                    )[0]
                redundancy /= len(selected_local)
                score = relevance[cand] - redundancy
                if score > best_score:
                    best_score = score
                    best_idx = cand
            if best_idx is None:
                break
        selected_local.append(best_idx)
        remaining.discard(best_idx)

    return [candidate_idx[idx] for idx in selected_local]


def _panel_selection_order(
    X_all: np.ndarray,
    y: np.ndarray,
    candidate_idx: list[int],
    gene_names: list[str],
    *,
    panel_method: str,
    model: str,
    n_estimators: int,
    needs_scaling: bool,
    n_splits: int,
    max_genes: int,
    cancel_check=None,
) -> dict[str, object]:
    if panel_method == "lasso":
        ordered_idx = _lasso_order(X_all, y, candidate_idx, max_genes=max_genes, cancel_check=cancel_check)
        method = "lasso"
    elif panel_method == "stability":
        ordered_idx = _stability_order(X_all, y, candidate_idx, max_genes=max_genes, cancel_check=cancel_check)
        method = "stability"
    elif panel_method == "mrmr":
        ordered_idx = _mrmr_order(X_all, y, candidate_idx, max_genes=max_genes, cancel_check=cancel_check)
        method = "mrmr"
    else:
        ordered_idx = _forward_order(
            X_all,
            y,
            candidate_idx,
            model=model,
            n_estimators=n_estimators,
            needs_scaling=needs_scaling,
            n_splits=n_splits,
            max_genes=max_genes,
            cancel_check=cancel_check,
        )
        method = "forward"

    inner_auc_curve, best_n = _build_inner_auc_curve(
        X_all,
        y,
        ordered_idx,
        gene_names,
        model=model,
        n_estimators=n_estimators,
        needs_scaling=needs_scaling,
        n_splits=n_splits,
        cancel_check=cancel_check,
    )
    return {
        "ordered_idx": ordered_idx,
        "inner_auc_curve": inner_auc_curve,
        "best_n": best_n,
        "method": method,
    }


def _aggregate_panel_summary(
    panel_orders: list[dict[str, object]],
    selection_curves: list[list[dict[str, float | int | str]]],
    heldout_curves: list[list[dict[str, float | int | str]]],
    gene_names: list[str],
    *,
    method: str,
) -> dict[str, object] | None:
    if not panel_orders:
        return None

    max_k = max(len(curve) for curve in heldout_curves)
    heldout_by_k = {k: [] for k in range(1, max_k + 1)}
    inner_by_k = {k: [] for k in range(1, max_k + 1)}
    gene_added_by_k: dict[int, list[str]] = {k: [] for k in range(1, max_k + 1)}

    for curve in heldout_curves:
        for point in curve:
            k = int(point["n_genes"])
            heldout_by_k[k].append(float(point["auc"]))
            gene_added_by_k[k].append(str(point["gene_added"]))
    for curve in selection_curves:
        for point in curve:
            inner_by_k[int(point["n_genes"])].append(float(point["auc"]))

    inner_curve_summary = []
    best_n = 1
    best_inner_auc = -1.0
    for k in range(1, max_k + 1):
        if not heldout_by_k[k]:
            continue
        mean_inner = float(np.mean(inner_by_k[k])) if inner_by_k[k] else 0.5
        mean_heldout = float(np.mean(heldout_by_k[k]))
        heldout_std = float(np.std(heldout_by_k[k]))
        gene_name = max(set(gene_added_by_k[k]), key=gene_added_by_k[k].count)
        inner_curve_summary.append(
            {
                "n_genes": k,
                "auc": mean_heldout,
                "std": heldout_std,
                "selection_auc": mean_inner,
                "selection_std": float(np.std(inner_by_k[k])) if inner_by_k[k] else 0.0,
                "gene_added": gene_name,
            }
        )
        if mean_inner > best_inner_auc:
            best_inner_auc = mean_inner
            best_n = k

    selected_gene_stats: dict[str, dict[str, float]] = {}
    for panel in panel_orders:
        ordered_idx = list(panel["ordered_idx"])
        for rank, gene_idx in enumerate(ordered_idx[:best_n], start=1):
            gene_name = gene_names[gene_idx]
            stats = selected_gene_stats.setdefault(gene_name, {"count": 0.0, "rank_sum": 0.0})
            stats["count"] += 1.0
            stats["rank_sum"] += float(rank)

    ranked_consensus = sorted(
        selected_gene_stats.items(),
        key=lambda item: (
            -item[1]["count"] / len(panel_orders),
            item[1]["rank_sum"] / item[1]["count"],
            item[0],
        ),
    )
    best_genes = [gene for gene, _ in ranked_consensus[:best_n]]
    selection_frequency = [
        {
            "gene": gene,
            "frequency": round(stats["count"] / len(panel_orders), 4),
            "mean_rank": round(stats["rank_sum"] / stats["count"], 4),
        }
        for gene, stats in ranked_consensus
    ]

    heldout_best = heldout_by_k.get(best_n, [])
    return {
        "best_genes": best_genes,
        "best_auc": round(float(np.mean(heldout_best)) if heldout_best else 0.5, 4),
        "auc_std": round(float(np.std(heldout_best)) if heldout_best else 0.0, 4),
        "n_genes": best_n,
        "auc_curve": inner_curve_summary,
        "method": method,
        "auc_note": "nested_outer_cv",
        "evaluation": "nested_outer_cv",
        "selection_frequency": selection_frequency,
    }


def run_biomarker_analysis(
    expression: np.ndarray,
    gene_names: list[str],
    sample_groups: dict[str, list[int]],
    n_top_genes: int = 20,
    n_estimators: int = 500,
    cv_folds: int = 5,
    model: str = "rf",
    panel_method: str = "forward",
    progress_callback=None,
    cancel_check=None,
) -> dict:
    """Full biomarker discovery pipeline with reviewer-safe nested CV."""
    raise_if_cancelled(cancel_check)

    def _progress(step, pct, msg):
        if progress_callback:
            progress_callback(step, pct, msg)

    model_display = _model_display_name(model)
    sample_indices = []
    labels = []
    group_names = sorted(sample_groups.keys())
    for group in group_names:
        sample_indices.extend(sample_groups[group])
        labels.extend([group] * len(sample_groups[group]))

    X_all = expression[:, sample_indices].T
    all_gene_names = list(gene_names)
    le = LabelEncoder()
    y = le.fit_transform(labels)
    n_all_genes = X_all.shape[1]
    prefilter_n = min(PREFILTER_MAX_GENES, n_all_genes)
    max_panel_genes = min(n_top_genes, MAX_PANEL_GENES)

    _progress("preprocessing", 5, f"Data prepared — using {model_display}")

    _, use_shap_tree, needs_scaling = _build_model(model, n_estimators, len(y))
    min_class_count = int(min(np.bincount(y)))
    n_splits = max(2, min(cv_folds, min_class_count))
    outer_cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)

    warnings = []
    roc_skipped_folds: list[tuple[int, str]] = []
    if n_splits < cv_folds:
        warnings.append(
            f"Requested {cv_folds} CV folds, but the smallest class has only {min_class_count} sample(s). "
            f"Using {n_splits} folds instead."
        )
    shap_fallback_used = False
    shap_fallback_folds: list[int] = []

    fold_importances = np.zeros(n_all_genes, dtype=np.float64)
    fold_shap_abs_sum = np.zeros(n_all_genes, dtype=np.float64)
    fold_shap_counts = np.zeros(n_all_genes, dtype=np.float64)
    fold_accuracies: list[float] = []
    sample_shap_abs = np.zeros((len(y), n_all_genes), dtype=np.float64)

    n_classes = len(le.classes_)
    y_bin_all = label_binarize(y, classes=range(n_classes))
    if n_classes == 2:
        y_bin_all = np.column_stack([1 - y_bin_all.ravel(), y_bin_all.ravel()])
    roc_mean_fpr = np.linspace(0, 1, 100)
    fold_roc_tprs: dict[int, list[np.ndarray]] = {i: [] for i in range(n_classes)}
    fold_roc_aucs: dict[int, list[float]] = {i: [] for i in range(n_classes)}

    panel_orders: list[dict[str, object]] = []
    selection_curves: list[list[dict[str, float | int | str]]] = []
    heldout_curves: list[list[dict[str, float | int | str]]] = []

    _progress("training", 15, f"Starting nested CV ({n_splits} folds) with {model_display}")

    for fold_i, (train_idx, test_idx) in enumerate(outer_cv.split(X_all, y)):
        raise_if_cancelled(cancel_check)
        pct = 15 + int(fold_i / n_splits * 50)
        _progress("training", pct, f"CV fold {fold_i + 1}/{n_splits}")

        X_train_all = X_all[train_idx]
        X_test_all = X_all[test_idx]
        y_train = y[train_idx]
        y_test = y[test_idx]

        train_var = np.nan_to_num(np.var(X_train_all, axis=0), nan=0.0)
        fold_var_idx = np.argsort(train_var)[-prefilter_n:]
        X_train = X_train_all[:, fold_var_idx]
        X_test = X_test_all[:, fold_var_idx]

        clf, _, _ = _build_model(model, n_estimators, len(y_train))
        if needs_scaling:
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)
        clf.fit(X_train, y_train)

        this_fi = _feature_importance_vector(clf)
        fold_importances[fold_var_idx] += this_fi
        fold_accuracies.append(float(np.mean(clf.predict(X_test) == y_test)))

        fold_topn_local = np.argsort(this_fi)[-n_top_genes:]
        try:
            clf_topn, _, _ = _build_model(model, min(100, n_estimators), len(y_train))
            X_tr_topn = X_train_all[:, fold_var_idx[fold_topn_local]]
            X_te_topn = X_test_all[:, fold_var_idx[fold_topn_local]]
            if needs_scaling:
                sc_topn = StandardScaler()
                X_tr_topn = sc_topn.fit_transform(X_tr_topn)
                X_te_topn = sc_topn.transform(X_te_topn)
            clf_topn.fit(X_tr_topn, y_train)
            y_prob_topn = _predict_proba_full(clf_topn, X_te_topn, n_classes)
            y_test_bin = y_bin_all[test_idx]
            for class_idx in range(n_classes):
                target = y_test_bin[:, class_idx]
                if target.min() == target.max():
                    continue
                try:
                    fpr_f, tpr_f, _ = roc_curve(target, y_prob_topn[:, class_idx])
                except ValueError:
                    continue
                fold_roc_tprs[class_idx].append(np.interp(roc_mean_fpr, fpr_f, tpr_f))
                fold_roc_tprs[class_idx][-1][0] = 0.0
                fold_roc_aucs[class_idx].append(auc(fpr_f, tpr_f))
        except Exception as exc:
            roc_skipped_folds.append((fold_i + 1, str(exc)))
            _progress("training", pct, f"ROC summary skipped fold {fold_i + 1} ({exc})")

        try:
            import shap

            if use_shap_tree:
                explainer = shap.TreeExplainer(clf)
            else:
                try:
                    explainer = shap.LinearExplainer(clf, X_train)
                except Exception:
                    bg_size = min(50, len(X_train))
                    bg_idx = np.random.RandomState(42).choice(len(X_train), bg_size, replace=False)
                    explainer = shap.KernelExplainer(clf.predict_proba, X_train[bg_idx])

            shap_values_fold = explainer.shap_values(X_test)
            if isinstance(shap_values_fold, list):
                shap_abs_fold = np.mean([np.abs(sv) for sv in shap_values_fold], axis=0)
            elif isinstance(shap_values_fold, np.ndarray) and shap_values_fold.ndim == 3:
                shap_abs_fold = np.mean(np.abs(shap_values_fold), axis=2)
            else:
                shap_abs_fold = np.abs(shap_values_fold)
            if shap_abs_fold.ndim == 1:
                shap_abs_fold = shap_abs_fold.reshape(1, -1)
            fold_shap_abs_sum[fold_var_idx] += np.sum(shap_abs_fold, axis=0)
            fold_shap_counts[fold_var_idx] += shap_abs_fold.shape[0]
            sample_shap_abs[np.ix_(test_idx, fold_var_idx)] = shap_abs_fold
        except Exception as exc:
            shap_fallback_used = True
            shap_fallback_folds.append(fold_i + 1)
            _progress("training", pct, f"SHAP fallback fold {fold_i + 1} ({exc})")
            fi = _feature_importance_vector(clf)
            for test_row in test_idx:
                sample_shap_abs[test_row, fold_var_idx] = fi
            fold_shap_abs_sum[fold_var_idx] += fi * len(test_idx)
            fold_shap_counts[fold_var_idx] += len(test_idx)

        if max_panel_genes > 0:
            panel_candidate_n = min(max(n_top_genes, max_panel_genes), len(fold_var_idx))
            candidate_local = np.argsort(this_fi)[-panel_candidate_n:][::-1]
            candidate_idx = fold_var_idx[candidate_local].tolist()

            _progress("optimal", 68 + int((fold_i + 1) / n_splits * 20), f"Nested panel selection {fold_i + 1}/{n_splits}")
            panel_result = _panel_selection_order(
                X_train_all,
                y_train,
                candidate_idx,
                all_gene_names,
                panel_method=panel_method,
                model=model,
                n_estimators=n_estimators,
                needs_scaling=needs_scaling,
                n_splits=n_splits,
                max_genes=max_panel_genes,
                cancel_check=cancel_check,
            )
            ordered_idx = list(panel_result["ordered_idx"])
            if ordered_idx:
                heldout_curve = _evaluate_prefix_auc_curve(
                    X_train_all,
                    y_train,
                    X_test_all,
                    y_test,
                    ordered_idx,
                    all_gene_names,
                    model=model,
                    n_estimators=n_estimators,
                    needs_scaling=needs_scaling,
                    n_classes=n_classes,
                    cancel_check=cancel_check,
                )
                panel_orders.append(panel_result)
                selection_curves.append(list(panel_result["inner_auc_curve"]))
                heldout_curves.append(heldout_curve)

    avg_importances = fold_importances / n_splits
    avg_shap_mean = fold_shap_abs_sum / np.maximum(fold_shap_counts, 1.0)
    accuracy = float(np.mean(fold_accuracies))

    _progress("shap", 76, "SHAP aggregated across folds")

    top_idx = np.argsort(avg_shap_mean)[-n_top_genes:][::-1]

    _progress("auc", 82, "Aggregating out-of-fold ROC curves")
    roc_data = []
    for class_idx, class_name in enumerate(le.classes_):
        if fold_roc_tprs[class_idx]:
            mean_tpr = np.mean(fold_roc_tprs[class_idx], axis=0)
            mean_tpr[-1] = 1.0
            roc_data.append(
                {
                    "group": str(class_name),
                    "fpr": roc_mean_fpr.tolist(),
                    "tpr": mean_tpr.tolist(),
                    "auc": float(np.mean(fold_roc_aucs[class_idx])),
                    "std": float(np.std(fold_roc_aucs[class_idx])),
                }
            )
        else:
            roc_data.append(
                {
                    "group": str(class_name),
                    "fpr": roc_mean_fpr.tolist(),
                    "tpr": roc_mean_fpr.tolist(),
                    "auc": 0.5,
                    "std": 0.0,
                }
            )

    _progress("optimal", 92, "Aggregating nested panel performance")
    optimal_combo = _aggregate_panel_summary(
        panel_orders,
        selection_curves,
        heldout_curves,
        all_gene_names,
        method=panel_method,
    )
    _progress("complete", 100, "Analysis complete")

    if shap_fallback_used:
        warnings.append(
            "SHAP could not be computed for one or more CV folds. "
            "Feature importance values were used as a fallback for those folds."
        )
    if roc_skipped_folds:
        skipped = ", ".join(f"fold {fold}: {reason}" for fold, reason in roc_skipped_folds)
        warnings.append(
            "One or more ROC summaries could not be computed from the fold-specific top-gene model. "
            f"Skipped {skipped}."
        )

    top_genes = [
        {
            "rank": rank + 1,
            "symbol": all_gene_names[idx],
            "importance": float(avg_importances[idx]),
            "shap_mean_abs": float(avg_shap_mean[idx]),
        }
        for rank, idx in enumerate(top_idx)
    ]
    shap_plot_data = [
        {
            "gene": all_gene_names[idx],
            "values": sample_shap_abs[:, idx].tolist(),
            "expression": X_all[:, idx].tolist(),
        }
        for idx in top_idx
    ]

    return {
        "top_genes": top_genes,
        "shap_plot_data": shap_plot_data,
        "roc_data": roc_data,
        "roc_evaluation": "out_of_fold",
        "accuracy": accuracy,
        "group_names": group_names,
        "sample_labels": labels,
        "n_samples": len(y),
        "optimal_combo": optimal_combo,
        "model": model_display,
        "panel_method": panel_method,
        "cv_folds_requested": cv_folds,
        "cv_folds_used": n_splits,
        "shap_fallback_used": shap_fallback_used,
        "shap_fallback_folds": shap_fallback_folds,
        "warnings": warnings,
        "panel_size_cap_genes": max_panel_genes,
        "performance_scope": {
            "panel_method": (
                "Compact panel metrics below use nested outer-CV. "
                "Top-gene SHAP ranking and ROC above remain out-of-fold summaries of the selected ML model."
            ),
            "panel_size_cap": (
                f"Compact panel selection currently evaluates at most {MAX_PANEL_GENES} genes "
                "within the chosen top-gene candidate pool."
            ),
        },
    }
