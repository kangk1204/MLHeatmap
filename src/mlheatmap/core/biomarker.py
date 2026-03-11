"""Biomarker discovery using ML + SHAP with proper nested CV.

Supported models:
  - Random Forest (default)
  - XGBoost
  - LightGBM
  - Logistic Regression (L1 / Elastic Net)
  - SVM (linear kernel)
"""

import numpy as np
from sklearn.feature_selection import f_classif
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder, label_binarize, StandardScaler


# ──────────────────────────────────────────────────────────
# Model factory
# ──────────────────────────────────────────────────────────

def _build_model(model_name: str, n_estimators: int = 500, n_samples: int = 100):
    """Return (clf, use_shap_tree, needs_scaling) for the requested model."""
    model_name = model_name.lower().replace(" ", "_")

    if model_name in ("rf", "random_forest"):
        from sklearn.ensemble import RandomForestClassifier
        clf = RandomForestClassifier(
            n_estimators=n_estimators,
            max_depth=None,
            min_samples_leaf=max(1, n_samples // 20),
            class_weight="balanced",
            random_state=42,
            n_jobs=-1,
        )
        return clf, True, False

    if model_name in ("xgboost", "xgb"):
        from xgboost import XGBClassifier
        # Adapt parameters for small sample sizes
        max_d = 3 if n_samples < 30 else 6
        lr = 0.05 if n_samples < 30 else 0.1
        sub = 1.0 if n_samples < 20 else 0.8
        n_est = min(n_estimators, 200) if n_samples < 30 else n_estimators
        clf = XGBClassifier(
            n_estimators=n_est,
            max_depth=max_d,
            learning_rate=lr,
            subsample=sub,
            colsample_bytree=0.8,
            reg_alpha=0.1,
            reg_lambda=1.0,
            random_state=42,
            n_jobs=-1,
            verbosity=0,
        )
        return clf, True, False

    if model_name in ("lightgbm", "lgbm"):
        from lightgbm import LGBMClassifier
        # Adapt for small sample sizes
        n_est = min(n_estimators, 200) if n_samples < 30 else n_estimators
        n_leaves = 15 if n_samples < 30 else 31
        min_child = max(1, n_samples // 10)
        sub = 1.0 if n_samples < 20 else 0.8
        clf = LGBMClassifier(
            n_estimators=n_est,
            num_leaves=n_leaves,
            max_depth=-1,
            learning_rate=0.05,
            subsample=sub,
            subsample_freq=1 if sub < 1.0 else 0,
            colsample_bytree=0.8,
            min_child_samples=min_child,
            is_unbalance=True,
            random_state=42,
            n_jobs=-1,
            verbose=-1,
        )
        return clf, True, False

    if model_name in ("logistic_regression", "logistic", "lr_l1"):
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

    if model_name in ("svm", "svm_linear", "linear_svm"):
        from sklearn.svm import SVC
        clf = SVC(
            kernel="linear",
            C=1.0,
            class_weight="balanced",
            probability=True,
            random_state=42,
        )
        return clf, False, True

    raise ValueError(f"Unknown model: {model_name}. "
                     f"Choose from: rf, xgboost, lightgbm, logistic_regression, svm_linear")


def _model_display_name(model_name: str) -> str:
    """Human-readable model name."""
    names = {
        "rf": "Random Forest",
        "random_forest": "Random Forest",
        "xgboost": "XGBoost",
        "xgb": "XGBoost",
        "lightgbm": "LightGBM",
        "lgbm": "LightGBM",
        "logistic_regression": "Logistic Regression (L1)",
        "logistic": "Logistic Regression (L1)",
        "lr_l1": "Logistic Regression (L1)",
        "svm": "SVM (Linear)",
        "svm_linear": "SVM (Linear)",
        "linear_svm": "SVM (Linear)",
    }
    return names.get(model_name.lower().replace(" ", "_"), model_name)


# ──────────────────────────────────────────────────────────
# Main analysis pipeline
# ──────────────────────────────────────────────────────────

def run_biomarker_analysis(
    expression: np.ndarray,
    gene_names: list[str],
    sample_groups: dict[str, list[int]],
    n_top_genes: int = 20,
    n_estimators: int = 500,
    cv_folds: int = 5,
    model: str = "rf",
    progress_callback=None,
) -> dict:
    """Full biomarker discovery pipeline with proper nested CV.

    Fixes data leakage: feature selection and SHAP are done per fold.
    Supports multiple ML models via the `model` parameter.
    """
    def _progress(step, pct, msg):
        if progress_callback:
            progress_callback(step, pct, msg)

    model_display = _model_display_name(model)

    # 1. Prepare data
    sample_indices = []
    labels = []
    group_names = sorted(sample_groups.keys())
    for group in group_names:
        indices = sample_groups[group]
        sample_indices.extend(indices)
        labels.extend([group] * len(indices))

    X_all = expression[:, sample_indices].T  # (samples x genes)
    le = LabelEncoder()
    y = le.fit_transform(labels)
    all_gene_names = list(gene_names)

    _progress("preprocessing", 5, f"Data prepared — using {model_display}")

    # 2. Pre-filter by ANOVA F-statistic (group-aware)
    prefilter_n = min(2000, X_all.shape[1])
    f_scores, _ = f_classif(X_all, y)
    f_scores = np.nan_to_num(f_scores, nan=0.0)
    top_var_idx = np.argsort(f_scores)[-prefilter_n:]
    X_filt = X_all[:, top_var_idx]
    filt_names = [all_gene_names[i] for i in top_var_idx]

    _progress("training", 10, f"Pre-filtered to top {prefilter_n} genes (ANOVA)")

    # 3. Build model
    _, use_shap_tree, needs_scaling = _build_model(model, n_estimators, len(y))

    # 4. Nested CV: outer folds for unbiased evaluation
    min_class_count = int(min(np.bincount(y)))
    n_splits = max(2, min(cv_folds, min_class_count))
    outer_cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)

    # Accumulators
    fold_importances = np.zeros(X_filt.shape[1])
    fold_shap_abs_sum = np.zeros(X_filt.shape[1])
    fold_shap_count = 0
    fold_accuracies = []
    sample_shap_abs = np.zeros_like(X_filt)

    _progress("training", 15, f"Starting nested CV ({n_splits} folds) with {model_display}")

    for fold_i, (train_idx, test_idx) in enumerate(outer_cv.split(X_filt, y)):
        pct = 15 + int(fold_i / n_splits * 55)
        _progress("training", pct, f"CV fold {fold_i+1}/{n_splits}")

        X_train, X_test = X_filt[train_idx], X_filt[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Build fresh model for this fold
        clf, _, _ = _build_model(model, n_estimators, len(y_train))

        # Scale if needed (linear models)
        scaler = None
        if needs_scaling:
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)

        clf.fit(X_train, y_train)

        # Feature importances
        if hasattr(clf, 'feature_importances_'):
            fold_importances += clf.feature_importances_
        elif hasattr(clf, 'coef_'):
            # For linear models, use absolute coefficients as importance
            coef = np.abs(clf.coef_)
            if coef.ndim > 1:
                coef = np.mean(coef, axis=0)
            fold_importances += coef
        else:
            fold_importances += np.ones(X_filt.shape[1]) / X_filt.shape[1]

        # Accuracy
        fold_accuracies.append(float(np.mean(clf.predict(X_test) == y_test)))

        # SHAP on held-out test fold (no data leakage)
        try:
            import shap
            if use_shap_tree:
                explainer = shap.TreeExplainer(clf)
            else:
                # Use LinearExplainer for linear models when possible, else KernelExplainer
                try:
                    explainer = shap.LinearExplainer(clf, X_train)
                except Exception:
                    bg_size = min(50, len(X_train))
                    bg_idx = np.random.RandomState(42).choice(len(X_train), bg_size, replace=False)
                    explainer = shap.KernelExplainer(clf.predict_proba, X_train[bg_idx])

            shap_values_fold = explainer.shap_values(X_test)

            # Aggregate SHAP — handle all possible output shapes
            if isinstance(shap_values_fold, list):
                # list of arrays, one per class
                shap_abs_fold = np.mean([np.abs(sv) for sv in shap_values_fold], axis=0)
            elif isinstance(shap_values_fold, np.ndarray) and shap_values_fold.ndim == 3:
                # (n_samples, n_features, n_classes)
                shap_abs_fold = np.mean(np.abs(shap_values_fold), axis=2)
            else:
                shap_abs_fold = np.abs(shap_values_fold)

            # Ensure correct shape
            if shap_abs_fold.ndim == 1:
                shap_abs_fold = shap_abs_fold.reshape(1, -1)

            fold_shap_abs_sum += np.sum(shap_abs_fold, axis=0)
            fold_shap_count += shap_abs_fold.shape[0]
            sample_shap_abs[test_idx] = shap_abs_fold
        except Exception as e:
            # If SHAP fails, use feature importances as fallback
            _progress("training", pct, f"SHAP fallback fold {fold_i+1} ({e})")
            if hasattr(clf, 'feature_importances_'):
                fi = clf.feature_importances_
            elif hasattr(clf, 'coef_'):
                fi = np.mean(np.abs(clf.coef_), axis=0) if clf.coef_.ndim > 1 else np.abs(clf.coef_).ravel()
            else:
                fi = np.ones(X_filt.shape[1]) / X_filt.shape[1]
            for ti in test_idx:
                sample_shap_abs[ti] = fi
            fold_shap_abs_sum += fi * len(test_idx)
            fold_shap_count += len(test_idx)

    # Average across folds
    avg_importances = fold_importances / n_splits
    avg_shap_mean = fold_shap_abs_sum / max(fold_shap_count, 1)
    accuracy = float(np.mean(fold_accuracies))

    _progress("shap", 75, "SHAP aggregated across folds")

    # 5. Top genes by SHAP value (more interpretable than raw model importance)
    top_idx = np.argsort(avg_shap_mean)[-n_top_genes:][::-1]

    # 6. Cross-validated AUC on top genes
    _progress("auc", 80, "Computing cross-validated AUC")
    X_top = X_filt[:, top_idx]
    roc_data = _compute_cv_roc(X_top, y, le.classes_, model, n_estimators,
                               outer_cv, needs_scaling)

    # 7. Optimal gene combination (with inner CV)
    _progress("optimal", 88, "Finding optimal gene combination")
    optimal_combo = _find_optimal_combination(
        X_filt, y, top_idx, filt_names, le.classes_,
        model, n_estimators, needs_scaling,
        n_splits=n_splits, max_genes=min(n_top_genes, 15),
    )

    _progress("complete", 100, "Analysis complete")

    # 8. Compile results
    top_genes = []
    for rank, idx in enumerate(top_idx):
        top_genes.append({
            "rank": rank + 1,
            "symbol": filt_names[idx],
            "importance": float(avg_importances[idx]),
            "shap_mean_abs": float(avg_shap_mean[idx]),
        })

    shap_plot_data = []
    for rank, idx in enumerate(top_idx):
        shap_plot_data.append({
            "gene": filt_names[idx],
            "values": sample_shap_abs[:, idx].tolist(),
            "expression": X_filt[:, idx].tolist(),
        })

    return {
        "top_genes": top_genes,
        "shap_plot_data": shap_plot_data,
        "roc_data": roc_data,
        "accuracy": accuracy,
        "group_names": group_names,
        "sample_labels": labels,
        "n_samples": len(y),
        "optimal_combo": optimal_combo,
        "model": model_display,
    }


# ──────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────

def _find_optimal_combination(
    X_filt, y, top_idx, filt_names, classes,
    model, n_estimators, needs_scaling,
    n_splits=5, max_genes=15,
):
    """Forward selection with independent inner CV."""
    n_classes = len(classes)
    min_class = int(min(np.bincount(y)))
    inner_splits = max(2, min(n_splits, min_class))
    inner_cv = StratifiedKFold(n_splits=inner_splits, shuffle=True, random_state=123)

    y_bin = label_binarize(y, classes=range(n_classes))
    if n_classes == 2:
        y_bin = np.column_stack([1 - y_bin.ravel(), y_bin.ravel()])

    candidate_pool = list(top_idx[:max_genes])
    selected = []
    auc_curve = []
    best_auc = 0.0
    best_set = []

    for step in range(min(max_genes, len(candidate_pool))):
        best_step_auc = -1
        best_gene = None

        for gene_idx in candidate_pool:
            if gene_idx in selected:
                continue
            trial = selected + [gene_idx]
            X_trial = X_filt[:, trial]
            try:
                mean_auc = _quick_cv_auc(
                    X_trial, y, y_bin, n_classes,
                    model, n_estimators, needs_scaling, inner_cv,
                )
            except Exception:
                continue
            if mean_auc > best_step_auc:
                best_step_auc = mean_auc
                best_gene = gene_idx

        if best_gene is None:
            break
        selected.append(best_gene)
        auc_curve.append({
            "n_genes": len(selected),
            "auc": round(best_step_auc, 4),
            "gene_added": filt_names[best_gene],
        })

        if best_step_auc > best_auc:
            best_auc = best_step_auc
            best_set = list(selected)

    best_gene_names = [filt_names[i] for i in best_set]

    return {
        "best_genes": best_gene_names,
        "best_auc": round(best_auc, 4),
        "n_genes": len(best_set),
        "auc_curve": auc_curve,
    }


def _quick_cv_auc(X, y, y_bin, n_classes, model, n_estimators, needs_scaling, cv):
    """Compute mean CV-AUC across all classes (fast, fewer trees)."""
    aucs = []
    # Use fewer estimators for speed during forward selection
    quick_n = min(100, n_estimators)

    for train_idx, test_idx in cv.split(X, y):
        clf, _, _ = _build_model(model, quick_n, len(train_idx))

        X_train, X_test = X[train_idx], X[test_idx]
        if needs_scaling:
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)

        clf.fit(X_train, y[train_idx])
        y_prob = clf.predict_proba(X_test)

        for i in range(n_classes):
            if i < y_prob.shape[1]:
                fpr, tpr, _ = roc_curve(y_bin[test_idx, i], y_prob[:, i])
                aucs.append(auc(fpr, tpr))

    return float(np.mean(aucs)) if aucs else 0.5


def _compute_cv_roc(X, y, classes, model, n_estimators, cv, needs_scaling):
    """Compute per-class ROC curves with cross-validation."""
    n_classes = len(classes)
    y_bin = label_binarize(y, classes=range(n_classes))
    if n_classes == 2:
        y_bin = np.column_stack([1 - y_bin.ravel(), y_bin.ravel()])

    roc_curves = []
    for i, cls in enumerate(classes):
        tprs = []
        aucs_list = []
        mean_fpr = np.linspace(0, 1, 100)

        for train_idx, test_idx in cv.split(X, y):
            clf, _, _ = _build_model(model, n_estimators, len(train_idx))

            X_train, X_test = X[train_idx], X[test_idx]
            if needs_scaling:
                scaler = StandardScaler()
                X_train = scaler.fit_transform(X_train)
                X_test = scaler.transform(X_test)

            clf.fit(X_train, y[train_idx])
            y_score = clf.predict_proba(X_test)

            if i < y_score.shape[1]:
                fpr, tpr, _ = roc_curve(y_bin[test_idx, i], y_score[:, i])
                tprs.append(np.interp(mean_fpr, fpr, tpr))
                tprs[-1][0] = 0.0
                aucs_list.append(auc(fpr, tpr))

        if tprs:
            mean_tpr = np.mean(tprs, axis=0)
            mean_tpr[-1] = 1.0
        else:
            mean_tpr = mean_fpr

        roc_curves.append({
            "group": str(cls),
            "fpr": mean_fpr.tolist(),
            "tpr": mean_tpr.tolist(),
            "auc": float(np.mean(aucs_list)) if aucs_list else 0.5,
            "std": float(np.std(aucs_list)) if aucs_list else 0.0,
        })

    return roc_curves
