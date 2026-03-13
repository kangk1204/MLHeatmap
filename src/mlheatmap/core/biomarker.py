"""Biomarker discovery using ML + SHAP with proper nested CV.

Supported models:
  - Random Forest (default)
  - XGBoost
  - LightGBM
  - Logistic Regression (L1 / Elastic Net)
  - SVM (linear kernel)
"""

import numpy as np
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder, label_binarize, StandardScaler

from mlheatmap.core.capabilities import normalize_model_name


# ──────────────────────────────────────────────────────────
# Model factory
# ──────────────────────────────────────────────────────────

def _build_model(model_name: str, n_estimators: int = 500, n_samples: int = 100):
    """Return (clf, use_shap_tree, needs_scaling) for the requested model."""
    model_name = normalize_model_name(model_name)

    if model_name in ("rf", "random_forest"):
        from sklearn.ensemble import RandomForestClassifier
        # Constrain tree depth for very small samples to reduce overfitting
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
        # Aggressively constrain for small sample sizes
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
        # Aggressively constrain for small sample sizes to prevent overfitting
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

    raise ValueError(f"Unknown model: {model_name}. "
                     f"Choose from: rf, xgboost, lightgbm, logistic_regression, svm_linear")


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
    panel_method: str = "forward",
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

    # 2. Variance prefilter is performed per fold (training data only)
    #    to prevent any test-fold information influencing the feature space.
    prefilter_n = min(2000, X_all.shape[1])
    n_all_genes = X_all.shape[1]

    _progress("training", 10, f"Preparing nested CV (prefilter top {prefilter_n} per fold)")

    # 3. Build model
    _, use_shap_tree, needs_scaling = _build_model(model, n_estimators, len(y))

    # 4. Nested CV: outer folds for unbiased evaluation
    min_class_count = int(min(np.bincount(y)))
    n_splits = max(2, min(cv_folds, min_class_count))
    outer_cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)

    # Accumulators (on FULL gene space — fold results mapped back)
    fold_importances = np.zeros(n_all_genes)
    fold_shap_abs_sum = np.zeros(n_all_genes)
    fold_shap_count = 0
    fold_accuracies = []
    sample_shap_abs = np.zeros((len(y), n_all_genes))

    # ── Out-of-fold ROC accumulators (truly unbiased) ──
    n_classes = len(le.classes_)
    y_bin_all = label_binarize(y, classes=range(n_classes))
    if n_classes == 2:
        y_bin_all = np.column_stack([1 - y_bin_all.ravel(), y_bin_all.ravel()])
    roc_mean_fpr = np.linspace(0, 1, 100)
    fold_roc_tprs: dict[int, list] = {i: [] for i in range(n_classes)}
    fold_roc_aucs: dict[int, list] = {i: [] for i in range(n_classes)}

    _progress("training", 15, f"Starting nested CV ({n_splits} folds) with {model_display}")

    for fold_i, (train_idx, test_idx) in enumerate(outer_cv.split(X_all, y)):
        pct = 15 + int(fold_i / n_splits * 55)
        _progress("training", pct, f"CV fold {fold_i+1}/{n_splits}")

        y_train, y_test = y[train_idx], y[test_idx]

        # Per-fold variance filter (training data only — no test leakage)
        train_var = np.var(X_all[train_idx], axis=0)
        train_var = np.nan_to_num(train_var, nan=0.0)
        fold_var_idx = np.argsort(train_var)[-prefilter_n:]

        X_train = X_all[train_idx][:, fold_var_idx]
        X_test = X_all[test_idx][:, fold_var_idx]

        # Build fresh model for this fold
        clf, _, _ = _build_model(model, n_estimators, len(y_train))

        # Scale if needed (linear models)
        scaler = None
        if needs_scaling:
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)

        clf.fit(X_train, y_train)

        # Feature importances (captured per-fold for nested top-N evaluation)
        if hasattr(clf, 'feature_importances_'):
            this_fi = clf.feature_importances_
        elif hasattr(clf, 'coef_'):
            coef = np.abs(clf.coef_)
            this_fi = np.mean(coef, axis=0) if coef.ndim > 1 else coef.ravel()
        else:
            this_fi = np.ones(len(fold_var_idx)) / len(fold_var_idx)
        fold_importances[fold_var_idx] += this_fi

        # Accuracy
        fold_accuracies.append(float(np.mean(clf.predict(X_test) == y_test)))

        # ── Per-fold top-N evaluation for unbiased ROC ──
        # Gene selection uses only training-fold importance → no test leakage
        fold_topn_idx = np.argsort(this_fi)[-n_top_genes:]
        try:
            clf_topn, _, _ = _build_model(model, min(100, n_estimators), len(y_train))
            X_tr_topn = X_all[train_idx][:, fold_var_idx[fold_topn_idx]]
            X_te_topn = X_all[test_idx][:, fold_var_idx[fold_topn_idx]]
            if needs_scaling:
                sc_topn = StandardScaler()
                X_tr_topn = sc_topn.fit_transform(X_tr_topn)
                X_te_topn = sc_topn.transform(X_te_topn)
            clf_topn.fit(X_tr_topn, y_train)
            y_prob_topn = clf_topn.predict_proba(X_te_topn)
            y_test_bin = y_bin_all[test_idx]
            for ci in range(n_classes):
                if ci < y_prob_topn.shape[1]:
                    fpr_f, tpr_f, _ = roc_curve(y_test_bin[:, ci], y_prob_topn[:, ci])
                    fold_roc_tprs[ci].append(np.interp(roc_mean_fpr, fpr_f, tpr_f))
                    fold_roc_tprs[ci][-1][0] = 0.0
                    fold_roc_aucs[ci].append(auc(fpr_f, tpr_f))
        except Exception:
            pass  # If top-N model fails for this fold, skip

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

            fold_shap_abs_sum[fold_var_idx] += np.sum(shap_abs_fold, axis=0)
            fold_shap_count += shap_abs_fold.shape[0]
            sample_shap_abs[np.ix_(test_idx, fold_var_idx)] = shap_abs_fold
        except Exception as e:
            # If SHAP fails, use feature importances as fallback
            _progress("training", pct, f"SHAP fallback fold {fold_i+1} ({e})")
            if hasattr(clf, 'feature_importances_'):
                fi = clf.feature_importances_
            elif hasattr(clf, 'coef_'):
                fi = np.mean(np.abs(clf.coef_), axis=0) if clf.coef_.ndim > 1 else np.abs(clf.coef_).ravel()
            else:
                fi = np.ones(len(fold_var_idx)) / len(fold_var_idx)
            for ti in test_idx:
                sample_shap_abs[ti, fold_var_idx] = fi
            fold_shap_abs_sum[fold_var_idx] += fi * len(test_idx)
            fold_shap_count += len(test_idx)

    # Average across folds
    avg_importances = fold_importances / n_splits
    avg_shap_mean = fold_shap_abs_sum / max(fold_shap_count, 1)
    accuracy = float(np.mean(fold_accuracies))

    _progress("shap", 75, "SHAP aggregated across folds")

    # 5. Top genes by SHAP value (more interpretable than raw model importance)
    top_idx = np.argsort(avg_shap_mean)[-n_top_genes:][::-1]

    # 6. Aggregate out-of-fold ROC (unbiased — each fold selected genes
    #    from training-only importance and evaluated on the held-out test set)
    _progress("auc", 80, "Aggregating out-of-fold ROC curves")
    roc_data = []
    for ci, cls in enumerate(le.classes_):
        if fold_roc_tprs[ci]:
            mean_tpr = np.mean(fold_roc_tprs[ci], axis=0)
            mean_tpr[-1] = 1.0
            roc_data.append({
                "group": str(cls),
                "fpr": roc_mean_fpr.tolist(),
                "tpr": mean_tpr.tolist(),
                "auc": float(np.mean(fold_roc_aucs[ci])),
                "std": float(np.std(fold_roc_aucs[ci])),
            })
        else:
            roc_data.append({
                "group": str(cls),
                "fpr": roc_mean_fpr.tolist(),
                "tpr": roc_mean_fpr.tolist(),
                "auc": 0.5,
                "std": 0.0,
            })

    # 7. Optimal gene combination (dispatched by panel_method)
    max_panel = min(n_top_genes, 15)
    combo_args = dict(
        X_filt=X_all, y=y, top_idx=top_idx, filt_names=all_gene_names,
        classes=le.classes_, model=model, n_estimators=n_estimators,
        needs_scaling=needs_scaling, n_splits=n_splits, max_genes=max_panel,
    )

    if panel_method == "lasso":
        _progress("optimal", 85, "Running LASSO panel selection")
        optimal_combo = _lasso_panel(**combo_args, progress_fn=_progress)
    elif panel_method == "stability":
        _progress("optimal", 85, "Running Stability Selection")
        optimal_combo = _stability_panel(**combo_args, progress_fn=_progress)
    elif panel_method == "mrmr":
        _progress("optimal", 85, "Running mRMR panel selection")
        optimal_combo = _mrmr_panel(**combo_args, progress_fn=_progress)
    else:
        _progress("optimal", 88, "Finding optimal gene combination")
        optimal_combo = _find_optimal_combination(**combo_args)

    _progress("complete", 100, "Analysis complete")

    # 8. Compile results
    top_genes = []
    for rank, idx in enumerate(top_idx):
        top_genes.append({
            "rank": rank + 1,
            "symbol": all_gene_names[idx],
            "importance": float(avg_importances[idx]),
            "shap_mean_abs": float(avg_shap_mean[idx]),
        })

    shap_plot_data = []
    for rank, idx in enumerate(top_idx):
        shap_plot_data.append({
            "gene": all_gene_names[idx],
            "values": sample_shap_abs[:, idx].tolist(),
            "expression": X_all[:, idx].tolist(),
        })

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
        "method": "forward",
        "auc_note": "cv_model_selection",
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


def _build_auc_curve(X_filt, y, ordered_idx, filt_names, classes,
                     model, n_estimators, needs_scaling, n_splits):
    """Build step-by-step AUC curve for an ordered list of gene indices."""
    n_classes = len(classes)
    min_class = int(min(np.bincount(y)))
    inner_splits = max(2, min(n_splits, min_class))
    inner_cv = StratifiedKFold(n_splits=inner_splits, shuffle=True, random_state=123)
    y_bin = label_binarize(y, classes=range(n_classes))
    if n_classes == 2:
        y_bin = np.column_stack([1 - y_bin.ravel(), y_bin.ravel()])

    auc_curve = []
    best_auc = 0.0
    best_n = 0
    for k in range(1, len(ordered_idx) + 1):
        X_sub = X_filt[:, ordered_idx[:k]]
        try:
            mean_auc = _quick_cv_auc(X_sub, y, y_bin, n_classes,
                                     model, n_estimators, needs_scaling, inner_cv)
        except Exception:
            mean_auc = 0.5
        auc_curve.append({
            "n_genes": k,
            "auc": round(mean_auc, 4),
            "gene_added": filt_names[ordered_idx[k - 1]],
        })
        if mean_auc > best_auc:
            best_auc = mean_auc
            best_n = k
    best_genes = [filt_names[ordered_idx[i]] for i in range(best_n)]
    return auc_curve, best_genes, round(best_auc, 4), best_n


# ── LASSO panel selection ──

def _lasso_panel(X_filt, y, top_idx, filt_names, classes,
                 model, n_estimators, needs_scaling,
                 n_splits=5, max_genes=15, progress_fn=None):
    """Select genes via L1-penalized logistic regression."""
    from sklearn.linear_model import LogisticRegression

    candidate_idx = list(top_idx[:max_genes])
    X_cand = X_filt[:, candidate_idx]

    # Scale for LASSO
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_cand)

    n_classes = len(classes)
    multi = "ovr" if n_classes > 2 else "auto"

    # Scan C values to find one that selects a reasonable number of genes
    best_c, best_sel = None, []
    for C_val in np.logspace(-3, 1, 20):
        try:
            lr = LogisticRegression(
                penalty="l1", C=C_val, solver="saga", max_iter=5000,
                multi_class=multi, class_weight="balanced", random_state=42,
            )
            lr.fit(X_scaled, y)
            coef = np.abs(lr.coef_)
            if coef.ndim > 1:
                coef = np.mean(coef, axis=0)
            else:
                coef = coef.ravel()
            nonzero = np.where(coef > 1e-10)[0]
            if 1 <= len(nonzero) <= max_genes:
                if len(nonzero) > len(best_sel):
                    best_c = C_val
                    best_sel = nonzero.tolist()
                    best_coef = coef
        except Exception:
            continue

    if progress_fn:
        progress_fn("optimal", 90, f"LASSO selected {len(best_sel)} genes")

    if not best_sel:
        # Fallback: use all candidates ranked by a moderate LASSO
        try:
            lr = LogisticRegression(
                penalty="l1", C=1.0, solver="saga", max_iter=5000,
                multi_class=multi, class_weight="balanced", random_state=42,
            )
            lr.fit(X_scaled, y)
            coef = np.abs(lr.coef_)
            if coef.ndim > 1:
                coef = np.mean(coef, axis=0)
            else:
                coef = coef.ravel()
            ranked = np.argsort(coef)[::-1][:max_genes]
            ordered_idx = [candidate_idx[i] for i in ranked]
        except Exception:
            ordered_idx = candidate_idx[:max_genes]
    else:
        # Rank selected genes by |coefficient|
        ranked = sorted(best_sel, key=lambda i: best_coef[i], reverse=True)
        ordered_idx = [candidate_idx[i] for i in ranked]

    # Build AUC curve
    auc_curve, best_genes, best_auc, best_n = _build_auc_curve(
        X_filt, y, ordered_idx, filt_names, classes,
        model, n_estimators, needs_scaling, n_splits,
    )

    if progress_fn:
        progress_fn("optimal", 97, "LASSO panel complete")

    return {
        "best_genes": best_genes,
        "best_auc": best_auc,
        "n_genes": best_n,
        "auc_curve": auc_curve,
        "method": "lasso",
        "auc_note": "cv_model_selection",
    }


# ── Stability Selection ──

def _stability_panel(X_filt, y, top_idx, filt_names, classes,
                     model, n_estimators, needs_scaling,
                     n_splits=5, max_genes=15, progress_fn=None):
    """Bootstrap LASSO for stable gene selection."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedShuffleSplit

    candidate_idx = list(top_idx[:max_genes])
    X_cand = X_filt[:, candidate_idx]
    n_cand = len(candidate_idx)
    n_classes = len(classes)
    multi = "ovr" if n_classes > 2 else "auto"

    n_bootstrap = 100
    freq = np.zeros(n_cand)

    splitter = StratifiedShuffleSplit(n_splits=n_bootstrap, test_size=0.2, random_state=42)

    for i, (train_idx, _) in enumerate(splitter.split(X_cand, y)):
        if progress_fn and i % 20 == 0:
            pct = 85 + int(i / n_bootstrap * 10)
            progress_fn("optimal", pct, f"Stability bootstrap {i}/{n_bootstrap}")

        X_sub = X_cand[train_idx]
        y_sub = y[train_idx]

        scaler = StandardScaler()
        X_sub_s = scaler.fit_transform(X_sub)

        try:
            lr = LogisticRegression(
                penalty="l1", C=0.5, solver="saga", max_iter=2000,
                multi_class=multi, class_weight="balanced", random_state=i,
            )
            lr.fit(X_sub_s, y_sub)
            coef = np.abs(lr.coef_)
            if coef.ndim > 1:
                coef = np.mean(coef, axis=0)
            else:
                coef = coef.ravel()
            freq[coef > 1e-10] += 1
        except Exception:
            continue

    freq /= n_bootstrap  # selection frequency [0, 1]

    if progress_fn:
        progress_fn("optimal", 95, "Computing stability ranking")

    # Adaptive threshold: start at 0.7, lower if too few genes selected
    threshold = 0.7
    while threshold >= 0.3:
        stable_mask = freq >= threshold
        if np.sum(stable_mask) >= 2:
            break
        threshold -= 0.1

    if np.sum(freq >= threshold) < 2:
        # Fallback: rank all genes with any selection by frequency
        ranked = np.argsort(freq)[::-1][:max_genes]
    else:
        # Take genes meeting the stability threshold, ranked by frequency
        stable_indices = np.where(freq >= threshold)[0]
        ranked = sorted(stable_indices, key=lambda i: freq[i], reverse=True)[:max_genes]

    ordered_idx = [candidate_idx[i] for i in ranked]

    # Build AUC curve
    auc_curve, best_genes, best_auc, best_n = _build_auc_curve(
        X_filt, y, ordered_idx, filt_names, classes,
        model, n_estimators, needs_scaling, n_splits,
    )

    if progress_fn:
        progress_fn("optimal", 97, "Stability Selection complete")

    return {
        "best_genes": best_genes,
        "best_auc": best_auc,
        "n_genes": best_n,
        "auc_curve": auc_curve,
        "method": "stability",
        "auc_note": "cv_model_selection",
    }


# ── mRMR (minimum Redundancy Maximum Relevance) ──

def _mrmr_panel(X_filt, y, top_idx, filt_names, classes,
                model, n_estimators, needs_scaling,
                n_splits=5, max_genes=15, progress_fn=None):
    """Greedy mRMR feature selection using mutual information."""
    from sklearn.feature_selection import mutual_info_classif, mutual_info_regression

    candidate_idx = list(top_idx[:max_genes])
    X_cand = X_filt[:, candidate_idx]
    n_cand = len(candidate_idx)

    # Precompute relevance: MI(gene, target)
    relevance = mutual_info_classif(X_cand, y, random_state=42)

    if progress_fn:
        progress_fn("optimal", 88, "mRMR: relevance computed")

    # Greedy selection
    selected_local = []  # indices into candidate_idx
    remaining = set(range(n_cand))

    for step in range(min(max_genes, n_cand)):
        if not remaining:
            break

        if step == 0:
            # First gene: highest relevance
            best_idx = max(remaining, key=lambda i: relevance[i])
        else:
            # Score = relevance - mean redundancy with selected set
            best_score = -np.inf
            best_idx = None
            for cand in remaining:
                # Redundancy: mean MI with already-selected genes
                redundancy = 0.0
                for sel in selected_local:
                    mi = mutual_info_regression(
                        X_cand[:, sel].reshape(-1, 1),
                        X_cand[:, cand],
                        random_state=42,
                    )[0]
                    redundancy += mi
                redundancy /= len(selected_local)
                score = relevance[cand] - redundancy
                if score > best_score:
                    best_score = score
                    best_idx = cand

            if best_idx is None:
                break

        selected_local.append(best_idx)
        remaining.discard(best_idx)

        if progress_fn and step % 3 == 0:
            pct = 88 + int(step / max_genes * 9)
            progress_fn("optimal", pct, f"mRMR step {step + 1}/{max_genes}")

    ordered_idx = [candidate_idx[i] for i in selected_local]

    # Build AUC curve
    auc_curve, best_genes, best_auc, best_n = _build_auc_curve(
        X_filt, y, ordered_idx, filt_names, classes,
        model, n_estimators, needs_scaling, n_splits,
    )

    if progress_fn:
        progress_fn("optimal", 97, "mRMR panel complete")

    return {
        "best_genes": best_genes,
        "best_auc": best_auc,
        "n_genes": best_n,
        "auc_curve": auc_curve,
        "method": "mrmr",
        "auc_note": "cv_model_selection",
    }
