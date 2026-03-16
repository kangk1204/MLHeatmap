/* Biomarker Panel Logic */
const Biomarker = {
    _activeTab: 'ml',
    _activeStream: null,
    _mlRunToken: 0,
    _degRunToken: 0,
    _pendingDegRequest: false,
    _degAbortController: null,
    VOLCANO_NS_RASTER_THRESHOLD: 4000,

    _escapeHtml(str) {
        const div = document.createElement('div');
        div.textContent = str ?? '';
        return div.innerHTML;
    },

    _hexToRgba(hex, opacity) {
        const cleaned = (hex || '').replace('#', '');
        if (cleaned.length !== 6) return `rgba(107,114,128,${opacity})`;
        const r = parseInt(cleaned.slice(0, 2), 16);
        const g = parseInt(cleaned.slice(2, 4), 16);
        const b = parseInt(cleaned.slice(4, 6), 16);
        return `rgba(${r}, ${g}, ${b}, ${opacity})`;
    },

    _buildVolcanoBackgroundImage(points, {
        xRange,
        yRange,
        color = '#6b7280',
        opacity = 0.28,
        size = 4,
        adaptiveCap = null,
    }) {
        if (!points || points.length === 0) return null;

        const canvas = document.createElement('canvas');
        canvas.width = 1400;
        canvas.height = 1000;
        const ctx = canvas.getContext('2d');
        if (!ctx) return null;

        const [xMin, xMax] = xRange;
        const [yMin, yMax] = yRange;
        const pad = 12;
        const plotWidth = canvas.width - pad * 2;
        const plotHeight = canvas.height - pad * 2;
        const fill = this._hexToRgba(color, opacity);

        ctx.clearRect(0, 0, canvas.width, canvas.height);
        ctx.fillStyle = fill;
        ctx.strokeStyle = fill;

        const toCanvasX = (value) => {
            const ratio = (value - xMin) / Math.max(xMax - xMin, 1e-9);
            return pad + ratio * plotWidth;
        };
        const toCanvasY = (value) => {
            const ratio = (yMax - value) / Math.max(yMax - yMin, 1e-9);
            return pad + ratio * plotHeight;
        };

        points.forEach((point) => {
            const px = toCanvasX(point.log2fc);
            const py = toCanvasY(point.display_neg_log10_p);
            const isCapped = adaptiveCap !== null && point.neg_log10_p > adaptiveCap;
            const radius = isCapped ? size * 0.8 : size * 0.55;

            if (isCapped) {
                ctx.beginPath();
                ctx.moveTo(px, py - radius);
                ctx.lineTo(px + radius, py);
                ctx.lineTo(px, py + radius);
                ctx.lineTo(px - radius, py);
                ctx.closePath();
                ctx.fill();
            } else {
                ctx.beginPath();
                ctx.arc(px, py, radius, 0, Math.PI * 2);
                ctx.fill();
            }
        });

        return canvas.toDataURL('image/png');
    },

    cancelPending() {
        const hadActiveStream = !!this._activeStream;
        const hadDegRequest = this._pendingDegRequest;
        this._mlRunToken += 1;
        this._degRunToken += 1;
        if (this._activeStream) {
            this._activeStream.close();
            this._activeStream = null;
        }
        if (this._degAbortController) {
            this._degAbortController.abort();
            this._degAbortController = null;
        }
        this._pendingDegRequest = false;
        if ((hadActiveStream || hadDegRequest) && App.state.sessionId) {
            API.cancelSession(App.state.sessionId).catch(() => {});
        }

        const mlBtn = document.getElementById('btn-run-biomarker');
        if (mlBtn) mlBtn.disabled = false;
        const degBtn = document.getElementById('btn-run-deg');
        if (degBtn) degBtn.disabled = false;

        const progressEl = document.getElementById('biomarker-progress');
        if (progressEl) progressEl.classList.add('hidden');
        if (hadDegRequest) App.hideLoading();
    },

    init() {
        document.getElementById('btn-run-biomarker').addEventListener('click', () => this.run());
        document.getElementById('btn-to-export').addEventListener('click', () => {
            App.goToPanel('heatmap');
            Heatmap.renderShapHeatmap();
        });

        // DEG buttons
        document.getElementById('btn-run-deg').addEventListener('click', () => this.runDeg());
        document.getElementById('btn-deg-to-heatmap').addEventListener('click', () => {
            App.goToPanel('heatmap');
            Heatmap.renderDegHeatmap();
        });

        // Tab switching
        document.querySelectorAll('.bio-tab').forEach(tab => {
            tab.addEventListener('click', () => this._switchTab(tab.dataset.tab));
        });

        // P-value type selector — update label dynamically
        const pvalType = document.getElementById('deg-pval-type');
        if (pvalType) {
            pvalType.addEventListener('change', () => {
                const label = document.getElementById('deg-pval-label');
                if (label) label.textContent = pvalType.value === 'raw' ? 'P-value Threshold' : 'FDR Threshold';
            });
        }

        // Model-specific parameter visibility
        const modelSelect = document.getElementById('model-select');
        if (modelSelect) {
            modelSelect.addEventListener('change', () => this._updateModelParams());
            this._updateModelParams();
        }
    },

    applyCapabilities(capabilities) {
        const modelSelect = document.getElementById('model-select');
        if (!modelSelect || !capabilities || !capabilities.models) return;

        const firstAvailable = Object.keys(capabilities.models).find(
            key => capabilities.models[key].available,
        );

        Array.from(modelSelect.options).forEach(option => {
            const info = capabilities.models[option.value];
            if (!info) return;
            option.disabled = !info.available;
            option.textContent = info.available
                ? info.label
                : `${info.label} (Docker/full)`;
            option.title = info.available ? '' : (info.unavailable_reason || '');
        });

        if (modelSelect.selectedOptions[0] && modelSelect.selectedOptions[0].disabled && firstAvailable) {
            modelSelect.value = firstAvailable;
        }

        this._updateModelParams();
    },

    _switchTab(tab) {
        this._activeTab = tab;
        document.querySelectorAll('.bio-tab').forEach(t => t.classList.toggle('active', t.dataset.tab === tab));
        document.getElementById('ml-controls').classList.toggle('hidden', tab !== 'ml');
        document.getElementById('deg-controls').classList.toggle('hidden', tab !== 'deg');
        // Show/hide results
        document.getElementById('biomarker-results').classList.toggle('hidden', tab !== 'ml' || !App.state.biomarkerResults);
        document.getElementById('deg-results').classList.toggle('hidden', tab !== 'deg' || !App.state.degResults);
        // Refresh reference group dropdown when switching to DEG tab
        if (tab === 'deg') this._populateRefGroupDropdown();
    },

    _populateRefGroupDropdown() {
        const sel = document.getElementById('deg-reference-group');
        if (!sel) return;
        const groups = Object.keys(App.state.groups || {});
        const prev = sel.value;
        sel.innerHTML = '';
        groups.forEach(g => {
            const opt = document.createElement('option');
            opt.value = g;
            opt.textContent = g;
            sel.appendChild(opt);
        });
        // Restore previous selection if still valid
        if (prev && groups.includes(prev)) {
            sel.value = prev;
        } else if (groups.length > 1) {
            // Default: second group as reference (common convention: control is often second)
            sel.value = groups[1];
        }
    },

    _updateModelParams() {
        const model = document.getElementById('model-select').value;
        const estimatorsGroup = document.getElementById('estimators-group');
        const estimatorsLabel = document.getElementById('estimators-label');

        if (estimatorsGroup) {
            // Tree-based models show estimators, others hide
            const isTreeBased = ['rf', 'xgboost', 'lightgbm'].includes(model);
            estimatorsGroup.style.display = isTreeBased ? '' : 'none';
            if (estimatorsLabel) {
                estimatorsLabel.textContent = model === 'rf' ? 'Trees' :
                    model === 'xgboost' ? 'Boosting Rounds' :
                    model === 'lightgbm' ? 'Boosting Rounds' : 'Trees';
            }
        }

        // Show sample-size warning for boosting models
        const warningEl = document.getElementById('model-warning');
        if (warningEl) {
            const nSamples = (App.state.sampleNames || []).length;
            const isBoosting = ['xgboost', 'lightgbm'].includes(model);
            if (isBoosting && nSamples > 0 && nSamples < 30) {
                warningEl.textContent = `Warning: ${model === 'lightgbm' ? 'LightGBM' : 'XGBoost'} may perform poorly with only ${nSamples} samples. Boosting models need 30+ samples for reliable results. Consider Random Forest or Logistic Regression (L1) instead.`;
                warningEl.classList.remove('hidden');
            } else {
                warningEl.classList.add('hidden');
            }
        }
    },

    run() {
        if (!App.state.sessionId) return App.showToast('No data loaded', 'error');

        const groups = App.state.groups || {};
        if (Object.keys(groups).length < 2) {
            App.showToast('Need at least 2 groups for biomarker analysis', 'error');
            return;
        }

        // Validate inputs before calling API
        const nTopGenes = parseInt(document.getElementById('n-top-genes').value);
        const nEstimators = parseInt(document.getElementById('n-estimators').value);
        const cvFolds = parseInt(document.getElementById('cv-folds').value);

        if (isNaN(nTopGenes) || nTopGenes < 5 || nTopGenes > 200) {
            App.showToast('Top Genes must be between 5 and 200', 'error');
            return;
        }
        if (isNaN(nEstimators) || nEstimators < 50 || nEstimators > 2000) {
            App.showToast('Trees/Rounds must be between 50 and 2000', 'error');
            return;
        }
        if (isNaN(cvFolds) || cvFolds < 2 || cvFolds > 10) {
            App.showToast('CV Folds must be between 2 and 10', 'error');
            return;
        }

        const btn = document.getElementById('btn-run-biomarker');
        btn.disabled = true;

        const progressEl = document.getElementById('biomarker-progress');
        const progressFill = document.getElementById('progress-fill');
        const progressText = document.getElementById('progress-text');
        const resultsEl = document.getElementById('biomarker-results');

        progressEl.classList.remove('hidden');
        resultsEl.classList.add('hidden');
        progressFill.style.width = '0%';

        const modelSelect = document.getElementById('model-select');
        const panelSelect = document.getElementById('panel-method-select');
        const runToken = ++this._mlRunToken;
        if (this._activeStream) this._activeStream.close();
        const es = API.biomarkerStream(App.state.sessionId, {
            nTopGenes,
            nEstimators,
            cvFolds,
            model: modelSelect ? modelSelect.value : 'rf',
            panelMethod: panelSelect ? panelSelect.value : 'forward',
        });
        this._activeStream = es;

        let gotResponse = false;
        const finishWithError = (message) => {
            if (runToken !== this._mlRunToken) return;
            gotResponse = true;
            es.close();
            this._activeStream = null;
            btn.disabled = false;
            progressEl.classList.add('hidden');
            App.showToast(message || 'Analysis failed', 'error');
        };

        es.addEventListener('progress', (e) => {
            if (runToken !== this._mlRunToken) return;
            gotResponse = true;
            const data = JSON.parse(e.data);
            progressFill.style.width = data.pct + '%';
            progressText.textContent = data.msg;
        });

        es.addEventListener('complete', (e) => {
            if (runToken !== this._mlRunToken) return;
            gotResponse = true;
            es.close();
            this._activeStream = null;
            btn.disabled = false;
            progressEl.classList.add('hidden');

            const data = JSON.parse(e.data);
            App.state.biomarkerResults = data;
            this.showResults(data);
            App.markStepCompleted('biomarker');
            App.showToast('Analysis complete — review results below', 'success');
        });

        es.addEventListener('app_error', (e) => {
            if (runToken !== this._mlRunToken) return;
            try {
                const data = JSON.parse(e.data);
                finishWithError(data.detail || 'Analysis failed');
            } catch {
                finishWithError('Analysis failed');
            }
        });

        es.onerror = () => {
            if (runToken !== this._mlRunToken) return;
            es.close();
            this._activeStream = null;
            btn.disabled = false;
            progressEl.classList.add('hidden');
            if (!gotResponse) {
                // Never received any SSE event — likely a 4xx/5xx HTTP error
                App.showToast('Analysis failed — check parameters and try again', 'error');
            }
        };
    },

    showResults(data) {
        document.getElementById('biomarker-results').classList.remove('hidden');

        const warningBox = document.getElementById('biomarker-analysis-warnings');
        if (warningBox) {
            const warnings = Array.isArray(data.warnings) ? data.warnings.filter(Boolean) : [];
            if (warnings.length > 0) {
                warningBox.innerHTML = warnings
                    .map(message => `<p>${this._escapeHtml(message)}</p>`)
                    .join('');
                warningBox.classList.remove('hidden');
            } else {
                warningBox.innerHTML = '';
                warningBox.classList.add('hidden');
            }
        }

        // Accuracy badge
        const acc = data.accuracy;
        const badge = document.getElementById('accuracy-badge');
        const color = acc >= 0.9 ? '#10b981' : acc >= 0.7 ? '#f59e0b' : '#ef4444';
        const modelName = data.model || 'Random Forest';
        badge.textContent = `${modelName} — Nested CV accuracy: ${(acc * 100).toFixed(1)}%`;
        badge.style.background = `${color}20`;
        badge.style.color = color;
        badge.style.border = `1px solid ${color}40`;

        // ROC evaluation label
        const rocTag = document.getElementById('roc-eval-tag');
        if (rocTag) {
            if (data.roc_evaluation === 'out_of_fold') {
                rocTag.textContent = 'Out-of-Fold';
                rocTag.className = 'eval-tag oof';
            } else {
                rocTag.textContent = 'CV Re-eval';
                rocTag.className = 'eval-tag internal';
            }
        }

        const scopeNote = document.getElementById('panel-method-scope-note');
        if (scopeNote) {
            const scope = data.performance_scope || {};
            scopeNote.textContent = scope.panel_method
                || 'Compact panel metrics below use nested outer-CV. SHAP, ROC and accuracy above reflect the selected ML model with out-of-fold evaluation.';
        }

        // SHAP plot
        this.plotShap(data);

        // AUC plot
        this.plotAUC(data);

        // Optimal combination
        if (data.optimal_combo) this.showOptimalCombo(data.optimal_combo);

        // Table
        this.populateTable(data.top_genes);
    },

    plotShap(data) {
        // Sort by SHAP value (ascending for horizontal bar — largest at top)
        const sorted = [...data.top_genes].sort((a, b) => a.shap_mean_abs - b.shap_mean_abs);
        const genes = sorted.map(g => g.symbol);
        const values = sorted.map(g => g.shap_mean_abs);

        const maxVal = Math.max(...values);
        const colors = values.map(v => {
            const ratio = v / maxVal;
            const r = Math.round(59 + ratio * (139 - 59));
            const g = Math.round(130 + ratio * (92 - 130));
            const b = Math.round(246 + ratio * (246 - 246));
            return `rgb(${r},${g},${b})`;
        });

        const trace = {
            x: values,
            y: genes,
            type: 'bar',
            orientation: 'h',
            marker: {
                color: colors,
                line: { width: 0 },
            },
            hovertemplate: '<b>%{y}</b><br>SHAP: %{x:.4f}<extra></extra>',
        };

        const layout = {
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            font: { family: 'Inter, sans-serif', color: '#94a3b8', size: 11 },
            xaxis: {
                title: { text: 'Mean |SHAP value|', font: { size: 12 } },
                gridcolor: 'rgba(255,255,255,0.04)',
                zeroline: false,
            },
            yaxis: {
                gridcolor: 'rgba(255,255,255,0.04)',
                tickfont: { size: 11 },
            },
            margin: { l: 100, r: 20, t: 10, b: 50 },
            bargap: 0.15,
        };

        Plotly.newPlot('shap-plot', [trace], layout, {
            responsive: true,
            displayModeBar: false,
        });
    },

    plotAUC(data) {
        if (!data.roc_data || data.roc_data.length === 0) return;

        const colors = ['#3b82f6', '#ef4444', '#10b981', '#f59e0b', '#8b5cf6'];
        const traces = [];

        data.roc_data.forEach((curve, i) => {
            const color = colors[i % colors.length];

            // Mean ROC
            traces.push({
                x: curve.fpr,
                y: curve.tpr,
                name: `${curve.group} (AUC=${curve.auc.toFixed(3)}±${curve.std.toFixed(3)})`,
                line: { color: color, width: 2.5 },
                mode: 'lines',
            });

            // Confidence band
            const tprUpper = curve.tpr.map(v => Math.min(1, v + curve.std));
            const tprLower = curve.tpr.map(v => Math.max(0, v - curve.std));

            traces.push({
                x: [...curve.fpr, ...curve.fpr.slice().reverse()],
                y: [...tprUpper, ...tprLower.reverse()],
                fill: 'toself',
                fillcolor: color + '15',
                line: { width: 0 },
                showlegend: false,
                hoverinfo: 'skip',
            });
        });

        // Diagonal
        traces.push({
            x: [0, 1],
            y: [0, 1],
            line: { color: 'rgba(255,255,255,0.15)', dash: 'dash', width: 1 },
            showlegend: false,
            hoverinfo: 'skip',
        });

        const layout = {
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            font: { family: 'Inter, sans-serif', color: '#94a3b8', size: 11 },
            xaxis: {
                title: { text: 'False Positive Rate', font: { size: 12 } },
                range: [-0.02, 1.02],
                gridcolor: 'rgba(255,255,255,0.04)',
                zeroline: false,
            },
            yaxis: {
                title: { text: 'True Positive Rate', font: { size: 12 } },
                range: [-0.02, 1.02],
                gridcolor: 'rgba(255,255,255,0.04)',
                zeroline: false,
            },
            legend: {
                x: 0.45,
                y: 0.05,
                bgcolor: 'rgba(0,0,0,0.3)',
                bordercolor: 'rgba(255,255,255,0.1)',
                borderwidth: 1,
                font: { size: 10 },
            },
            margin: { l: 60, r: 20, t: 10, b: 60 },
        };

        Plotly.newPlot('auc-plot', traces, layout, {
            responsive: true,
            displayModeBar: false,
        });
    },

    showOptimalCombo(combo) {
        // Build SHAP rank lookup from top_genes
        const topGenes = App.state.biomarkerResults ? App.state.biomarkerResults.top_genes : [];
        const shapRankMap = {};
        const shapSorted = [...topGenes].sort((a, b) => b.shap_mean_abs - a.shap_mean_abs);
        shapSorted.forEach((g, i) => { shapRankMap[g.symbol] = i + 1; });

        // Method description (dynamic based on method used)
        const descEl = document.getElementById('optimal-method-desc');
        const methodDescs = {
            forward: 'Forward selection runs inside each training fold. Genes are added in the order that maximizes inner-CV AUC on the training split only.',
            lasso: 'LASSO uses sparse logistic models on each training fold, then nested CV reports held-out panel performance separately from fold-local selection.',
            stability: 'Stability Selection bootstraps LASSO inside each training fold. Consensus genes below are ranked by cross-fold selection frequency and mean rank.',
            mrmr: 'mRMR selects genes inside each training fold by balancing relevance and redundancy. The panel summary below reports nested outer-CV held-out AUC.',
        };
        descEl.textContent = methodDescs[combo.method] || methodDescs.forward;

        // Summary badge
        const summary = document.getElementById('optimal-combo-summary');
        const aucPct = (combo.best_auc * 100).toFixed(1);
        const color = combo.best_auc >= 0.9 ? '#10b981' : combo.best_auc >= 0.7 ? '#f59e0b' : '#ef4444';
        summary.textContent = '';
        const badge = document.createElement('div');
        badge.className = 'optimal-badge';
        badge.style.cssText = `background:${color}15;border:1px solid ${color}40;color:${color}`;
        badge.textContent = `Nested outer-CV panel AUC: ${aucPct}% ± ${((combo.auc_std || 0) * 100).toFixed(1)}% with ${combo.n_genes} gene${combo.n_genes > 1 ? 's' : ''}`;
        summary.appendChild(badge);

        // Panel AUC note — distinguish from the out-of-fold ROC above
        const noteEl = document.getElementById('panel-auc-note');
        if (noteEl) {
            if (combo.auc_note === 'nested_outer_cv') {
                noteEl.textContent = 'Panel size is chosen by inner cross-validation on each training fold. '
                    + 'The plotted panel AUC below is the mean held-out AUC from the outer folds, and the consensus genes are aggregated from fold-local panels.';
                noteEl.classList.remove('hidden');
            } else {
                noteEl.classList.add('hidden');
            }
        }

        // AUC vs #genes curve (with SHAP rank in tooltip)
        const curve = combo.auc_curve;
        const trace = {
            x: curve.map(c => c.n_genes),
            y: curve.map(c => c.auc),
            mode: 'lines+markers',
            line: { color: '#8b5cf6', width: 2.5 },
            marker: { size: 8, color: '#8b5cf6', line: { color: '#fff', width: 1 } },
            text: curve.map(c => {
                const sr = shapRankMap[c.gene_added];
                return sr ? `+${c.gene_added} (SHAP #${sr})` : `+${c.gene_added}`;
            }),
            customdata: curve.map(c => [c.selection_auc || null, c.std || 0]),
            hovertemplate: '<b>%{x} genes</b><br>Held-out AUC: %{y:.4f}<br>Selection CV-AUC: %{customdata[0]:.4f}<br>Outer-fold SD: %{customdata[1]:.4f}<br>Added: %{text}<extra></extra>',
        };
        // Highlight best point
        let bestIdx = curve.findIndex(c => c.n_genes === combo.n_genes);
        if (bestIdx < 0) bestIdx = curve.reduce((bi, c, i, arr) => c.auc > arr[bi].auc ? i : bi, 0);
        const bestTrace = {
            x: [curve[bestIdx].n_genes],
            y: [curve[bestIdx].auc],
            mode: 'markers',
            marker: { size: 14, color: '#10b981', symbol: 'star', line: { color: '#fff', width: 1.5 } },
            showlegend: false,
            hoverinfo: 'skip',
        };

        const layout = {
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            font: { family: 'Inter, sans-serif', color: '#94a3b8', size: 11 },
            xaxis: {
                title: { text: 'Number of Genes', font: { size: 12 } },
                gridcolor: 'rgba(255,255,255,0.04)',
                dtick: 1,
            },
            yaxis: {
                title: { text: 'Mean held-out AUC', font: { size: 12 } },
                gridcolor: 'rgba(255,255,255,0.04)',
                zeroline: false,
                range: [Math.max(0, Math.min(...curve.map(c => c.auc)) - 0.05), 1.02],
            },
            margin: { l: 60, r: 20, t: 10, b: 50 },
        };

        Plotly.newPlot('optimal-auc-curve', [trace, bestTrace], layout, {
            responsive: true,
            displayModeBar: false,
        });

        // Gene list chips with SHAP rank
        const listEl = document.getElementById('optimal-gene-list');
        listEl.innerHTML = '<h4 style="margin:0 0 8px;color:var(--text-primary)">Optimal Gene Set:</h4>';
        const chipContainer = document.createElement('div');
        chipContainer.style.cssText = 'display:flex;flex-wrap:wrap;gap:6px;';

        combo.best_genes.forEach((gene, i) => {
            const chip = document.createElement('span');
            chip.className = 'optimal-gene-chip';
            chip.textContent = `${i + 1}. ${gene}`;

            // Show SHAP rank
            const shapRank = shapRankMap[gene];
            if (shapRank) {
                const tag = document.createElement('span');
                tag.className = 'shap-rank-tag';
                tag.textContent = `SHAP #${shapRank}`;
                chip.appendChild(tag);
            }

            const freqInfo = Array.isArray(combo.selection_frequency)
                ? combo.selection_frequency.find(item => item.gene === gene)
                : null;
            if (freqInfo) {
                const tag = document.createElement('span');
                tag.className = 'shap-rank-tag';
                tag.textContent = `${(freqInfo.frequency * 100).toFixed(0)}% folds`;
                chip.appendChild(tag);
            }

            chipContainer.appendChild(chip);
        });
        listEl.appendChild(chipContainer);
    },

    populateTable(topGenes) {
        const tbody = document.querySelector('#biomarker-table tbody');
        tbody.innerHTML = '';

        // Sort by SHAP value (descending) and re-rank
        const sorted = [...topGenes].sort((a, b) => b.shap_mean_abs - a.shap_mean_abs);
        sorted.forEach((gene, i) => { gene.rank = i + 1; });

        sorted.forEach(gene => {
            const tr = document.createElement('tr');

            const tdRank = document.createElement('td');
            tdRank.style.cssText = 'font-weight:600;color:var(--text-accent)';
            tdRank.textContent = gene.rank;

            const tdSymbol = document.createElement('td');
            tdSymbol.style.fontWeight = '600';
            tdSymbol.textContent = gene.symbol;

            const tdImp = document.createElement('td');
            tdImp.textContent = gene.importance.toFixed(4);

            const tdShap = document.createElement('td');
            tdShap.textContent = gene.shap_mean_abs.toFixed(4);

            tr.appendChild(tdRank);
            tr.appendChild(tdSymbol);
            tr.appendChild(tdImp);
            tr.appendChild(tdShap);
            tbody.appendChild(tr);
        });
    },

    // =============================================
    // DEG Analysis
    // =============================================
    async runDeg() {
        if (!App.state.sessionId) return App.showToast('No data loaded', 'error');

        const groups = App.state.groups || {};
        if (Object.keys(groups).length !== 2) {
            App.showToast('DEG requires exactly 2 groups', 'error');
            return;
        }

        const btn = document.getElementById('btn-run-deg');
        btn.disabled = true;
        const runToken = ++this._degRunToken;
        this._pendingDegRequest = true;
        if (this._degAbortController) this._degAbortController.abort();
        const controller = new AbortController();
        this._degAbortController = controller;
        App.showLoading('Running DEG analysis...');

        try {
            const useRaw = document.getElementById('deg-pval-type').value === 'raw';
            const refGroup = document.getElementById('deg-reference-group')?.value || '';
            const data = await API.runDeg(App.state.sessionId, {
                method: document.getElementById('deg-method-select').value,
                log2fcThreshold: parseFloat(document.getElementById('deg-fc-threshold').value),
                pvalueThreshold: parseFloat(document.getElementById('deg-pval-threshold').value),
                useRawPvalue: useRaw,
                referenceGroup: refGroup,
                signal: controller.signal,
            });
            if (runToken !== this._degRunToken) return;

            App.state.degResults = data;
            this.showDegResults(data);
            App.markStepCompleted('biomarker');
            const cg = data.comparison_group || '';
            const toastUp = cg ? `${data.summary.n_up} up in ${cg}` : `${data.summary.n_up} up`;
            App.showToast(`DEG complete: ${toastUp}, ${data.summary.n_down} down`, 'success');
        } catch (err) {
            if (runToken !== this._degRunToken) return;
            if (err.name === 'AbortError') return;
            App.showToast(err.message, 'error');
        } finally {
            if (this._degAbortController === controller) this._degAbortController = null;
            if (runToken !== this._degRunToken) return;
            this._pendingDegRequest = false;
            btn.disabled = false;
            App.hideLoading();
        }
    },

    showDegResults(data) {
        document.getElementById('deg-results').classList.remove('hidden');
        document.getElementById('biomarker-results').classList.add('hidden');

        const isRawP = data.pvalue_type === 'raw';
        const pLabel = isRawP ? 'Raw P-value' : 'FDR';
        const cutoff = data.thresholds.pvalue;
        const comp = data.comparison_group || '';
        const ref = data.reference_group || '';
        const upLabel = comp ? `Up in ${comp}` : 'Up-regulated';
        const downLabel = comp ? `Down in ${comp}` : 'Down-regulated';

        // Summary cards
        const summaryEl = document.getElementById('deg-summary');
        summaryEl.innerHTML = '';
        const cards = [
            { className: 'deg-stat up', count: String(data.summary.n_up), label: upLabel },
            { className: 'deg-stat down', count: String(data.summary.n_down), label: downLabel },
            { className: 'deg-stat ns', count: String(data.summary.n_not_significant), label: 'Not Significant' },
            { className: 'deg-stat info', count: `${pLabel} < ${cutoff}`, label: `vs ${ref} (ref)`, small: true },
        ];
        cards.forEach(card => {
            const wrapper = document.createElement('div');
            wrapper.className = card.className;

            const count = document.createElement('span');
            count.className = 'count';
            if (card.small) count.style.fontSize = '0.85rem';
            count.textContent = card.count;

            const label = document.createElement('span');
            label.className = 'label';
            label.textContent = card.label;

            wrapper.appendChild(count);
            wrapper.appendChild(label);
            summaryEl.appendChild(wrapper);
        });

        const degNote = document.getElementById('deg-results-note');
        if (degNote) {
            const basis = (data.effect_size_basis || 'normalized_expression').replaceAll('_', ' ');
            const norm = data.normalization_method || 'normalized';
            degNote.textContent = `Exploratory DEG on ${norm} data. log2FC is computed from ${basis} on the linear scale; the statistical test remains ${data.method}.`;
            degNote.classList.remove('hidden');
        }

        // Volcano plot
        this.plotVolcano(data);

        // DEG table
        this.populateDegTable(data.table_results || data.results, comp);
    },

    plotVolcano(data) {
        const results = data.plot_results || data.results;
        const fcThresh = data.thresholds.log2fc;
        const pThresh = data.thresholds.pvalue;
        const negLog10PThresh = -Math.log10(pThresh);
        const isRawP = data.pvalue_type === 'raw';
        const pLabel = isRawP ? 'P-value' : 'FDR';
        const comp = data.comparison_group || '';
        const ref = data.reference_group || '';
        const compSafe = this._escapeHtml(comp);
        const refSafe = this._escapeHtml(ref);
        const fcAxisLabel = (compSafe && refSafe) ? `log2FC (${compSafe} / ${refSafe})` : 'log2 Fold Change';

        // Adaptive Y-cap: prevent extreme -log10(p) values from
        // squashing all points at the top of the plot.
        // Strategy: find the "useful" range where most differentiation
        // occurs, then cap above it.  Points above the cap are shown as
        // diamond markers at the cap line (hover still shows true value).
        const allNlpRaw = results.map(r => r.neg_log10_p);
        const maxRawNlp = Math.max(...allNlpRaw);
        const sortedNlp = [...allNlpRaw].sort((a, b) => a - b);
        const n = sortedNlp.length;
        const p50 = sortedNlp[Math.floor(n * 0.50)] || 1;
        const p90 = sortedNlp[Math.floor(n * 0.90)] || 1;
        // Base cap: focus on spreading out the majority of data
        const baseCap = Math.max(negLog10PThresh * 3, p90 * 1.8, 10);
        // Hard scientific ceiling — beyond ~50, "very significant" is enough
        const adaptiveCap = Math.min(baseCap, 50);
        const needsCap = maxRawNlp > adaptiveCap * 1.05;

        // Clamp displayed y-values; keep originals in hover
        const clampY = v => Math.min(v, adaptiveCap);
        const decoratedResults = results.map(r => ({
            ...r,
            display_neg_log10_p: clampY(r.neg_log10_p),
        }));

        // Separate by direction after display coordinates are available
        const up = decoratedResults.filter(r => r.direction === 'up');
        const down = decoratedResults.filter(r => r.direction === 'down');
        const ns = decoratedResults.filter(r => r.direction === 'ns');

        const makeInteractiveTrace = (subset, name, color, size, opacity) => ({
            x: subset.map(r => r.log2fc),
            y: subset.map(r => r.display_neg_log10_p),
            text: subset.map(r => r.gene),
            customdata: subset.map(r => r.neg_log10_p),
            type: 'scattergl',
            mode: 'markers',
            name: `${name} (${subset.length})`,
            marker: {
                color,
                size: subset.map(r => r.neg_log10_p > adaptiveCap ? size + 2 : size),
                opacity,
                symbol: subset.map(r => r.neg_log10_p > adaptiveCap ? 'diamond' : 'circle'),
            },
            hovertemplate: `<b>%{text}</b><br>log2FC: %{x:.3f}<br>-log10(${pLabel}): %{customdata:.2f}<extra></extra>`,
        });

        const makeBackgroundTrace = (subset, name, color, size, opacity) => ({
            x: subset.map(r => r.log2fc),
            y: subset.map(r => r.display_neg_log10_p),
            type: 'scattergl',
            mode: 'markers',
            name: `${name} (${subset.length})`,
            marker: {
                color,
                size,
                opacity,
                line: { width: 0 },
            },
            hoverinfo: 'skip',
        });

        const upName = compSafe ? `Up in ${compSafe}` : 'Up';
        const downName = compSafe ? `Down in ${compSafe}` : 'Down';
        const allFc = decoratedResults.map(r => r.log2fc);
        const maxFc = Math.max(Math.abs(Math.min(...allFc)), Math.abs(Math.max(...allFc)), fcThresh + 1);
        const maxNlp = needsCap ? adaptiveCap : Math.max(...allNlpRaw, negLog10PThresh + 1);
        const yPadding = maxNlp * 0.12;
        const xRange = [-maxFc * 1.1, maxFc * 1.1];
        const yRange = [-0.5, maxNlp + yPadding];
        const shouldRasterizeNs = ns.length >= this.VOLCANO_NS_RASTER_THRESHOLD;

        const traces = [];
        if (shouldRasterizeNs) {
            traces.push({
                x: [null],
                y: [null],
                type: 'scatter',
                mode: 'markers',
                name: `NS (${ns.length})`,
                marker: {
                    color: '#6b7280',
                    size: 6,
                    opacity: 0.28,
                    line: { width: 0 },
                },
                hoverinfo: 'skip',
                showlegend: true,
            });
        } else {
            traces.push(makeBackgroundTrace(ns, 'NS', '#6b7280', 4, 0.28));
        }
        traces.push(makeInteractiveTrace(up, upName, '#ef4444', 7, 0.85));
        traces.push(makeInteractiveTrace(down, downName, '#3b82f6', 7, 0.85));

        // Top gene labels
        const sigGenes = [...up, ...down].sort((a, b) => {
            const pa = isRawP ? a.pvalue : a.adj_pvalue;
            const pb = isRawP ? b.pvalue : b.adj_pvalue;
            return pa - pb;
        }).slice(0, 10);
        if (sigGenes.length > 0) {
            traces.push({
                x: sigGenes.map(r => r.log2fc),
                y: sigGenes.map(r => r.display_neg_log10_p),
                text: sigGenes.map(r => r.gene),
                mode: 'text',
                textposition: 'top center',
                textfont: { size: 9, color: '#e2e8f0' },
                showlegend: false,
                hoverinfo: 'skip',
            });
        }

        const layout = {
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            font: { family: 'Inter, sans-serif', color: '#94a3b8', size: 11 },
            xaxis: {
                title: { text: fcAxisLabel, font: { size: 12 } },
                gridcolor: 'rgba(255,255,255,0.04)',
                zeroline: false,
                range: xRange,
            },
            yaxis: {
                title: { text: `-log10(${pLabel})`, font: { size: 12 } },
                gridcolor: 'rgba(255,255,255,0.04)',
                zeroline: false,
                range: yRange,
            },
            legend: { x: 0.02, y: 0.98, bgcolor: 'rgba(0,0,0,0.3)', font: { size: 10 } },
            margin: { l: 60, r: 20, t: 10, b: 60 },
            images: [],
            shapes: [
                // Vertical FC thresholds
                { type: 'line', x0: fcThresh, x1: fcThresh, y0: 0, y1: 1, yref: 'paper', line: { color: 'rgba(255,255,255,0.15)', dash: 'dash', width: 1 } },
                { type: 'line', x0: -fcThresh, x1: -fcThresh, y0: 0, y1: 1, yref: 'paper', line: { color: 'rgba(255,255,255,0.15)', dash: 'dash', width: 1 } },
                // Horizontal p-value threshold
                { type: 'line', x0: 0, x1: 1, xref: 'paper', y0: negLog10PThresh, y1: negLog10PThresh, line: { color: 'rgba(255,255,255,0.15)', dash: 'dash', width: 1 } },
                // Y-axis cap indicator (only when capping is active)
                ...(needsCap ? [{
                    type: 'line', x0: 0, x1: 1, xref: 'paper',
                    y0: adaptiveCap, y1: adaptiveCap,
                    line: { color: 'rgba(251,191,36,0.3)', dash: 'dot', width: 1 },
                }] : []),
            ],
            annotations: needsCap ? [{
                x: 1, xref: 'paper', y: adaptiveCap, yanchor: 'bottom',
                text: `capped (\u25C6 = higher)`,
                showarrow: false, font: { size: 9, color: 'rgba(251,191,36,0.6)' },
            }] : [],
        };

        if (shouldRasterizeNs) {
            const nsImage = this._buildVolcanoBackgroundImage(ns, {
                xRange,
                yRange,
                color: '#6b7280',
                opacity: 0.28,
                size: 4,
                adaptiveCap,
            });
            if (nsImage) {
                layout.images.push({
                    source: nsImage,
                    xref: 'x',
                    yref: 'y',
                    x: xRange[0],
                    y: yRange[1],
                    sizex: xRange[1] - xRange[0],
                    sizey: yRange[1] - yRange[0],
                    sizing: 'stretch',
                    opacity: 1,
                    layer: 'below',
                });
            }
        }

        Plotly.newPlot('volcano-plot', traces, layout, {
            responsive: true, displayModeBar: false,
        });
    },

    populateDegTable(results, compGroup) {
        const tbody = document.querySelector('#deg-table tbody');
        tbody.innerHTML = '';

        // Show top 50 significant genes
        let sig = results.filter(r => r.direction !== 'ns').slice(0, 50);
        if (sig.length === 0) {
            // Fallback: show top 50 by p-value
            sig.push(...results.slice(0, 50));
        }

        const upText = compGroup ? `▲ ${compGroup}` : '▲ Up';
        const downText = compGroup ? `▼ ${compGroup}` : '▼ Down';
        sig.forEach((gene, i) => {
            const tr = document.createElement('tr');
            const dirText = gene.direction === 'up' ? upText : gene.direction === 'down' ? downText : '— NS';

            const tdRank = document.createElement('td');
            tdRank.style.cssText = 'font-weight:600;color:var(--text-accent)';
            tdRank.textContent = String(i + 1);

            const tdGene = document.createElement('td');
            tdGene.style.fontWeight = '600';
            tdGene.textContent = gene.gene;

            const tdLog2fc = document.createElement('td');
            tdLog2fc.textContent = gene.log2fc.toFixed(3);

            const tdP = document.createElement('td');
            tdP.textContent = gene.pvalue.toExponential(2);

            const tdAdjP = document.createElement('td');
            tdAdjP.textContent = gene.adj_pvalue.toExponential(2);

            const tdDir = document.createElement('td');
            tdDir.className = `dir-${gene.direction}`;
            tdDir.textContent = dirText;

            [tdRank, tdGene, tdLog2fc, tdP, tdAdjP, tdDir].forEach(td => tr.appendChild(td));
            tbody.appendChild(tr);
        });
    },
};
