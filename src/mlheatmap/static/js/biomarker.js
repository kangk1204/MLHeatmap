/* Biomarker Panel Logic */
const Biomarker = {
    _activeTab: 'ml',

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

    _switchTab(tab) {
        this._activeTab = tab;
        document.querySelectorAll('.bio-tab').forEach(t => t.classList.toggle('active', t.dataset.tab === tab));
        document.getElementById('ml-controls').classList.toggle('hidden', tab !== 'ml');
        document.getElementById('deg-controls').classList.toggle('hidden', tab !== 'deg');
        // Show/hide results
        document.getElementById('biomarker-results').classList.toggle('hidden', tab !== 'ml' || !App.state.biomarkerResults);
        document.getElementById('deg-results').classList.toggle('hidden', tab !== 'deg' || !App.state.degResults);
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
        const es = API.biomarkerStream(App.state.sessionId, {
            nTopGenes: parseInt(document.getElementById('n-top-genes').value),
            nEstimators: parseInt(document.getElementById('n-estimators').value),
            cvFolds: parseInt(document.getElementById('cv-folds').value),
            model: modelSelect ? modelSelect.value : 'rf',
            panelMethod: panelSelect ? panelSelect.value : 'forward',
        });

        es.addEventListener('progress', (e) => {
            const data = JSON.parse(e.data);
            progressFill.style.width = data.pct + '%';
            progressText.textContent = data.msg;
        });

        es.addEventListener('complete', (e) => {
            es.close();
            btn.disabled = false;
            progressEl.classList.add('hidden');

            const data = JSON.parse(e.data);
            App.state.biomarkerResults = data;
            this.showResults(data);
            App.markStepCompleted('biomarker');
            App.showToast('Analysis complete — review results below', 'success');
        });

        es.addEventListener('error', (e) => {
            es.close();
            btn.disabled = false;
            progressEl.classList.add('hidden');

            if (e.data) {
                const data = JSON.parse(e.data);
                App.showToast(data.detail || 'Analysis failed', 'error');
            } else {
                App.showToast('Connection lost', 'error');
            }
        });

        es.onerror = () => {
            es.close();
            btn.disabled = false;
            progressEl.classList.add('hidden');
            App.showToast('Analysis failed', 'error');
        };
    },

    showResults(data) {
        document.getElementById('biomarker-results').classList.remove('hidden');

        // Accuracy badge
        const acc = data.accuracy;
        const badge = document.getElementById('accuracy-badge');
        const color = acc >= 0.9 ? '#10b981' : acc >= 0.7 ? '#f59e0b' : '#ef4444';
        const modelName = data.model || 'Random Forest';
        badge.textContent = `${modelName} — CV Accuracy: ${(acc * 100).toFixed(1)}%`;
        badge.style.background = `${color}20`;
        badge.style.color = color;
        badge.style.border = `1px solid ${color}40`;

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
            forward: 'Forward selection from SHAP top candidates: at each step, the gene yielding the highest CV-AUC is added. Selection order may differ from SHAP ranking.',
            lasso: 'LASSO (L1-penalized logistic regression) selects genes with non-zero coefficients, ranked by absolute coefficient magnitude. Naturally produces sparse panels.',
            stability: 'Stability Selection bootstraps LASSO 100× on random 80% subsamples. Genes selected in ≥70% of iterations are deemed stable. Ranked by selection frequency.',
            mrmr: 'mRMR (minimum Redundancy Maximum Relevance) greedily selects genes maximizing mutual information with the target while minimizing redundancy with already-selected genes.',
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
        badge.textContent = `Best AUC: ${aucPct}% with ${combo.n_genes} gene${combo.n_genes > 1 ? 's' : ''}`;
        summary.appendChild(badge);

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
            hovertemplate: '<b>%{x} genes</b><br>AUC: %{y:.4f}<br>Added: %{text}<extra></extra>',
        };
        // Highlight best point
        let bestIdx = curve.findIndex(c => c.auc === combo.best_auc);
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
                title: { text: 'Mean CV-AUC', font: { size: 12 } },
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
        App.showLoading('Running DEG analysis...');

        try {
            const useRaw = document.getElementById('deg-pval-type').value === 'raw';
            const data = await API.runDeg(App.state.sessionId, {
                method: document.getElementById('deg-method-select').value,
                log2fcThreshold: parseFloat(document.getElementById('deg-fc-threshold').value),
                pvalueThreshold: parseFloat(document.getElementById('deg-pval-threshold').value),
                useRawPvalue: useRaw,
            });

            App.state.degResults = data;
            this.showDegResults(data);
            App.markStepCompleted('biomarker');
            App.showToast(`DEG complete: ${data.summary.n_up} up, ${data.summary.n_down} down`, 'success');
        } catch (err) {
            App.showToast(err.message, 'error');
        } finally {
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

        // Summary cards
        const summaryEl = document.getElementById('deg-summary');
        summaryEl.innerHTML = `
            <div class="deg-stat up"><span class="count">${data.summary.n_up}</span><span class="label">Up-regulated</span></div>
            <div class="deg-stat down"><span class="count">${data.summary.n_down}</span><span class="label">Down-regulated</span></div>
            <div class="deg-stat ns"><span class="count">${data.summary.n_not_significant}</span><span class="label">Not Significant</span></div>
            <div class="deg-stat info"><span class="count" style="font-size:0.85rem">${pLabel} < ${cutoff}</span><span class="label">Cutoff</span></div>
        `;

        // Volcano plot
        this.plotVolcano(data);

        // DEG table
        this.populateDegTable(data.results);
    },

    plotVolcano(data) {
        const results = data.results;
        const fcThresh = data.thresholds.log2fc;
        const pThresh = data.thresholds.pvalue;
        const negLog10PThresh = -Math.log10(pThresh);
        const isRawP = data.pvalue_type === 'raw';
        const pLabel = isRawP ? 'P-value' : 'FDR';

        // Separate by direction
        const up = results.filter(r => r.direction === 'up');
        const down = results.filter(r => r.direction === 'down');
        const ns = results.filter(r => r.direction === 'ns');

        const makeTrace = (subset, name, color, size, opacity) => ({
            x: subset.map(r => r.log2fc),
            y: subset.map(r => r.neg_log10_p),
            text: subset.map(r => r.gene),
            mode: 'markers',
            name: `${name} (${subset.length})`,
            marker: { color, size, opacity },
            hovertemplate: `<b>%{text}</b><br>log2FC: %{x:.3f}<br>-log10(${pLabel}): %{y:.2f}<extra></extra>`,
        });

        const traces = [
            makeTrace(ns, 'NS', '#6b7280', 4, 0.45),
            makeTrace(up, 'Up', '#ef4444', 7, 0.85),
            makeTrace(down, 'Down', '#3b82f6', 7, 0.85),
        ];

        // Top gene labels
        const sigGenes = [...up, ...down].sort((a, b) => {
            const pa = isRawP ? a.pvalue : a.adj_pvalue;
            const pb = isRawP ? b.pvalue : b.adj_pvalue;
            return pa - pb;
        }).slice(0, 10);
        if (sigGenes.length > 0) {
            traces.push({
                x: sigGenes.map(r => r.log2fc),
                y: sigGenes.map(r => r.neg_log10_p),
                text: sigGenes.map(r => r.gene),
                mode: 'text',
                textposition: 'top center',
                textfont: { size: 9, color: '#e2e8f0' },
                showlegend: false,
                hoverinfo: 'skip',
            });
        }

        // Compute x range
        const allFc = results.map(r => r.log2fc);
        const maxFc = Math.max(Math.abs(Math.min(...allFc)), Math.abs(Math.max(...allFc)), fcThresh + 1);

        const layout = {
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            font: { family: 'Inter, sans-serif', color: '#94a3b8', size: 11 },
            xaxis: {
                title: { text: 'log2 Fold Change', font: { size: 12 } },
                gridcolor: 'rgba(255,255,255,0.04)',
                zeroline: false,
                range: [-maxFc * 1.1, maxFc * 1.1],
            },
            yaxis: {
                title: { text: `-log10(${pLabel})`, font: { size: 12 } },
                gridcolor: 'rgba(255,255,255,0.04)',
                zeroline: false,
            },
            legend: { x: 0.02, y: 0.98, bgcolor: 'rgba(0,0,0,0.3)', font: { size: 10 } },
            margin: { l: 60, r: 20, t: 10, b: 60 },
            shapes: [
                // Vertical FC thresholds
                { type: 'line', x0: fcThresh, x1: fcThresh, y0: 0, y1: 1, yref: 'paper', line: { color: 'rgba(255,255,255,0.15)', dash: 'dash', width: 1 } },
                { type: 'line', x0: -fcThresh, x1: -fcThresh, y0: 0, y1: 1, yref: 'paper', line: { color: 'rgba(255,255,255,0.15)', dash: 'dash', width: 1 } },
                // Horizontal p-value threshold
                { type: 'line', x0: 0, x1: 1, xref: 'paper', y0: negLog10PThresh, y1: negLog10PThresh, line: { color: 'rgba(255,255,255,0.15)', dash: 'dash', width: 1 } },
            ],
        };

        Plotly.newPlot('volcano-plot', traces, layout, {
            responsive: true, displayModeBar: false,
        });
    },

    populateDegTable(results) {
        const tbody = document.querySelector('#deg-table tbody');
        tbody.innerHTML = '';

        // Show top 50 significant genes
        const sig = results.filter(r => r.direction !== 'ns').slice(0, 50);
        if (sig.length === 0) {
            // Fallback: show top 50 by p-value
            sig.push(...results.slice(0, 50));
        }

        sig.forEach((gene, i) => {
            const tr = document.createElement('tr');
            tr.innerHTML = `
                <td style="font-weight:600;color:var(--text-accent)">${i + 1}</td>
                <td style="font-weight:600">${gene.gene}</td>
                <td>${gene.log2fc.toFixed(3)}</td>
                <td>${gene.pvalue.toExponential(2)}</td>
                <td>${gene.adj_pvalue.toExponential(2)}</td>
                <td class="dir-${gene.direction}">${gene.direction === 'up' ? '▲ Up' : gene.direction === 'down' ? '▼ Down' : '— NS'}</td>
            `;
            tbody.appendChild(tr);
        });
    },
};
