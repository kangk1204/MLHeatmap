/* Biomarker Panel Logic */
const Biomarker = {
    init() {
        document.getElementById('btn-run-biomarker').addEventListener('click', () => this.run());
        document.getElementById('btn-to-export').addEventListener('click', () => {
            App.goToPanel('heatmap');
            Heatmap.renderShapHeatmap();
        });

        // Model-specific parameter visibility
        const modelSelect = document.getElementById('model-select');
        if (modelSelect) {
            modelSelect.addEventListener('change', () => this._updateModelParams());
            this._updateModelParams();
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
        const es = API.biomarkerStream(App.state.sessionId, {
            nTopGenes: parseInt(document.getElementById('n-top-genes').value),
            nEstimators: parseInt(document.getElementById('n-estimators').value),
            cvFolds: parseInt(document.getElementById('cv-folds').value),
            model: modelSelect ? modelSelect.value : 'rf',
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

        // AUC vs #genes curve
        const curve = combo.auc_curve;
        const trace = {
            x: curve.map(c => c.n_genes),
            y: curve.map(c => c.auc),
            mode: 'lines+markers',
            line: { color: '#8b5cf6', width: 2.5 },
            marker: { size: 8, color: '#8b5cf6', line: { color: '#fff', width: 1 } },
            text: curve.map(c => `+${c.gene_added}`),
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

        // Gene list chips
        const listEl = document.getElementById('optimal-gene-list');
        listEl.innerHTML = '<h4 style="margin:0 0 8px;color:var(--text-primary)">Optimal Gene Set:</h4>';
        const chipContainer = document.createElement('div');
        chipContainer.style.cssText = 'display:flex;flex-wrap:wrap;gap:6px;';
        combo.best_genes.forEach((gene, i) => {
            const chip = document.createElement('span');
            chip.className = 'optimal-gene-chip';
            chip.textContent = `${i + 1}. ${gene}`;
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
};
