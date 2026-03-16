/* Export Panel Logic */
const Export = {
    init() {
        document.querySelectorAll('.export-formats .btn').forEach(btn => {
            btn.addEventListener('click', () => this.download(btn.dataset.type));
        });
    },

    /** Show / hide export cards based on which analyses have been run. */
    refresh() {
        const hasDeg = !!App.state.degResults;
        const hasMl = !!App.state.biomarkerResults;

        const volcanoCard = document.querySelector('.export-card[data-type="volcano_png"]');
        if (volcanoCard) volcanoCard.style.display = hasDeg ? '' : 'none';

        const shapCard = document.querySelector('.export-card[data-type="shap_png"]');
        if (shapCard) shapCard.style.display = hasMl ? '' : 'none';

        const aucCard = document.querySelector('.export-card[data-type="auc_png"]');
        if (aucCard) aucCard.style.display = hasMl ? '' : 'none';

        const optimalCard = document.querySelector('.export-card[data-type="optimal_combo_png"]');
        if (optimalCard) optimalCard.style.display = hasMl ? '' : 'none';
    },

    download(type) {
        if (!App.state.sessionId) {
            App.showToast('No data loaded', 'error');
            return;
        }

        const needsDeg = type.startsWith('volcano_');
        const needsMl = type.startsWith('shap_') || type.startsWith('auc_') || type.startsWith('optimal_combo_');
        if (needsDeg && !App.state.degResults) {
            App.showToast('Run DEG analysis first', 'error');
            return;
        }
        if (needsMl && !App.state.biomarkerResults) {
            App.showToast('Run biomarker analysis first', 'error');
            return;
        }

        if (type === 'heatmap_png' || type === 'heatmap_svg') {
            this._exportHeatmapClientSide(type);
            return;
        }

        if (type === 'shap_png' || type === 'shap_svg') {
            this._exportPlotlyElement('shap-plot', type, 'shap_plot');
            return;
        }

        if (type === 'auc_png' || type === 'auc_svg') {
            this._exportPlotlyElement('auc-plot', type, 'roc_curve');
            return;
        }

        if (type === 'optimal_combo_png' || type === 'optimal_combo_svg') {
            this._exportPlotlyElement('optimal-auc-curve', type, 'optimal_gene_combination_curve');
            return;
        }

        if (type === 'volcano_png' || type === 'volcano_svg') {
            this._exportPlotlyElement('volcano-plot', type, 'volcano_plot');
            return;
        }

        if (type === 'results_excel') {
            const dpi = document.getElementById('export-dpi').value;
            const url = API.exportUrl(App.state.sessionId, type, dpi);
            const anchor = document.createElement('a');
            anchor.href = url;
            anchor.download = '';
            document.body.appendChild(anchor);
            anchor.click();
            document.body.removeChild(anchor);
            App.showToast('Downloading results workbook...', 'info');
        }
    },

    _exportHeatmapClientSide(type) {
        const serverViewer = document.getElementById('server-heatmap-viewer');
        if (serverViewer) {
            this._downloadServerHeatmap(type);
            return;
        }

        const plotEl = document.getElementById('heatmap-plot');
        if (!plotEl || !plotEl.data || plotEl.data.length === 0) {
            App.showToast('Render heatmap first', 'error');
            return;
        }

        const dpiSelect = document.getElementById('export-dpi');
        const dpi = dpiSelect ? parseInt(dpiSelect.value, 10) : 300;
        const scale = dpi / 96;
        const fmt = type === 'heatmap_svg' ? 'svg' : 'png';
        const currentLayout = plotEl.layout || {};
        const width = currentLayout.width || plotEl.offsetWidth || 1200;
        const height = currentLayout.height || 800;

        Plotly.downloadImage(plotEl, {
            format: fmt,
            width,
            height,
            scale,
            filename: `heatmap_${new Date().toISOString().slice(0, 10)}`,
            setBackground: function() { return '#0a0a12'; },
        }).then(() => {
            App.showToast(`Heatmap ${fmt.toUpperCase()} exported`, 'success');
        }).catch(err => {
            App.showToast(`Export failed: ${err.message}`, 'error');
        });
    },

    async _downloadServerHeatmap(type) {
        try {
            const dpiSelect = document.getElementById('export-dpi');
            const clusterRowsEl = document.getElementById('cluster-rows');
            const clusterColsEl = document.getElementById('cluster-cols');
            const imageUrl = await API.getHeatmapImage(App.state.sessionId, {
                topN: parseInt(
                    document.getElementById('topn-input').value || document.getElementById('topn-slider').value,
                    10,
                ),
                distance: document.getElementById('distance-select').value,
                linkage: document.getElementById('linkage-select').value,
                colorScale: document.getElementById('colorscale-select').value,
                clusterRows: clusterRowsEl ? clusterRowsEl.checked : true,
                clusterCols: clusterColsEl ? clusterColsEl.checked : true,
                dpi: dpiSelect ? parseInt(dpiSelect.value, 10) : 300,
                fmt: type === 'heatmap_svg' ? 'svg' : 'png',
            });

            const anchor = document.createElement('a');
            anchor.href = imageUrl;
            anchor.download = type === 'heatmap_svg' ? 'heatmap.svg' : 'heatmap.png';
            document.body.appendChild(anchor);
            anchor.click();
            document.body.removeChild(anchor);
            window.setTimeout(() => URL.revokeObjectURL(imageUrl), 1000);
            App.showToast(`Heatmap ${type.endsWith('_svg') ? 'SVG' : 'PNG'} exported`, 'success');
        } catch (err) {
            App.showToast(`Export failed: ${err.message}`, 'error');
        }
    },

    _exportPlotlyElement(elementId, type, filenameBase) {
        const plotEl = document.getElementById(elementId);
        if (!plotEl || !plotEl.data || plotEl.data.length === 0) {
            App.showToast('Generate the plot first', 'error');
            return;
        }

        const dpiSelect = document.getElementById('export-dpi');
        const dpi = dpiSelect ? parseInt(dpiSelect.value, 10) : 300;
        const scale = dpi / 96;
        const fmt = type.endsWith('_svg') ? 'svg' : 'png';
        const layout = plotEl.layout || {};
        const width = layout.width || plotEl.offsetWidth || 1200;
        const height = layout.height || 800;

        Plotly.downloadImage(plotEl, {
            format: fmt,
            width,
            height,
            scale,
            filename: `${filenameBase}_${new Date().toISOString().slice(0, 10)}`,
            setBackground: function() { return '#0a0a12'; },
        }).then(() => {
            App.showToast(`${filenameBase} ${fmt.toUpperCase()} exported`, 'success');
        }).catch(err => {
            App.showToast(`Export failed: ${err.message}`, 'error');
        });
    },
};
