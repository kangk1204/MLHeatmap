/* Export Panel Logic */
const Export = {
    init() {
        document.querySelectorAll('.export-formats .btn').forEach(btn => {
            btn.addEventListener('click', () => this.download(btn.dataset.type));
        });
    },

    /** Show / hide export cards based on which analyses have been run. */
    refresh() {
        const hasDeg  = !!App.state.degResults;
        const hasMl   = !!App.state.biomarkerResults;

        // Volcano — requires DEG
        const volcanoCard = document.querySelector('.export-card[data-type="volcano_png"]');
        if (volcanoCard) volcanoCard.style.display = hasDeg ? '' : 'none';

        // SHAP — requires ML
        const shapCard = document.querySelector('.export-card[data-type="shap_png"]');
        if (shapCard) shapCard.style.display = hasMl ? '' : 'none';

        // ROC — requires ML
        const aucCard = document.querySelector('.export-card[data-type="auc_png"]');
        if (aucCard) aucCard.style.display = hasMl ? '' : 'none';
    },

    download(type) {
        if (!App.state.sessionId) {
            App.showToast('No data loaded', 'error');
            return;
        }

        // Guard: block exports for data that doesn't exist
        const needsDeg = type.startsWith('volcano_');
        const needsMl  = type.startsWith('shap_') || type.startsWith('auc_');
        if (needsDeg && !App.state.degResults) {
            App.showToast('Run DEG analysis first', 'error');
            return;
        }
        if (needsMl && !App.state.biomarkerResults) {
            App.showToast('Run biomarker analysis first', 'error');
            return;
        }

        // Client-side export for heatmap (preserves exact visual appearance)
        if (type === 'heatmap_png' || type === 'heatmap_svg') {
            this._exportHeatmapClientSide(type);
            return;
        }

        const dpi = document.getElementById('export-dpi').value;
        const url = API.exportUrl(App.state.sessionId, type, dpi);

        const a = document.createElement('a');
        a.href = url;
        a.download = '';
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);

        App.showToast(`Downloading ${type}...`, 'info');
    },

    _exportHeatmapClientSide(type) {
        const plotEl = document.getElementById('heatmap-plot');
        if (!plotEl || !plotEl.data || plotEl.data.length === 0) {
            App.showToast('Render heatmap first', 'error');
            return;
        }

        const dpiSelect = document.getElementById('export-dpi');
        const dpi = dpiSelect ? parseInt(dpiSelect.value) : 300;
        const scale = dpi / 96;  // 96 DPI is screen resolution
        const fmt = type === 'heatmap_svg' ? 'svg' : 'png';

        // Get current layout dimensions
        const currentLayout = plotEl.layout || {};
        const width = currentLayout.width || plotEl.offsetWidth || 1200;
        const height = currentLayout.height || 800;

        // Export with dark background matching the UI
        Plotly.downloadImage(plotEl, {
            format: fmt,
            width: width,
            height: height,
            scale: scale,
            filename: `heatmap_${new Date().toISOString().slice(0, 10)}`,
            setBackground: function() { return '#0a0a12'; },
        }).then(() => {
            App.showToast(`Heatmap ${fmt.toUpperCase()} exported`, 'success');
        }).catch(err => {
            App.showToast(`Export failed: ${err.message}`, 'error');
        });
    },
};
