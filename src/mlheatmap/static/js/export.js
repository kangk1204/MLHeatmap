/* Export Panel Logic */
const Export = {
    init() {
        document.querySelectorAll('.export-formats .btn').forEach(btn => {
            btn.addEventListener('click', () => this.download(btn.dataset.type));
        });
    },

    download(type) {
        if (!App.state.sessionId) {
            App.showToast('No data loaded', 'error');
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
