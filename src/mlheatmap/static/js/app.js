/* Main Application */
const App = {
    state: {
        sessionId: null,
        sampleNames: [],
        species: 'unknown',
        idType: 'unknown',
        groups: {},
        excludedSamples: [],
        currentPanel: 'upload',
    },

    completedSteps: new Set(),

    init() {
        Upload.init();
        Heatmap.init();
        Groups.init();
        Biomarker.init();
        Export.init();

        // Logo click → full reset
        document.querySelector('.logo').addEventListener('click', () => this.resetAll());

        // Sidebar navigation
        document.querySelectorAll('.nav-step').forEach(btn => {
            btn.addEventListener('click', () => {
                const panel = btn.dataset.panel;
                this.goToPanel(panel);
            });
        });

        // Gene mapping button
        document.getElementById('btn-map-genes').addEventListener('click', () => this.mapGenes());
        document.getElementById('btn-to-normalize').addEventListener('click', () => this.goToPanel('normalize'));

        // Normalize
        document.querySelectorAll('.method-card').forEach(card => {
            card.addEventListener('click', () => {
                document.querySelectorAll('.method-card').forEach(c => c.classList.remove('selected'));
                card.classList.add('selected');
            });
        });
        document.getElementById('btn-normalize').addEventListener('click', () => this.normalize());
        document.getElementById('btn-to-groups').addEventListener('click', () => this.goToPanel('groups'));
    },

    _activatePanel(panelId) {
        document.querySelectorAll('.panel').forEach(p => p.classList.remove('active'));
        document.querySelectorAll('.nav-step').forEach(s => s.classList.remove('active'));

        document.getElementById(`panel-${panelId}`).classList.add('active');
        const navBtn = document.querySelector(`.nav-step[data-panel="${panelId}"]`);
        if (navBtn) navBtn.classList.add('active');

        this.state.currentPanel = panelId;

        if (panelId === 'groups' && typeof Groups !== 'undefined') {
            Groups.populate(this.state.sampleNames);
        }

        // Refresh export card visibility when entering export panel
        if (panelId === 'export' && typeof Export !== 'undefined') {
            Export.refresh();
        }
    },

    async goToPanel(panelId) {
        if (panelId === this.state.currentPanel) return;

        const leavingGroups = this.state.currentPanel === 'groups' && panelId !== 'groups';
        if (leavingGroups && typeof Groups !== 'undefined' && typeof Groups.persistIfDirty === 'function') {
            const requiresValidGroups = ['biomarker', 'heatmap', 'export'].includes(panelId);
            const saved = await Groups.persistIfDirty({ requireValid: requiresValidGroups });
            if (!saved && requiresValidGroups) return;
        }

        this._activatePanel(panelId);
    },

    markStepCompleted(stepId) {
        this.completedSteps.add(stepId);
        const navBtn = document.querySelector(`.nav-step[data-panel="${stepId}"]`);
        if (navBtn) navBtn.classList.add('completed');
    },

    clearCompletedSteps(stepIds) {
        stepIds.forEach(stepId => {
            this.completedSteps.delete(stepId);
            const navBtn = document.querySelector(`.nav-step[data-panel="${stepId}"]`);
            if (navBtn) navBtn.classList.remove('completed');
        });
    },

    invalidateAnalysisState({ clearNormalization = false } = {}) {
        this.state.biomarkerResults = null;
        this.state.degResults = null;
        if (clearNormalization) this.state.totalGenes = null;
        if (typeof Biomarker !== 'undefined' && typeof Biomarker.cancelPending === 'function') {
            Biomarker.cancelPending();
        }
        if (typeof Heatmap !== 'undefined' && typeof Heatmap.cancelPending === 'function') {
            Heatmap.cancelPending();
        }

        const clearPlot = (id) => {
            const el = document.getElementById(id);
            if (!el) return;
            Plotly.purge(el);
            el.innerHTML = '';
        };

        ['heatmap-plot', 'shap-plot', 'auc-plot', 'optimal-auc-curve', 'volcano-plot'].forEach(clearPlot);
        if (clearNormalization) clearPlot('dist-plot');

        ['biomarker-results', 'deg-results'].forEach(id => {
            const el = document.getElementById(id);
            if (el) el.classList.add('hidden');
        });
        if (clearNormalization) {
            const normResult = document.getElementById('norm-result');
            if (normResult) normResult.classList.add('hidden');
        }

        const aucNote = document.getElementById('panel-auc-note');
        if (aucNote) aucNote.classList.add('hidden');
        const rocTag = document.getElementById('roc-eval-tag');
        if (rocTag) {
            rocTag.textContent = '';
            rocTag.className = 'eval-tag';
        }
        const renderBadge = document.getElementById('render-mode-badge');
        if (renderBadge) renderBadge.style.display = 'none';

        if (typeof Heatmap !== 'undefined') {
            if (Heatmap._serverImageUrl) {
                URL.revokeObjectURL(Heatmap._serverImageUrl);
                Heatmap._serverImageUrl = null;
            }
            const maxGenes = clearNormalization ? 60000 : (this.state.totalGenes || Heatmap._maxGenes || 60000);
            Heatmap._isShapMode = false;
            Heatmap._isDegMode = false;
            Heatmap._maxGenes = maxGenes;
            Heatmap._activeMax = maxGenes;

            const topnSlider = document.getElementById('topn-slider');
            const topnInput = document.getElementById('topn-input');
            const currentValue = clearNormalization
                ? 500
                : Math.min(
                    parseInt(topnInput?.value || topnSlider?.value || '500', 10) || 500,
                    maxGenes,
                );
            if (topnSlider) {
                topnSlider.max = maxGenes;
                topnSlider.step = typeof Heatmap._getTopNStep === 'function'
                    ? Heatmap._getTopNStep(currentValue)
                    : (maxGenes > 10000 ? 1000 : maxGenes > 5000 ? 500 : 50);
                topnSlider.value = currentValue;
            }
            if (topnInput) {
                topnInput.max = maxGenes;
                topnInput.value = currentValue;
            }
        }

        this.clearCompletedSteps(clearNormalization ? ['normalize', 'heatmap', 'biomarker'] : ['heatmap', 'biomarker']);

        if (typeof Export !== 'undefined') {
            Export.refresh();
        }
    },

    async mapGenes() {
        if (!this.state.sessionId) return this.showToast('Upload data first', 'error');

        const species = document.querySelector('input[name="species"]:checked').value;
        const idType = document.getElementById('id-type-select').value;

        this.showLoading('Mapping gene IDs...');
        try {
            const result = await API.mapGenes(this.state.sessionId, species, idType);
            this.invalidateAnalysisState({ clearNormalization: true });

            document.getElementById('mapped-count').textContent = result.mapped_count.toLocaleString();
            document.getElementById('unmapped-count').textContent = result.unmapped_count.toLocaleString();
            document.getElementById('total-count').textContent = result.total.toLocaleString();

            document.getElementById('mapping-result').classList.remove('hidden');
            this.markStepCompleted('mapping');

            if (result.unmapped_count > result.total * 0.1) {
                this.showToast(`Warning: ${result.unmapped_count} genes unmapped (${(result.unmapped_count / result.total * 100).toFixed(1)}%)`, 'error');
            } else {
                this.showToast('Gene mapping complete', 'success');
            }
        } catch (err) {
            this.showToast(err.message, 'error');
        } finally {
            this.hideLoading();
        }
    },

    async normalize() {
        if (!this.state.sessionId) return this.showToast('Upload data first', 'error');

        const selected = document.querySelector('.method-card.selected');
        const method = selected ? selected.dataset.method : 'deseq2';

        this.showLoading('Normalizing...');
        try {
            const result = await API.normalize(this.state.sessionId, method);
            this.invalidateAnalysisState();

            document.getElementById('norm-min').textContent = result.stats.min.toFixed(2);
            document.getElementById('norm-median').textContent = result.stats.median.toFixed(2);
            document.getElementById('norm-max').textContent = result.stats.max.toFixed(2);

            // Distribution plot
            if (result.distribution_sample) {
                const trace = {
                    x: result.distribution_sample,
                    type: 'histogram',
                    nbinsx: 50,
                    marker: {
                        color: 'rgba(59, 130, 246, 0.6)',
                        line: { color: 'rgba(59, 130, 246, 0.8)', width: 1 },
                    },
                };
                const layout = {
                    paper_bgcolor: 'rgba(0,0,0,0)',
                    plot_bgcolor: 'rgba(0,0,0,0)',
                    font: { family: 'Inter, sans-serif', color: '#94a3b8', size: 11 },
                    xaxis: {
                        title: { text: 'Normalized value', font: { size: 12 } },
                        gridcolor: 'rgba(255,255,255,0.04)',
                    },
                    yaxis: {
                        title: { text: 'Frequency', font: { size: 12 } },
                        gridcolor: 'rgba(255,255,255,0.04)',
                    },
                    margin: { l: 60, r: 20, t: 10, b: 50 },
                    bargap: 0.02,
                };
                Plotly.newPlot('dist-plot', [trace], layout, { responsive: true, displayModeBar: false });
            }

            // Store total gene count for slider capping
            if (result.shape && result.shape.length >= 1) {
                this.state.totalGenes = result.shape[0];
                Heatmap.updateTopNMax(result.shape[0]);
            }

            document.getElementById('norm-result').classList.remove('hidden');
            this.markStepCompleted('normalize');
            this.showToast(`Normalization complete (${method})`, 'success');
        } catch (err) {
            this.showToast(err.message, 'error');
        } finally {
            this.hideLoading();
        }
    },

    resetAll() {
        if (typeof Biomarker !== 'undefined' && typeof Biomarker.cancelPending === 'function') {
            Biomarker.cancelPending();
        }
        if (typeof Heatmap !== 'undefined' && typeof Heatmap.cancelPending === 'function') {
            Heatmap.cancelPending();
        }

        // Clear state
        this.state.sessionId = null;
        this.state.sampleNames = [];
        this.state.species = 'unknown';
        this.state.idType = 'unknown';
        this.state.groups = {};
        this.state.excludedSamples = [];
        this.state.biomarkerResults = null;
        this.state.degResults = null;
        this.state.totalGenes = null;
        this.completedSteps.clear();

        // Reset sidebar step indicators
        document.querySelectorAll('.nav-step').forEach(s => s.classList.remove('completed'));

        // Hide result sections
        document.querySelectorAll('#upload-result, #mapping-result, #norm-result, #biomarker-results').forEach(el => el.classList.add('hidden'));
        document.getElementById('upload-zone').style.display = '';

        // Clear plots
        ['heatmap-plot', 'dist-plot', 'shap-plot', 'auc-plot', 'optimal-auc-curve', 'volcano-plot'].forEach(id => {
            const el = document.getElementById(id);
            if (el) Plotly.purge(el);
        });

        // Reset heatmap mode flags and slider max
        Heatmap._isShapMode = false;
        Heatmap._isDegMode = false;
        Heatmap._maxGenes = 60000;
        Heatmap._activeMax = 60000;
        const topnSlider = document.getElementById('topn-slider');
        const topnInput = document.getElementById('topn-input');
        if (topnSlider) { topnSlider.max = 60000; topnSlider.value = 500; }
        if (topnInput) { topnInput.max = 60000; topnInput.value = 500; }

        // Hide DEG results
        const degResults = document.getElementById('deg-results');
        if (degResults) degResults.classList.add('hidden');

        // Reset groups
        Groups.groupCount = 0;
        Groups.selectedSamples.clear();
        if (typeof Groups.markClean === 'function') Groups.markClean();
        document.getElementById('groups-area').innerHTML = '';
        document.getElementById('sample-pool').innerHTML = '';

        // Go to upload
        this._activatePanel('upload');
    },

    showLoading(text = 'Processing...') {
        document.getElementById('loading-text').textContent = text;
        document.getElementById('loading-overlay').classList.remove('hidden');
    },

    hideLoading() {
        document.getElementById('loading-overlay').classList.add('hidden');
    },

    showToast(message, type = 'info') {
        const container = document.getElementById('toast-container');
        const toast = document.createElement('div');
        toast.className = `toast ${type}`;
        toast.textContent = message;
        container.appendChild(toast);

        setTimeout(() => {
            toast.style.opacity = '0';
            toast.style.transform = 'translateY(12px)';
            toast.style.transition = 'all 0.3s ease';
            setTimeout(() => toast.remove(), 300);
        }, 4000);
    },
};

// Initialize when DOM is ready
document.addEventListener('DOMContentLoaded', () => App.init());
