/* Main Application */
const App = {
    state: {
        sessionId: null,
        sampleNames: [],
        species: 'unknown',
        idType: 'unknown',
        groups: {},
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
        document.getElementById('btn-to-groups').addEventListener('click', () => {
            this.goToPanel('groups');
            Groups.populate(this.state.sampleNames);
        });
    },

    goToPanel(panelId) {
        document.querySelectorAll('.panel').forEach(p => p.classList.remove('active'));
        document.querySelectorAll('.nav-step').forEach(s => s.classList.remove('active'));

        document.getElementById(`panel-${panelId}`).classList.add('active');
        const navBtn = document.querySelector(`.nav-step[data-panel="${panelId}"]`);
        if (navBtn) navBtn.classList.add('active');

        this.state.currentPanel = panelId;
    },

    markStepCompleted(stepId) {
        this.completedSteps.add(stepId);
        const navBtn = document.querySelector(`.nav-step[data-panel="${stepId}"]`);
        if (navBtn) navBtn.classList.add('completed');
    },

    async mapGenes() {
        if (!this.state.sessionId) return this.showToast('Upload data first', 'error');

        const species = document.querySelector('input[name="species"]:checked').value;
        const idType = document.getElementById('id-type-select').value;

        this.showLoading('Mapping gene IDs...');
        try {
            const result = await API.mapGenes(this.state.sessionId, species, idType);

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
        // Clear state
        this.state.sessionId = null;
        this.state.sampleNames = [];
        this.state.species = 'unknown';
        this.state.idType = 'unknown';
        this.state.groups = {};
        this.completedSteps.clear();

        // Reset sidebar step indicators
        document.querySelectorAll('.nav-step').forEach(s => s.classList.remove('completed'));

        // Hide result sections
        document.querySelectorAll('#upload-result, #mapping-result, #norm-result, #biomarker-results').forEach(el => el.classList.add('hidden'));
        document.getElementById('upload-zone').style.display = '';

        // Clear plots
        ['heatmap-plot', 'dist-plot', 'shap-plot', 'auc-plot', 'optimal-auc-curve'].forEach(id => {
            const el = document.getElementById(id);
            if (el) Plotly.purge(el);
        });

        // Reset groups
        Groups.groupCount = 0;
        Groups.selectedSamples.clear();
        document.getElementById('groups-area').innerHTML = '';
        document.getElementById('sample-pool').innerHTML = '';

        // Go to upload
        this.goToPanel('upload');
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
