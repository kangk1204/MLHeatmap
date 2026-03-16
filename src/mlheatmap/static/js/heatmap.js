/* Heatmap Panel Logic - SOTA-quality with dendrograms & group color bar */
const Heatmap = {
    _isShapMode: false,
    _isDegMode: false,
    _requestToken: 0,
    _pendingRequest: false,
    _abortController: null,
    _serverImageUrl: null,  // Track blob URL for cleanup
    _maxGenes: 60000,  // Updated after normalization
    _activeMax: 60000, // Current mode's effective max (changes per mode)
    SERVER_RENDER_THRESHOLD: 5000,
    GROUP_COLORS: ['#3b82f6', '#ef4444', '#10b981', '#f59e0b', '#8b5cf6',
                   '#ec4899', '#06b6d4', '#84cc16', '#f97316', '#6366f1'],

    _escapeHtml(str) {
        const div = document.createElement('div');
        div.textContent = str;
        return div.innerHTML;
    },

    _buildTopLeftTitleAnnotation(text) {
        return {
            text: this._escapeHtml(text),
            xref: 'paper',
            yref: 'paper',
            x: 0,
            y: 1.04,
            showarrow: false,
            font: { size: 12, color: '#cbd5e1' },
            xanchor: 'left',
            yanchor: 'bottom',
            align: 'left',
        };
    },

    cancelPending() {
        const hadPending = this._pendingRequest;
        this._requestToken += 1;
        this._pendingRequest = false;
        if (this._abortController) {
            this._abortController.abort();
            this._abortController = null;
        }
        if (hadPending && App.state.sessionId) {
            API.cancelSession(App.state.sessionId).catch(() => {});
        }
        if (hadPending) App.hideLoading();
    },

    _newAbortController() {
        if (this._abortController) this._abortController.abort();
        this._abortController = new AbortController();
        return this._abortController;
    },

    /** Consolidate dendrogram branches into a single trace using null separators. */
    _buildDendroTrace(dendroData, {type, xaxis, yaxis}) {
        const { icoord, dcoord } = dendroData;
        const allX = [];
        const allY = [];
        for (let i = 0; i < icoord.length; i++) {
            if (type === 'col') {
                // Column dendrogram: x = mapped icoord, y = raw dcoord
                for (let j = 0; j < icoord[i].length; j++) allX.push((icoord[i][j] - 5) / 10);
                allY.push(...dcoord[i]);
            } else {
                // Row dendrogram: x = -dcoord, y = mapped icoord
                for (let j = 0; j < dcoord[i].length; j++) allX.push(-dcoord[i][j]);
                for (let j = 0; j < icoord[i].length; j++) allY.push((icoord[i][j] - 5) / 10);
            }
            allX.push(null);
            allY.push(null);
        }
        return {
            x: allX, y: allY,
            mode: 'lines', line: { color: 'rgba(148, 163, 184, 0.6)', width: 1.2 },
            xaxis, yaxis, showlegend: false, hoverinfo: 'skip',
        };
    },

    _getTopNStep(val) {
        if (this._isShapMode) return 1;
        if (this._isDegMode) return 5;
        if (val > 10000) return 1000;
        if (val > 5000) return 500;
        if (val > 1000) return 100;
        if (val > 200) return 25;
        return 1;
    },

    init() {
        const slider = document.getElementById('topn-slider');
        const topnInput = document.getElementById('topn-input');

        const syncTopN = (val) => {
            const effectiveMax = this._activeMax || this._maxGenes || 60000;
            const clamped = Math.max(10, Math.min(effectiveMax, val));
            slider.step = this._getTopNStep(clamped);
            slider.value = clamped;
            topnInput.value = clamped;

            // Render mode badge
            const badge = document.getElementById('render-mode-badge');
            if (badge) {
                if (!this._isShapMode && !this._isDegMode && clamped > this.SERVER_RENDER_THRESHOLD) {
                    badge.textContent = 'Server Render';
                    badge.style.display = 'inline-block';
                } else {
                    badge.style.display = 'none';
                }
            }
        };

        slider.addEventListener('input', () => {
            syncTopN(parseInt(slider.value, 10) || 10);
        });

        topnInput.addEventListener('change', () => {
            syncTopN(parseInt(topnInput.value, 10) || 10);
        });

        document.getElementById('btn-render-heatmap').addEventListener('click', () => {
            if (this._isDegMode && App.state.degResults) {
                this.renderDegHeatmap();
            } else if (this._isShapMode && App.state.biomarkerResults) {
                this.renderShapHeatmap();
            } else {
                this.render();
            }
        });
        document.getElementById('btn-to-biomarker').addEventListener('click', () => App.goToPanel('export'));

        // Dynamic heatmap size controls
        const widthSlider = document.getElementById('heatmap-width-slider');
        const heightSlider = document.getElementById('heatmap-height-slider');
        const widthInput = document.getElementById('heatmap-width-input');
        const heightInput = document.getElementById('heatmap-height-input');
        const widthLabel = document.getElementById('heatmap-width-value');
        const heightLabel = document.getElementById('heatmap-height-value');

        const syncWidth = (rawValue, { resize = true } = {}) => {
            const parsed = Number.isFinite(rawValue) ? rawValue : parseInt(rawValue, 10);
            const value = Math.max(0, Math.min(2400, Number.isFinite(parsed) ? parsed : 0));
            if (widthSlider) widthSlider.value = value;
            if (widthInput) widthInput.value = value;
            if (widthLabel) widthLabel.textContent = value === 0 ? 'Auto' : `${value}px`;
            if (resize) this._resizePlot();
        };

        const syncHeight = (rawValue, { resize = true } = {}) => {
            const parsed = Number.isFinite(rawValue) ? rawValue : parseInt(rawValue, 10);
            const value = Math.max(400, Math.min(3000, Number.isFinite(parsed) ? parsed : 800));
            if (heightSlider) heightSlider.value = value;
            if (heightInput) heightInput.value = value;
            if (heightLabel) heightLabel.textContent = `${value}px`;
            if (resize) this._resizePlot();
        };

        if (widthSlider) {
            widthSlider.addEventListener('input', () => {
                syncWidth(parseInt(widthSlider.value, 10));
            });
        }
        if (heightSlider) {
            heightSlider.addEventListener('input', () => {
                syncHeight(parseInt(heightSlider.value, 10));
            });
        }
        if (widthInput) {
            widthInput.addEventListener('input', () => {
                syncWidth(widthInput.value);
            });
        }
        if (heightInput) {
            heightInput.addEventListener('input', () => {
                syncHeight(heightInput.value);
            });
        }

        syncWidth(widthSlider ? parseInt(widthSlider.value, 10) : 0, { resize: false });
        syncHeight(heightSlider ? parseInt(heightSlider.value, 10) : 800, { resize: false });
    },

    updateTopNMax(totalGenes) {
        const safeMax = Math.max(totalGenes, 10);
        this._maxGenes = safeMax;
        this._activeMax = safeMax;
        const slider = document.getElementById('topn-slider');
        const topnInput = document.getElementById('topn-input');
        if (slider) {
            slider.max = safeMax;
            slider.step = this._getTopNStep(parseInt(slider.value, 10) || 10);
        }
        if (topnInput) topnInput.max = safeMax;
        // Clamp current values if needed
        if (slider && parseInt(slider.value, 10) > safeMax) {
            slider.value = safeMax;
        }
        if (topnInput && parseInt(topnInput.value, 10) > safeMax) {
            topnInput.value = safeMax;
        }
    },

    _resizePlot() {
        // Skip if showing server-rendered image
        if (document.getElementById('server-heatmap-viewer')) return;
        const plotEl = document.getElementById('heatmap-plot');
        if (!plotEl || !plotEl.data || plotEl.data.length === 0) return;

        const widthSlider = document.getElementById('heatmap-width-slider');
        const heightSlider = document.getElementById('heatmap-height-slider');

        const wVal = widthSlider ? parseInt(widthSlider.value) : 0;
        const hVal = heightSlider ? parseInt(heightSlider.value) : 800;

        // For auto width, use the container's actual width
        const containerWidth = plotEl.parentElement ? plotEl.parentElement.clientWidth : 0;
        const newWidth = wVal > 0 ? wVal : (containerWidth > 0 ? containerWidth : undefined);

        Plotly.relayout(plotEl, {
            width: newWidth,
            height: hVal,
        });
    },

    async render() {
        if (!App.state.sessionId) return App.showToast('No data loaded', 'error');
        this._isShapMode = false;
        this._isDegMode = false;

        // Restore slider to full gene range
        if (this._maxGenes) {
            this._activeMax = this._maxGenes;
            const slider = document.getElementById('topn-slider');
            const topnInput = document.getElementById('topn-input');
            if (slider) {
                slider.max = this._maxGenes;
                slider.step = this._getTopNStep(parseInt(topnInput?.value || slider.value || '10', 10));
            }
            if (topnInput) topnInput.max = this._maxGenes;
        }

        const topN = parseInt(document.getElementById('topn-input').value || document.getElementById('topn-slider').value, 10);

        // Auto-switch to server-side rendering for large gene sets
        if (topN > this.SERVER_RENDER_THRESHOLD) {
            return this.renderServerSide(topN, ++this._requestToken);
        }

        const requestToken = ++this._requestToken;
        this._pendingRequest = true;
        const controller = this._newAbortController();
        App.showLoading('Computing heatmap...');
        try {
            const clusterRowsEl = document.getElementById('cluster-rows');
            const clusterColsEl = document.getElementById('cluster-cols');
            const data = await API.getHeatmap(App.state.sessionId, {
                topN: topN,
                distance: document.getElementById('distance-select').value,
                linkage: document.getElementById('linkage-select').value,
                colorScale: document.getElementById('colorscale-select').value,
                clusterRows: clusterRowsEl ? clusterRowsEl.checked : true,
                clusterCols: clusterColsEl ? clusterColsEl.checked : true,
                signal: controller.signal,
            });
            if (requestToken !== this._requestToken) return;

            // Restore header for regular heatmap
            const header = document.querySelector('#panel-heatmap .panel-header h1');
            const subtitle = document.querySelector('#panel-heatmap .panel-header .panel-subtitle');
            if (header) header.textContent = 'Heatmap';
            if (subtitle) subtitle.textContent = 'Interactive clustered heatmap of expression data';

            this.plotHeatmap(data);
            App.markStepCompleted('heatmap');
        } catch (err) {
            if (requestToken !== this._requestToken) return;
            if (err.name === 'AbortError') return;
            App.showToast(err.message, 'error');
        } finally {
            if (this._abortController === controller) this._abortController = null;
            if (requestToken !== this._requestToken) return;
            this._pendingRequest = false;
            App.hideLoading();
        }
    },

    async renderServerSide(topN, requestToken = ++this._requestToken) {
        this._pendingRequest = true;
        const controller = this._newAbortController();
        App.showLoading(`Rendering ${topN.toLocaleString()} genes on server...`);
        try {
            // Clean up previous blob URL
            if (this._serverImageUrl) {
                URL.revokeObjectURL(this._serverImageUrl);
                this._serverImageUrl = null;
            }

            const clusterRowsEl = document.getElementById('cluster-rows');
            const clusterColsEl = document.getElementById('cluster-cols');
            const imageUrl = await API.getHeatmapImage(App.state.sessionId, {
                topN: topN,
                distance: document.getElementById('distance-select').value,
                linkage: document.getElementById('linkage-select').value,
                colorScale: document.getElementById('colorscale-select').value,
                clusterRows: clusterRowsEl ? clusterRowsEl.checked : true,
                clusterCols: clusterColsEl ? clusterColsEl.checked : true,
                dpi: 150,
                signal: controller.signal,
            });
            if (requestToken !== this._requestToken) {
                URL.revokeObjectURL(imageUrl);
                return;
            }

            this._serverImageUrl = imageUrl;

            // Update header
            const header = document.querySelector('#panel-heatmap .panel-header h1');
            const subtitle = document.querySelector('#panel-heatmap .panel-header .panel-subtitle');
            if (header) header.textContent = 'Heatmap';
            if (subtitle) subtitle.textContent = `Top ${topN.toLocaleString()} most variable genes (server-rendered)`;

            this.displayServerImage(imageUrl, topN);
            App.markStepCompleted('heatmap');
        } catch (err) {
            if (requestToken !== this._requestToken) return;
            if (err.name === 'AbortError') return;
            App.showToast('Server render failed: ' + err.message, 'error');
        } finally {
            if (this._abortController === controller) this._abortController = null;
            if (requestToken !== this._requestToken) return;
            this._pendingRequest = false;
            App.hideLoading();
        }
    },

    displayServerImage(imageUrl, topN) {
        const container = document.getElementById('heatmap-plot');
        // Clear any existing Plotly plot
        if (container.data) Plotly.purge(container);

        container.innerHTML = `
            <div class="server-heatmap-viewer" id="server-heatmap-viewer">
                <div class="server-heatmap-toolbar">
                    <button class="btn btn-xs btn-outline" id="btn-zoom-in" title="Zoom In">
                        <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="11" cy="11" r="8"/><line x1="21" y1="21" x2="16.65" y2="16.65"/><line x1="11" y1="8" x2="11" y2="14"/><line x1="8" y1="11" x2="14" y2="11"/></svg>
                    </button>
                    <button class="btn btn-xs btn-outline" id="btn-zoom-out" title="Zoom Out">
                        <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="11" cy="11" r="8"/><line x1="21" y1="21" x2="16.65" y2="16.65"/><line x1="8" y1="11" x2="14" y2="11"/></svg>
                    </button>
                    <button class="btn btn-xs btn-outline" id="btn-zoom-reset" title="Reset Zoom">Reset</button>
                    <button class="btn btn-xs btn-outline" id="btn-zoom-fit" title="Fit to View">Fit</button>
                    <button class="btn btn-xs btn-outline" id="btn-download-server-img" title="Download Image">
                        <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4"/><polyline points="7 10 12 15 17 10"/><line x1="12" y1="15" x2="12" y2="3"/></svg>
                        Download
                    </button>
                    <span class="server-heatmap-info">${topN.toLocaleString()} genes · server-rendered</span>
                </div>
                <div class="server-heatmap-scroll" id="server-heatmap-scroll">
                    <img src="${imageUrl}" id="server-heatmap-img" alt="Clustered Heatmap" />
                </div>
            </div>
        `;

        let scale = 1;
        const img = document.getElementById('server-heatmap-img');
        const scroll = document.getElementById('server-heatmap-scroll');

        const applyScale = () => {
            img.style.transform = `scale(${scale})`;
            img.style.transformOrigin = 'top left';
        };

        document.getElementById('btn-zoom-in').onclick = () => {
            scale = Math.min(scale * 1.3, 5);
            applyScale();
        };
        document.getElementById('btn-zoom-out').onclick = () => {
            scale = Math.max(scale / 1.3, 0.1);
            applyScale();
        };
        document.getElementById('btn-zoom-reset').onclick = () => {
            scale = 1;
            applyScale();
        };
        document.getElementById('btn-zoom-fit').onclick = () => {
            // Fit image to scroll container width
            if (img.naturalWidth > 0) {
                scale = scroll.clientWidth / img.naturalWidth;
                applyScale();
            }
        };
        document.getElementById('btn-download-server-img').onclick = () => {
            const a = document.createElement('a');
            a.href = imageUrl;
            a.download = `heatmap_${topN}genes.png`;
            a.click();
        };

        // Mouse wheel zoom
        scroll.addEventListener('wheel', (e) => {
            e.preventDefault();
            const delta = e.deltaY > 0 ? 0.9 : 1.1;
            scale = Math.min(Math.max(scale * delta, 0.1), 5);
            applyScale();
        }, { passive: false });

        // Auto-fit on load (use requestAnimationFrame to ensure layout is ready)
        img.onload = () => {
            requestAnimationFrame(() => {
                const sw = scroll.clientWidth || container.clientWidth || 800;
                if (img.naturalWidth > 0 && sw > 0) {
                    scale = Math.min(sw / img.naturalWidth, 1);
                    applyScale();
                }
            });
        };
    },

    plotHeatmap(data) {
        const { z, x, y, groups, color_scale, row_dendrogram, col_dendrogram } = data;

        const nGenes = y.length;
        const nSamples = x.length;

        // Group mapping
        const groupNames = Object.keys(groups || {});
        const hasGroups = groupNames.length >= 2;
        const sampleGroupMap = {};
        if (hasGroups) {
            groupNames.forEach((name, i) => {
                (groups[name] || []).forEach(sample => {
                    sampleGroupMap[sample] = name;
                });
            });
        }

        const hasRowDendro = row_dendrogram && row_dendrogram.icoord && row_dendrogram.icoord.length > 0;
        const hasColDendro = col_dendrogram && col_dendrogram.icoord && col_dendrogram.icoord.length > 0;

        const traces = [];

        // Layout domain proportions
        const rowDendroWidth = hasRowDendro ? 0.08 : 0;
        const geneNameGap = hasRowDendro ? 0.10 : 0;  // space for gene names between dendro and heatmap
        const colDendroHeight = hasColDendro ? 0.1 : 0;
        const groupBarHeight = hasGroups ? 0.03 : 0;
        const legendSpace = hasGroups ? 0.04 : 0;  // space above for group legend text
        const colorbarSpace = 0.08;

        // Domains — dendro → gene names gap → heatmap
        const xRowDendro = [0, rowDendroWidth > 0 ? rowDendroWidth - 0.01 : 0];
        const heatmapLeft = rowDendroWidth + geneNameGap;
        const xHeatmap = [heatmapLeft, 1 - colorbarSpace];
        const yHeatmap = [0, 1 - colDendroHeight - groupBarHeight - legendSpace];
        const yGroupBar = [1 - colDendroHeight - groupBarHeight - legendSpace + 0.005, 1 - colDendroHeight - legendSpace - 0.005];
        const yColDendro = [1 - colDendroHeight - legendSpace, 1 - legendSpace];
        // Legend annotation will go at y ~ 1 - legendSpace/2 (in the reserved top band)

        // Numeric x/y for heatmap
        const xNums = Array.from({length: nSamples}, (_, i) => i);
        const yNums = Array.from({length: nGenes}, (_, i) => i);

        // =============================================
        // 1. MAIN HEATMAP
        // =============================================
        const hoverTexts = z.map((row, gi) =>
            row.map((val, si) => {
                const group = sampleGroupMap[x[si]] || '';
                return `<b>${y[gi]}</b><br>${x[si]}${group ? ` (${group})` : ''}<br>Z-score: ${val.toFixed(2)}`;
            })
        );

        traces.push({
            z: z,
            x: xNums,
            y: yNums,
            type: 'heatmap',
            colorscale: this._getColorScale(color_scale),
            zmid: 0, zmin: -3, zmax: 3,
            hovertext: hoverTexts,
            hovertemplate: '%{hovertext}<extra></extra>',
            xaxis: 'x', yaxis: 'y',
            showscale: true,
            colorbar: {
                title: { text: 'Z-score', font: { color: '#e2e8f0', size: 11 } },
                tickfont: { color: '#94a3b8', size: 10 },
                thickness: 14, len: 0.4, y: 0.3, x: 1.01,
                outlinewidth: 0,
                tickvals: [-3, -2, -1, 0, 1, 2, 3],
            },
        });

        // =============================================
        // 2. COLUMN DENDROGRAM (top) — single consolidated trace
        // =============================================
        if (hasColDendro) {
            traces.push(this._buildDendroTrace(col_dendrogram, {type: 'col', xaxis: 'x', yaxis: 'y2'}));
        }

        // =============================================
        // 3. ROW DENDROGRAM (left of gene names) — single consolidated trace
        // =============================================
        if (hasRowDendro) {
            traces.push(this._buildDendroTrace(row_dendrogram, {type: 'row', xaxis: 'x3', yaxis: 'y6'}));
        }

        // =============================================
        // 4. GROUP COLOR BAR
        // =============================================
        if (hasGroups) {
            const nGroups = groupNames.length;
            const groupZ = [xNums.map(si => {
                const gi = groupNames.indexOf(sampleGroupMap[x[si]]);
                return gi >= 0 ? gi + 0.5 : -0.5;
            })];

            const groupColorscale = [];
            for (let i = 0; i < nGroups; i++) {
                groupColorscale.push([i / nGroups, this.GROUP_COLORS[i % this.GROUP_COLORS.length]]);
                groupColorscale.push([(i + 1) / nGroups, this.GROUP_COLORS[i % this.GROUP_COLORS.length]]);
            }

            traces.push({
                z: groupZ,
                x: xNums,
                y: [0],
                type: 'heatmap',
                colorscale: groupColorscale,
                zmin: 0, zmax: nGroups,
                showscale: false,
                xaxis: 'x', yaxis: 'y4',
                hovertemplate: '<b>%{customdata}</b><extra></extra>',
                customdata: [x.map(s => sampleGroupMap[s] || 'Unassigned')],
            });
        }

        // =============================================
        // LAYOUT
        // =============================================
        // Use slider value if available, otherwise auto-calculate
        const heightSlider = document.getElementById('heatmap-height-slider');
        const widthSlider = document.getElementById('heatmap-width-slider');
        const sliderHeight = heightSlider ? parseInt(heightSlider.value) : 0;
        const sliderWidth = widthSlider ? parseInt(widthSlider.value) : 0;
        const totalHeight = sliderHeight > 0 ? sliderHeight : Math.max(550, Math.min(2400, nGenes * 6 + 280));
        const sampleFontSize = Math.max(6, Math.min(11, 500 / nSamples));

        // Compute how many pixels per gene row we have
        const plotAreaHeight = totalHeight * (yHeatmap[1] - yHeatmap[0]) * 0.85;
        const pxPerGene = plotAreaHeight / nGenes;

        // Determine gene label visibility: show every Nth label so they don't overlap
        // Each label needs ~10px minimum vertical space
        const minPxPerLabel = 10;
        const labelStep = pxPerGene >= minPxPerLabel ? 1 : Math.ceil(minPxPerLabel / pxPerGene);
        const geneFontSize = Math.max(5, Math.min(11, pxPerGene * 0.85));

        // Filter tick labels: show every Nth gene
        const geneTickVals = [];
        const geneTickText = [];
        for (let i = 0; i < nGenes; i++) {
            if (i % labelStep === 0) {
                geneTickVals.push(yNums[i]);
                geneTickText.push(y[i]);
            }
        }

        // Compute left margin based on longest visible gene name
        const maxGeneLen = geneTickText.reduce((max, g) => Math.max(max, g.length), 0);
        const leftMargin = Math.max(80, Math.min(160, maxGeneLen * geneFontSize * 0.65 + 10));

        const layout = {
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            font: { family: 'Inter, sans-serif', color: '#94a3b8', size: 10 },
            height: totalHeight, width: sliderWidth > 0 ? sliderWidth : null,
            margin: { l: hasRowDendro ? 20 : leftMargin, r: 80, t: hasGroups ? 40 : 25, b: 10 },

            // Main heatmap X (shared by col dendro & group bar)
            xaxis: {
                domain: xHeatmap,
                tickvals: xNums,
                ticktext: x,
                tickangle: -45,
                tickfont: { size: sampleFontSize, color: '#94a3b8' },
                side: 'bottom',
                showgrid: false, zeroline: false,
                range: [-0.5, nSamples - 0.5],
            },

            // Main heatmap Y — gene names positioned at heatmap left edge
            yaxis: {
                domain: yHeatmap,
                anchor: hasRowDendro ? 'free' : undefined,
                position: hasRowDendro ? heatmapLeft : undefined,
                tickvals: geneTickVals,
                ticktext: geneTickText,
                tickfont: { size: geneFontSize, color: '#94a3b8' },
                autorange: 'reversed',
                showgrid: false, zeroline: false,
                range: [-0.5, nGenes - 0.5],
            },

            // Col dendrogram Y (top, own scale)
            yaxis2: {
                domain: yColDendro,
                showticklabels: false, showgrid: false,
                zeroline: false, showline: false,
            },

            // Row dendrogram X (left, own scale)
            xaxis3: {
                domain: xRowDendro,
                showticklabels: false, showgrid: false,
                zeroline: false, showline: false,
                autorange: true,
            },

            // Group color bar Y
            yaxis4: {
                domain: yGroupBar,
                showticklabels: false, showgrid: false,
                zeroline: false, showline: false,
                fixedrange: true,
            },

            // Row dendrogram Y (separate axis, mirrors heatmap y range)
            yaxis6: {
                domain: yHeatmap,
                anchor: 'x3',
                range: [-0.5, nGenes - 0.5],
                autorange: 'reversed',
                showticklabels: false, showgrid: false,
                zeroline: false, showline: false,
            },

            annotations: [],
        };

        // Group legend — placed in reserved top band above dendro/group bar
        if (hasGroups) {
            const parts = groupNames.map((name, i) => {
                const c = this.GROUP_COLORS[i % this.GROUP_COLORS.length];
                return `<span style="color:${c}">\u25a0</span> ${this._escapeHtml(name)}`;
            });
            layout.annotations.push({
                text: parts.join('&nbsp;&nbsp;&nbsp;'),
                xref: 'paper', yref: 'paper',
                x: (xHeatmap[0] + xHeatmap[1]) / 2,
                y: 1 - legendSpace / 2,
                showarrow: false,
                font: { size: 11, color: '#cbd5e1' },
                xanchor: 'center', yanchor: 'middle',
            });
        }

        // Gene count
        if (data.n_total_genes && data.n_shown_genes) {
            layout.annotations.push({
                text: `Showing ${data.n_shown_genes} of ${data.n_total_genes} genes (most variable)`,
                xref: 'paper', yref: 'paper',
                x: 0, y: -0.01,
                showarrow: false,
                font: { size: 10, color: '#64748b' },
                xanchor: 'left', yanchor: 'top',
            });
        }

        Plotly.newPlot('heatmap-plot', traces, layout, {
            responsive: true,
            displayModeBar: true,
            modeBarButtonsToRemove: ['lasso2d', 'select2d'],
            displaylogo: false,
            toImageButtonOptions: { format: 'png', height: 1600, width: 2000, scale: 2 },
        }).then(() => this._resizePlot());
    },

    async renderShapHeatmap() {
        if (!App.state.sessionId) return App.showToast('No data loaded', 'error');

        const bioResults = App.state.biomarkerResults;
        if (!bioResults) return App.showToast('Run biomarker analysis first', 'error');
        this._isShapMode = true;
        this._isDegMode = false;
        const requestToken = ++this._requestToken;
        this._pendingRequest = true;
        const controller = this._newAbortController();

        // Clamp slider to SHAP gene range
        const maxShapGenes = bioResults.top_genes ? bioResults.top_genes.length : 20;
        this._activeMax = maxShapGenes;
        const slider = document.getElementById('topn-slider');
        const topnInput = document.getElementById('topn-input');
        const sliderVal = parseInt(topnInput?.value || slider.value, 10);
        const nTopGenes = Math.min(sliderVal, maxShapGenes);
        // Sync slider UI to actual value used
        slider.max = maxShapGenes;
        slider.value = nTopGenes;
        slider.step = this._getTopNStep(nTopGenes);
        if (topnInput) { topnInput.max = maxShapGenes; topnInput.value = nTopGenes; }

        App.showLoading('Computing SHAP heatmap...');
        try {
            const clusterRowsEl = document.getElementById('cluster-rows');
            const clusterColsEl = document.getElementById('cluster-cols');
            const data = await API.getShapHeatmap(App.state.sessionId, {
                topN: nTopGenes,
                distance: document.getElementById('distance-select').value,
                linkage: document.getElementById('linkage-select').value,
                colorScale: document.getElementById('colorscale-select').value,
                clusterRows: clusterRowsEl ? clusterRowsEl.checked : false,
                clusterCols: clusterColsEl ? clusterColsEl.checked : true,
                signal: controller.signal,
            });
            if (requestToken !== this._requestToken) return;

            // Update panel header to indicate SHAP mode
            const header = document.querySelector('#panel-heatmap .panel-header h1');
            const subtitle = document.querySelector('#panel-heatmap .panel-header .panel-subtitle');
            if (header) header.textContent = 'SHAP Heatmap';
            if (subtitle) subtitle.textContent = `Top ${nTopGenes} genes by ML feature importance (SHAP). Unlike DEG heatmaps, group separation is not guaranteed.`;

            this.plotShapHeatmap(data);
            App.markStepCompleted('heatmap');
        } catch (err) {
            if (requestToken !== this._requestToken) return;
            if (err.name === 'AbortError') return;
            App.showToast('SHAP heatmap failed: ' + err.message, 'error');
        } finally {
            if (this._abortController === controller) this._abortController = null;
            if (requestToken !== this._requestToken) return;
            this._pendingRequest = false;
            App.hideLoading();
        }
    },

    plotShapHeatmap(data) {
        const { z, x, y, groups, color_scale, row_dendrogram, col_dendrogram, shap_values, model } = data;

        const nGenes = y.length;
        const nSamples = x.length;

        // Group mapping
        const groupNames = Object.keys(groups || {});
        const hasGroups = groupNames.length >= 2;
        const sampleGroupMap = {};
        if (hasGroups) {
            groupNames.forEach((name, i) => {
                (groups[name] || []).forEach(sample => {
                    sampleGroupMap[sample] = name;
                });
            });
        }

        const hasRowDendro = row_dendrogram && row_dendrogram.icoord && row_dendrogram.icoord.length > 0;
        const hasColDendro = col_dendrogram && col_dendrogram.icoord && col_dendrogram.icoord.length > 0;

        const traces = [];

        // Layout proportions — dendro → gene names → heatmap + SHAP bar on right
        const rowDendroWidth = hasRowDendro ? 0.07 : 0;
        const geneNameGap = hasRowDendro ? 0.10 : 0;  // space for gene names
        const shapBarWidth = 0.18;  // right side SHAP bar
        const colDendroHeight = hasColDendro ? 0.10 : 0;
        const groupBarHeight = hasGroups ? 0.03 : 0;
        const legendSpace = hasGroups ? 0.04 : 0;
        const gapBetween = 0.02;

        const xRowDendro = [0, rowDendroWidth > 0 ? rowDendroWidth - 0.01 : 0];
        const heatmapLeft = rowDendroWidth + geneNameGap;
        const xHeatmap = [heatmapLeft, 1 - shapBarWidth - gapBetween];
        const xShapBar = [1 - shapBarWidth + 0.01, 1];
        const yHeatmap = [0, 1 - colDendroHeight - groupBarHeight - legendSpace];
        const yGroupBar = [1 - colDendroHeight - groupBarHeight - legendSpace + 0.005, 1 - colDendroHeight - legendSpace - 0.005];
        const yColDendro = [1 - colDendroHeight - legendSpace, 1 - legendSpace];

        const xNums = Array.from({length: nSamples}, (_, i) => i);
        const yNums = Array.from({length: nGenes}, (_, i) => i);

        // =============================================
        // 1. MAIN HEATMAP
        // =============================================
        const hoverTexts = z.map((row, gi) =>
            row.map((val, si) => {
                const group = sampleGroupMap[x[si]] || '';
                const shapStr = shap_values && shap_values[gi] !== undefined ? `\nSHAP: ${shap_values[gi].toFixed(4)}` : '';
                return `<b>${y[gi]}</b><br>${x[si]}${group ? ` (${group})` : ''}<br>Z-score: ${val.toFixed(2)}${shapStr}`;
            })
        );

        traces.push({
            z: z, x: xNums, y: yNums,
            type: 'heatmap',
            colorscale: this._getColorScale(color_scale),
            zmid: 0, zmin: -3, zmax: 3,
            hovertext: hoverTexts,
            hovertemplate: '%{hovertext}<extra></extra>',
            xaxis: 'x', yaxis: 'y',
            showscale: true,
            colorbar: {
                title: { text: 'Z-score', font: { color: '#e2e8f0', size: 10 } },
                tickfont: { color: '#94a3b8', size: 9 },
                thickness: 12, len: 0.3,
                y: 0.15, yanchor: 'middle',
                x: 1.02, xanchor: 'left',
                orientation: 'v',
                outlinewidth: 0,
                tickvals: [-3, 0, 3],
            },
        });

        // =============================================
        // 2. COLUMN DENDROGRAM — single consolidated trace
        // =============================================
        if (hasColDendro) {
            traces.push(this._buildDendroTrace(col_dendrogram, {type: 'col', xaxis: 'x', yaxis: 'y2'}));
        }

        // =============================================
        // 3. ROW DENDROGRAM — single consolidated trace
        // =============================================
        if (hasRowDendro) {
            traces.push(this._buildDendroTrace(row_dendrogram, {type: 'row', xaxis: 'x3', yaxis: 'y6'}));
        }

        // =============================================
        // 4. GROUP COLOR BAR
        // =============================================
        if (hasGroups) {
            const nGroups = groupNames.length;
            const groupZ = [xNums.map(si => {
                const gi = groupNames.indexOf(sampleGroupMap[x[si]]);
                return gi >= 0 ? gi + 0.5 : -0.5;
            })];

            const groupColorscale = [];
            for (let i = 0; i < nGroups; i++) {
                groupColorscale.push([i / nGroups, this.GROUP_COLORS[i % this.GROUP_COLORS.length]]);
                groupColorscale.push([(i + 1) / nGroups, this.GROUP_COLORS[i % this.GROUP_COLORS.length]]);
            }

            traces.push({
                z: groupZ, x: xNums, y: [0],
                type: 'heatmap',
                colorscale: groupColorscale,
                zmin: 0, zmax: nGroups,
                showscale: false,
                xaxis: 'x', yaxis: 'y4',
                hovertemplate: '<b>%{customdata}</b><extra></extra>',
                customdata: [x.map(s => sampleGroupMap[s] || 'Unassigned')],
            });
        }

        // =============================================
        // 5. SHAP BAR PLOT (right side, aligned with genes)
        // Server returns shap_values already reordered to match y (clustered order)
        // =============================================
        if (shap_values && shap_values.length > 0) {
            const shapBars = shap_values;
            const maxShap = Math.max(...shapBars, 0.001);

            traces.push({
                x: shapBars,
                y: yNums,
                type: 'bar',
                orientation: 'h',
                xaxis: 'x5',
                yaxis: 'y',
                marker: {
                    color: shapBars.map(v => {
                        const ratio = v / maxShap;
                        const r = Math.round(59 + ratio * (139 - 59));
                        const g = Math.round(130 + ratio * (92 - 130));
                        const b = 246;
                        return `rgb(${r},${g},${b})`;
                    }),
                },
                hovertemplate: '<b>%{customdata}</b><br>SHAP: %{x:.4f}<extra></extra>',
                customdata: y,
                showlegend: false,
            });
        }

        // =============================================
        // LAYOUT
        // =============================================
        const heightSlider = document.getElementById('heatmap-height-slider');
        const widthSlider = document.getElementById('heatmap-width-slider');
        const sliderHeight = heightSlider ? parseInt(heightSlider.value) : 0;
        const sliderWidth = widthSlider ? parseInt(widthSlider.value) : 0;
        const totalHeight = sliderHeight > 0 ? sliderHeight : Math.max(550, Math.min(2400, nGenes * 18 + 280));
        const sampleFontSize = Math.max(6, Math.min(11, 500 / nSamples));

        const plotAreaHeight = totalHeight * (yHeatmap[1] - yHeatmap[0]) * 0.85;
        const pxPerGene = plotAreaHeight / nGenes;
        const geneFontSize = Math.max(7, Math.min(12, pxPerGene * 0.85));

        // Show all gene labels for SHAP heatmap (usually ≤100 genes)
        const geneTickVals = yNums;
        const geneTickText = y;

        const maxGeneLen = geneTickText.reduce((max, g) => Math.max(max, g.length), 0);
        const leftMargin = Math.max(80, Math.min(160, maxGeneLen * geneFontSize * 0.65 + 10));

        const layout = {
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            font: { family: 'Inter, sans-serif', color: '#94a3b8', size: 10 },
            height: totalHeight, width: sliderWidth > 0 ? sliderWidth : null,
            margin: { l: hasRowDendro ? 20 : leftMargin, r: 80, t: hasGroups ? 68 : 54, b: 60 },

            xaxis: {
                domain: xHeatmap,
                tickvals: xNums, ticktext: x,
                tickangle: -45,
                tickfont: { size: sampleFontSize, color: '#94a3b8' },
                side: 'bottom',
                showgrid: false, zeroline: false,
                range: [-0.5, nSamples - 0.5],
            },

            // Gene names positioned at heatmap left edge (right of dendrogram)
            yaxis: {
                domain: yHeatmap,
                anchor: hasRowDendro ? 'free' : undefined,
                position: hasRowDendro ? heatmapLeft : undefined,
                tickvals: geneTickVals, ticktext: geneTickText,
                tickfont: { size: geneFontSize, color: '#94a3b8' },
                autorange: 'reversed',
                showgrid: false, zeroline: false,
                range: [-0.5, nGenes - 0.5],
            },

            yaxis2: {
                domain: yColDendro,
                showticklabels: false, showgrid: false,
                zeroline: false, showline: false,
            },

            xaxis3: {
                domain: xRowDendro,
                showticklabels: false, showgrid: false,
                zeroline: false, showline: false,
                autorange: true,
            },

            yaxis4: {
                domain: yGroupBar,
                showticklabels: false, showgrid: false,
                zeroline: false, showline: false,
                fixedrange: true,
            },

            // SHAP bar x-axis (right side)
            xaxis5: {
                domain: xShapBar,
                showgrid: false, zeroline: false,
                showticklabels: true,
                tickfont: { size: 8, color: '#64748b' },
                title: { text: 'SHAP', font: { size: 9, color: '#94a3b8' } },
                side: 'bottom',
            },

            // Row dendrogram Y (separate axis, mirrors heatmap y range)
            yaxis6: {
                domain: yHeatmap,
                anchor: 'x3',
                range: [-0.5, nGenes - 0.5],
                autorange: 'reversed',
                showticklabels: false, showgrid: false,
                zeroline: false, showline: false,
            },

            annotations: [],
        };

        // Group legend
        if (hasGroups) {
            const parts = groupNames.map((name, i) => {
                const c = this.GROUP_COLORS[i % this.GROUP_COLORS.length];
                return `<span style="color:${c}">\u25a0</span> ${this._escapeHtml(name)}`;
            });
            layout.annotations.push({
                text: parts.join('&nbsp;&nbsp;&nbsp;'),
                xref: 'paper', yref: 'paper',
                x: (xHeatmap[0] + xHeatmap[1]) / 2,
                y: 1 - legendSpace / 2,
                showarrow: false,
                font: { size: 11, color: '#cbd5e1' },
                xanchor: 'center', yanchor: 'middle',
            });
        }

        layout.annotations.push(this._buildTopLeftTitleAnnotation(`${model || 'ML'} SHAP Top ${nGenes} genes`));

        Plotly.newPlot('heatmap-plot', traces, layout, {
            responsive: true,
            displayModeBar: true,
            modeBarButtonsToRemove: ['lasso2d', 'select2d'],
            displaylogo: false,
            toImageButtonOptions: { format: 'png', height: 1600, width: 2000, scale: 2 },
        }).then(() => this._resizePlot());
    },

    // =============================================
    // DEG HEATMAP MODE
    // =============================================
    async renderDegHeatmap() {
        if (!App.state.sessionId) return App.showToast('No data loaded', 'error');
        if (!App.state.degResults) return App.showToast('Run DEG analysis first', 'error');

        this._isShapMode = false;
        this._isDegMode = true;
        const requestToken = ++this._requestToken;
        this._pendingRequest = true;
        const controller = this._newAbortController();

        const DEG_MAX_GENES = 500;
        this._activeMax = DEG_MAX_GENES;
        const degSlider = document.getElementById('topn-slider');
        const degInput = document.getElementById('topn-input');
        const rawSlider = parseInt(degInput?.value || degSlider.value, 10) || 30;
        const nTopGenes = Math.min(rawSlider, DEG_MAX_GENES);
        // Sync slider UI to DEG range
        degSlider.max = DEG_MAX_GENES;
        degSlider.value = nTopGenes;
        degSlider.step = this._getTopNStep(nTopGenes);
        if (degInput) { degInput.max = DEG_MAX_GENES; degInput.value = nTopGenes; }
        if (rawSlider > DEG_MAX_GENES) {
            App.showToast(`DEG heatmap capped at ${DEG_MAX_GENES} genes`, 'info');
        }

        App.showLoading('Computing DEG heatmap...');
        try {
            const clusterRowsEl = document.getElementById('cluster-rows');
            const clusterColsEl = document.getElementById('cluster-cols');
            const data = await API.getDegHeatmap(App.state.sessionId, {
                topN: nTopGenes,
                distance: document.getElementById('distance-select').value,
                linkage: document.getElementById('linkage-select').value,
                colorScale: document.getElementById('colorscale-select').value,
                clusterRows: clusterRowsEl ? clusterRowsEl.checked : true,
                clusterCols: clusterColsEl ? clusterColsEl.checked : true,
                signal: controller.signal,
            });
            if (requestToken !== this._requestToken) return;

            const header = document.querySelector('#panel-heatmap .panel-header h1');
            const subtitle = document.querySelector('#panel-heatmap .panel-header .panel-subtitle');
            if (header) header.textContent = 'DEG Heatmap';
            const pType = App.state.degResults?.pvalue_type === 'raw' ? 'raw p-value' : 'adjusted p-value';
            if (subtitle) subtitle.textContent = `Top ${nTopGenes} DEGs by ${pType} (${data.method || 'Wilcoxon'} test)`;

            this.plotDegHeatmap(data);
            App.markStepCompleted('heatmap');
        } catch (err) {
            if (requestToken !== this._requestToken) return;
            if (err.name === 'AbortError') return;
            App.showToast('DEG heatmap failed: ' + err.message, 'error');
        } finally {
            if (this._abortController === controller) this._abortController = null;
            if (requestToken !== this._requestToken) return;
            this._pendingRequest = false;
            App.hideLoading();
        }
    },

    plotDegHeatmap(data) {
        const { z, x, y, groups, color_scale, row_dendrogram, col_dendrogram, log2fc_values, neglog10p_values } = data;

        const nGenes = y.length;
        const nSamples = x.length;

        const groupNames = Object.keys(groups || {});
        const hasGroups = groupNames.length >= 2;
        const sampleGroupMap = {};
        if (hasGroups) {
            groupNames.forEach((name, i) => {
                (groups[name] || []).forEach(sample => { sampleGroupMap[sample] = name; });
            });
        }

        const hasRowDendro = row_dendrogram && row_dendrogram.icoord && row_dendrogram.icoord.length > 0;
        const hasColDendro = col_dendrogram && col_dendrogram.icoord && col_dendrogram.icoord.length > 0;

        const traces = [];

        // Layout — same structure as SHAP heatmap, but with log2FC bar on right
        const rowDendroWidth = hasRowDendro ? 0.07 : 0;
        const geneNameGap = hasRowDendro ? 0.10 : 0;
        const fcBarWidth = 0.14;
        const colDendroHeight = hasColDendro ? 0.10 : 0;
        const groupBarHeight = hasGroups ? 0.03 : 0;
        const legendSpace = hasGroups ? 0.04 : 0;
        const gapBetween = 0.02;

        const xRowDendro = [0, rowDendroWidth > 0 ? rowDendroWidth - 0.01 : 0];
        const heatmapLeft = rowDendroWidth + geneNameGap;
        const xHeatmap = [heatmapLeft, 1 - fcBarWidth - gapBetween];
        const xFcBar = [1 - fcBarWidth + 0.01, 1];
        const yHeatmap = [0, 1 - colDendroHeight - groupBarHeight - legendSpace];
        const yGroupBar = [1 - colDendroHeight - groupBarHeight - legendSpace + 0.005, 1 - colDendroHeight - legendSpace - 0.005];
        const yColDendro = [1 - colDendroHeight - legendSpace, 1 - legendSpace];

        const xNums = Array.from({length: nSamples}, (_, i) => i);
        const yNums = Array.from({length: nGenes}, (_, i) => i);

        // 1. MAIN HEATMAP
        const hoverTexts = z.map((row, gi) =>
            row.map((val, si) => {
                const group = sampleGroupMap[x[si]] || '';
                const fc = log2fc_values ? `\nlog2FC: ${log2fc_values[gi].toFixed(3)}` : '';
                return `<b>${y[gi]}</b><br>${x[si]}${group ? ` (${group})` : ''}<br>Z-score: ${val.toFixed(2)}${fc}`;
            })
        );

        traces.push({
            z: z, x: xNums, y: yNums,
            type: 'heatmap',
            colorscale: this._getColorScale(color_scale),
            zmid: 0, zmin: -3, zmax: 3,
            hovertext: hoverTexts,
            hovertemplate: '%{hovertext}<extra></extra>',
            xaxis: 'x', yaxis: 'y',
            showscale: true,
            colorbar: {
                title: { text: 'Z-score', font: { color: '#e2e8f0', size: 10 } },
                tickfont: { color: '#94a3b8', size: 9 },
                thickness: 12, len: 0.3,
                y: 0.15, yanchor: 'middle',
                x: 1.02, xanchor: 'left',
                orientation: 'v',
                outlinewidth: 0,
                tickvals: [-3, 0, 3],
            },
        });

        // 2. COLUMN DENDROGRAM — single consolidated trace
        if (hasColDendro) {
            traces.push(this._buildDendroTrace(col_dendrogram, {type: 'col', xaxis: 'x', yaxis: 'y2'}));
        }

        // 3. ROW DENDROGRAM — single consolidated trace
        if (hasRowDendro) {
            traces.push(this._buildDendroTrace(row_dendrogram, {type: 'row', xaxis: 'x3', yaxis: 'y6'}));
        }

        // 4. GROUP COLOR BAR
        if (hasGroups) {
            const nGroups = groupNames.length;
            const groupZ = [xNums.map(si => {
                const gi = groupNames.indexOf(sampleGroupMap[x[si]]);
                return gi >= 0 ? gi + 0.5 : -0.5;
            })];
            const groupColorscale = [];
            for (let i = 0; i < nGroups; i++) {
                groupColorscale.push([i / nGroups, this.GROUP_COLORS[i % this.GROUP_COLORS.length]]);
                groupColorscale.push([(i + 1) / nGroups, this.GROUP_COLORS[i % this.GROUP_COLORS.length]]);
            }
            traces.push({
                z: groupZ, x: xNums, y: [0], type: 'heatmap',
                colorscale: groupColorscale, zmin: 0, zmax: nGroups, showscale: false,
                xaxis: 'x', yaxis: 'y4',
                hovertemplate: '<b>%{customdata}</b><extra></extra>',
                customdata: [x.map(s => sampleGroupMap[s] || 'Unassigned')],
            });
        }

        // 5. LOG2FC BAR (right side)
        if (log2fc_values && log2fc_values.length > 0) {
            const maxAbsFc = Math.max(...log2fc_values.map(Math.abs), 0.1);
            traces.push({
                x: log2fc_values, y: yNums,
                type: 'bar', orientation: 'h',
                xaxis: 'x5', yaxis: 'y',
                marker: {
                    color: log2fc_values.map(v => v > 0 ? '#ef4444' : '#3b82f6'),
                },
                hovertemplate: '<b>%{customdata}</b><br>log2FC: %{x:.3f}<extra></extra>',
                customdata: y,
                showlegend: false,
            });
        }

        // LAYOUT
        const heightSlider = document.getElementById('heatmap-height-slider');
        const widthSlider = document.getElementById('heatmap-width-slider');
        const sliderHeight = heightSlider ? parseInt(heightSlider.value) : 0;
        const sliderWidth = widthSlider ? parseInt(widthSlider.value) : 0;
        const totalHeight = sliderHeight > 0 ? sliderHeight : Math.max(550, Math.min(2400, nGenes * 18 + 280));
        const sampleFontSize = Math.max(6, Math.min(11, 500 / nSamples));
        const plotAreaHeight = totalHeight * (yHeatmap[1] - yHeatmap[0]) * 0.85;
        const pxPerGene = plotAreaHeight / nGenes;
        const geneFontSize = Math.max(7, Math.min(12, pxPerGene * 0.85));
        const geneTickVals = yNums;
        const geneTickText = y;
        const maxGeneLen = geneTickText.reduce((max, g) => Math.max(max, g.length), 0);
        const leftMargin = Math.max(80, Math.min(160, maxGeneLen * geneFontSize * 0.65 + 10));

        const layout = {
            paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
            font: { family: 'Inter, sans-serif', color: '#94a3b8', size: 10 },
            height: totalHeight, width: sliderWidth > 0 ? sliderWidth : null,
            margin: { l: hasRowDendro ? 20 : leftMargin, r: 80, t: hasGroups ? 68 : 54, b: 60 },

            xaxis: { domain: xHeatmap, tickvals: xNums, ticktext: x, tickangle: -45,
                tickfont: { size: sampleFontSize, color: '#94a3b8' }, side: 'bottom',
                showgrid: false, zeroline: false, range: [-0.5, nSamples - 0.5] },
            yaxis: { domain: yHeatmap,
                anchor: hasRowDendro ? 'free' : undefined, position: hasRowDendro ? heatmapLeft : undefined,
                tickvals: geneTickVals, ticktext: geneTickText,
                tickfont: { size: geneFontSize, color: '#94a3b8' },
                autorange: 'reversed', showgrid: false, zeroline: false, range: [-0.5, nGenes - 0.5] },
            yaxis2: { domain: yColDendro, showticklabels: false, showgrid: false, zeroline: false, showline: false },
            xaxis3: { domain: xRowDendro, showticklabels: false, showgrid: false, zeroline: false, showline: false, autorange: true },
            yaxis4: { domain: yGroupBar, showticklabels: false, showgrid: false, zeroline: false, showline: false, fixedrange: true },
            xaxis5: { domain: xFcBar, showgrid: false, zeroline: true, zerolinecolor: 'rgba(255,255,255,0.2)',
                showticklabels: true, tickfont: { size: 8, color: '#64748b' },
                title: { text: 'log2FC', font: { size: 9, color: '#94a3b8' } }, side: 'bottom' },
            yaxis6: { domain: yHeatmap, anchor: 'x3', range: [-0.5, nGenes - 0.5],
                autorange: 'reversed', showticklabels: false, showgrid: false, zeroline: false, showline: false },

            annotations: [],
        };

        if (hasGroups) {
            const parts = groupNames.map((name, i) => {
                const c = this.GROUP_COLORS[i % this.GROUP_COLORS.length];
                return `<span style="color:${c}">\u25a0</span> ${this._escapeHtml(name)}`;
            });
            layout.annotations.push({ text: parts.join('&nbsp;&nbsp;&nbsp;'),
                xref: 'paper', yref: 'paper', x: (xHeatmap[0] + xHeatmap[1]) / 2, y: 1 - legendSpace / 2,
                showarrow: false, font: { size: 11, color: '#cbd5e1' }, xanchor: 'center', yanchor: 'middle' });
        }

        layout.annotations.push(this._buildTopLeftTitleAnnotation(`DEG Top ${nGenes} genes`));

        Plotly.newPlot('heatmap-plot', traces, layout, {
            responsive: true, displayModeBar: true,
            modeBarButtonsToRemove: ['lasso2d', 'select2d'], displaylogo: false,
            toImageButtonOptions: { format: 'png', height: 1600, width: 2000, scale: 2 },
        }).then(() => this._resizePlot());
    },

    _getColorScale(name) {
        const scales = {
            'RdBu_r': [
                [0, '#2166ac'], [0.1, '#4393c3'], [0.2, '#92c5de'],
                [0.3, '#d1e5f0'], [0.45, '#f7f7f7'], [0.55, '#f7f7f7'],
                [0.7, '#fddbc7'], [0.8, '#f4a582'], [0.9, '#d6604d'],
                [1, '#b2182b'],
            ],
            'RdYlBu_r': [
                [0, '#313695'], [0.1, '#4575b4'], [0.2, '#74add1'],
                [0.3, '#abd9e9'], [0.4, '#e0f3f8'], [0.5, '#ffffbf'],
                [0.6, '#fee090'], [0.7, '#fdae61'], [0.8, '#f46d43'],
                [0.9, '#d73027'], [1, '#a50026'],
            ],
            'Viridis': 'Viridis',
            'Inferno': 'Inferno',
            'Plasma': 'Plasma',
        };
        return scales[name] || scales['RdBu_r'];
    },
};
