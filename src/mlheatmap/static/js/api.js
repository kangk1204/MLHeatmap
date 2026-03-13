/* API Client */
const API = {
    baseUrl: '/api/v1',

    async getCapabilities() {
        const res = await fetch(`${this.baseUrl}/capabilities`);
        if (!res.ok) throw new Error((await res.json()).error || 'Capabilities request failed');
        return res.json();
    },

    async upload(file) {
        const formData = new FormData();
        formData.append('file', file);
        const res = await fetch(`${this.baseUrl}/upload`, { method: 'POST', body: formData });
        if (!res.ok) throw new Error((await res.json()).error || 'Upload failed');
        return res.json();
    },

    async mapGenes(sessionId, species, idType) {
        const res = await fetch(`${this.baseUrl}/gene-mapping`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ session_id: sessionId, species, id_type: idType }),
        });
        if (!res.ok) throw new Error((await res.json()).error || 'Mapping failed');
        return res.json();
    },

    async normalize(sessionId, method) {
        const res = await fetch(`${this.baseUrl}/normalize`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ session_id: sessionId, method }),
        });
        if (!res.ok) throw new Error((await res.json()).error || 'Normalization failed');
        return res.json();
    },

    async getHeatmap(sessionId, opts = {}) {
        const params = new URLSearchParams({
            session_id: sessionId,
            top_n: opts.topN || 500,
            distance: opts.distance || 'correlation',
            linkage: opts.linkage || 'average',
            color_scale: opts.colorScale || 'RdBu_r',
            cluster_rows: opts.clusterRows !== undefined ? opts.clusterRows : true,
            cluster_cols: opts.clusterCols !== undefined ? opts.clusterCols : true,
        });
        const res = await fetch(`${this.baseUrl}/heatmap?${params}`);
        if (!res.ok) throw new Error((await res.json()).error || 'Heatmap failed');
        return res.json();
    },

    async setGroups(sessionId, groups) {
        const res = await fetch(`${this.baseUrl}/groups`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ session_id: sessionId, groups }),
        });
        if (!res.ok) throw new Error((await res.json()).error || 'Groups failed');
        return res.json();
    },

    async excludeSamples(sessionId, samples) {
        const res = await fetch(`${this.baseUrl}/groups/exclude`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ session_id: sessionId, samples }),
        });
        if (!res.ok) throw new Error((await res.json()).error || 'Exclude failed');
        return res.json();
    },

    async includeSamples(sessionId, samples) {
        const res = await fetch(`${this.baseUrl}/groups/include`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ session_id: sessionId, samples }),
        });
        if (!res.ok) throw new Error((await res.json()).error || 'Include failed');
        return res.json();
    },

    biomarkerStream(sessionId, opts = {}) {
        const params = new URLSearchParams({
            session_id: sessionId,
            n_top_genes: opts.nTopGenes || 20,
            n_estimators: opts.nEstimators || 500,
            cv_folds: opts.cvFolds || 5,
            model: opts.model || 'rf',
            panel_method: opts.panelMethod || 'forward',
        });
        return new EventSource(`${this.baseUrl}/biomarker/stream?${params}`);
    },

    async getHeatmapImage(sessionId, opts = {}) {
        const params = new URLSearchParams({
            session_id: sessionId,
            top_n: opts.topN || 500,
            distance: opts.distance || 'correlation',
            linkage: opts.linkage || 'average',
            color_scale: opts.colorScale || 'RdBu_r',
            cluster_rows: opts.clusterRows !== undefined ? opts.clusterRows : true,
            cluster_cols: opts.clusterCols !== undefined ? opts.clusterCols : true,
            fmt: opts.fmt || 'png',
            dpi: opts.dpi || 150,
        });
        const res = await fetch(`${this.baseUrl}/heatmap/render?${params}`);
        if (!res.ok) {
            const err = await res.json();
            throw new Error(err.error || 'Server render failed');
        }
        const blob = await res.blob();
        return URL.createObjectURL(blob);
    },

    async getShapHeatmap(sessionId, opts = {}) {
        const params = new URLSearchParams({
            session_id: sessionId,
            top_n: opts.topN || 20,
            distance: opts.distance || 'correlation',
            linkage: opts.linkage || 'average',
            color_scale: opts.colorScale || 'RdBu_r',
            cluster_rows: opts.clusterRows !== undefined ? opts.clusterRows : true,
            cluster_cols: opts.clusterCols !== undefined ? opts.clusterCols : true,
        });
        const res = await fetch(`${this.baseUrl}/heatmap/shap?${params}`);
        if (!res.ok) throw new Error((await res.json()).error || 'SHAP heatmap failed');
        return res.json();
    },

    async runDeg(sessionId, opts = {}) {
        const params = new URLSearchParams({
            session_id: sessionId,
            method: opts.method || 'wilcoxon',
            log2fc_threshold: opts.log2fcThreshold || 1.0,
            pvalue_threshold: opts.pvalueThreshold || 0.05,
            use_raw_pvalue: opts.useRawPvalue || false,
        });
        if (opts.referenceGroup) {
            params.set('reference_group', opts.referenceGroup);
        }
        const res = await fetch(`${this.baseUrl}/biomarker/deg?${params}`);
        if (!res.ok) throw new Error((await res.json()).error || 'DEG analysis failed');
        return res.json();
    },

    async getDegHeatmap(sessionId, opts = {}) {
        const params = new URLSearchParams({
            session_id: sessionId,
            top_n: opts.topN || 30,
            distance: opts.distance || 'correlation',
            linkage: opts.linkage || 'average',
            color_scale: opts.colorScale || 'RdBu_r',
            cluster_rows: opts.clusterRows !== undefined ? opts.clusterRows : true,
            cluster_cols: opts.clusterCols !== undefined ? opts.clusterCols : true,
        });
        const res = await fetch(`${this.baseUrl}/heatmap/deg?${params}`);
        if (!res.ok) throw new Error((await res.json()).error || 'DEG heatmap failed');
        return res.json();
    },

    exportUrl(sessionId, type, dpi = 300) {
        const params = new URLSearchParams({ session_id: sessionId, type, dpi });
        return `${this.baseUrl}/export?${params}`;
    },
};
