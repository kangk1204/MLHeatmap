/* API Client */
const API = {
    baseUrl: '/api/v1',

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

    exportUrl(sessionId, type, dpi = 300) {
        const params = new URLSearchParams({ session_id: sessionId, type, dpi });
        return `${this.baseUrl}/export?${params}`;
    },
};
