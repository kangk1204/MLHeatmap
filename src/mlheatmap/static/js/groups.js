/* Groups Panel Logic - Enhanced for bulk operations */
const Groups = {
    groupCount: 0,
    groupColors: ['#3b82f6', '#ef4444', '#10b981', '#f59e0b', '#8b5cf6', '#ec4899',
                  '#06b6d4', '#84cc16', '#f97316', '#6366f1'],
    draggedSample: null,
    selectedSamples: new Set(),
    lastClickedIdx: -1,  // for shift-click range select
    allSamples: [],       // ordered list of all sample names

    init() {
        document.getElementById('btn-add-group').addEventListener('click', () => this.addGroup());
        document.getElementById('btn-auto-group').addEventListener('click', () => this.autoDetectGroups());
        document.getElementById('btn-to-heatmap').addEventListener('click', () => this.saveAndContinue());

        // Search filter
        const searchInput = document.getElementById('sample-search');
        if (searchInput) {
            searchInput.addEventListener('input', (e) => this.filterSamples(e.target.value));
        }

        // Select all / deselect all
        const selectAllBtn = document.getElementById('btn-select-all');
        if (selectAllBtn) {
            selectAllBtn.addEventListener('click', () => this.selectAllVisible());
        }
        const deselectBtn = document.getElementById('btn-deselect-all');
        if (deselectBtn) {
            deselectBtn.addEventListener('click', () => this.deselectAll());
        }

        // Regex auto-group
        const regexBtn = document.getElementById('btn-regex-group');
        if (regexBtn) {
            regexBtn.addEventListener('click', () => this.regexGroupDialog());
        }
    },

    populate(sampleNames) {
        this.allSamples = [...sampleNames];
        this.selectedSamples.clear();
        this.lastClickedIdx = -1;

        const pool = document.getElementById('sample-pool');
        pool.innerHTML = '';
        sampleNames.forEach((name, i) => {
            pool.appendChild(this.createChip(name, true, i));
        });
        this.updatePoolCount();

        // Auto-create 2 groups
        if (this.groupCount === 0) {
            this.addGroup('Group A');
            this.addGroup('Group B');
        }

        // Clear search
        const searchInput = document.getElementById('sample-search');
        if (searchInput) searchInput.value = '';
    },

    createChip(name, showRemove = false, poolIdx = -1) {
        const chip = document.createElement('div');
        chip.className = 'sample-chip';
        chip.draggable = true;
        chip.dataset.sample = name;
        if (poolIdx >= 0) chip.dataset.poolIdx = poolIdx;

        const label = document.createElement('span');
        label.className = 'chip-label';
        label.textContent = name;
        chip.appendChild(label);

        if (showRemove) {
            const btn = document.createElement('span');
            btn.className = 'remove-btn';
            btn.textContent = '\u00d7';
            btn.addEventListener('click', (e) => {
                e.stopPropagation();
                this.toggleExclude(chip);
            });
            chip.appendChild(btn);
        }

        // Click to select/deselect (only in pool)
        chip.addEventListener('click', (e) => {
            if (e.target.classList.contains('remove-btn')) return;
            const inPool = chip.closest('#sample-pool');
            if (!inPool) return;

            const idx = parseInt(chip.dataset.poolIdx, 10);

            if (e.shiftKey && this.lastClickedIdx >= 0) {
                // Shift-click: range select
                const poolChips = [...document.querySelectorAll('#sample-pool .sample-chip:not(.filter-hidden)')];
                const indices = poolChips.map(c => parseInt(c.dataset.poolIdx, 10));
                const lastPos = indices.indexOf(this.lastClickedIdx);
                const curPos = indices.indexOf(idx);
                if (lastPos >= 0 && curPos >= 0) {
                    const start = Math.min(lastPos, curPos);
                    const end = Math.max(lastPos, curPos);
                    for (let i = start; i <= end; i++) {
                        const sName = poolChips[i].dataset.sample;
                        this.selectedSamples.add(sName);
                        poolChips[i].classList.add('selected');
                    }
                }
            } else if (e.ctrlKey || e.metaKey) {
                // Ctrl/Cmd-click: toggle
                if (this.selectedSamples.has(name)) {
                    this.selectedSamples.delete(name);
                    chip.classList.remove('selected');
                } else {
                    this.selectedSamples.add(name);
                    chip.classList.add('selected');
                }
            } else {
                // Normal click: toggle single
                if (this.selectedSamples.has(name)) {
                    this.selectedSamples.delete(name);
                    chip.classList.remove('selected');
                } else {
                    this.selectedSamples.add(name);
                    chip.classList.add('selected');
                }
            }
            this.lastClickedIdx = idx;
            this.updateSelectionCount();
        });

        // Drag events
        chip.addEventListener('dragstart', (e) => {
            // If dragging a selected chip, drag all selected
            if (this.selectedSamples.has(name) && this.selectedSamples.size > 1) {
                e.dataTransfer.setData('text/plain', JSON.stringify([...this.selectedSamples]));
                e.dataTransfer.effectAllowed = 'move';
            } else {
                e.dataTransfer.setData('text/plain', JSON.stringify([name]));
                e.dataTransfer.effectAllowed = 'move';
            }
            this.draggedSample = name;
            chip.classList.add('dragging');
        });

        chip.addEventListener('dragend', () => {
            chip.classList.remove('dragging');
            this.draggedSample = null;
        });

        // Mark as selected if in selection set
        if (this.selectedSamples.has(name)) {
            chip.classList.add('selected');
        }

        return chip;
    },

    filterSamples(query) {
        const lower = query.toLowerCase();
        const chips = document.querySelectorAll('#sample-pool .sample-chip');
        chips.forEach(chip => {
            const name = chip.dataset.sample.toLowerCase();
            if (name.includes(lower)) {
                chip.classList.remove('filter-hidden');
            } else {
                chip.classList.add('filter-hidden');
            }
        });
    },

    selectAllVisible() {
        const chips = document.querySelectorAll('#sample-pool .sample-chip:not(.filter-hidden):not(.excluded)');
        chips.forEach(chip => {
            this.selectedSamples.add(chip.dataset.sample);
            chip.classList.add('selected');
        });
        this.updateSelectionCount();
    },

    deselectAll() {
        this.selectedSamples.clear();
        document.querySelectorAll('#sample-pool .sample-chip.selected').forEach(c => c.classList.remove('selected'));
        this.updateSelectionCount();
    },

    updateSelectionCount() {
        const badge = document.getElementById('selection-count');
        if (badge) {
            const n = this.selectedSamples.size;
            badge.textContent = n > 0 ? `${n} selected` : '';
            badge.classList.toggle('hidden', n === 0);
        }
    },

    addGroup(name) {
        const idx = this.groupCount++;
        const color = this.groupColors[idx % this.groupColors.length];
        if (!name) name = `Group ${String.fromCharCode(65 + idx)}`;

        const col = document.createElement('div');
        col.className = 'group-column';
        col.dataset.groupIdx = idx;

        const escapedName = name.replace(/&/g, '&amp;').replace(/"/g, '&quot;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
        col.innerHTML = `
            <div class="group-header">
                <div class="group-color" style="background:${color}"></div>
                <input class="group-name-input" value="${escapedName}" data-group-idx="${idx}">
                <span class="group-sample-count" data-group-idx="${idx}">0</span>
                <button class="group-delete-btn" data-group-idx="${idx}">\u00d7</button>
            </div>
            <div class="group-assign-bar">
                <button class="btn btn-sm btn-assign" data-group-idx="${idx}" title="Assign selected samples to this group">
                    <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5"><path d="M5 12h14M12 5l7 7-7 7"/></svg>
                    Assign Selected
                </button>
            </div>
            <div class="group-dropzone" data-group-idx="${idx}"></div>
        `;

        // Assign selected button
        col.querySelector('.btn-assign').addEventListener('click', () => {
            this.assignSelectedToGroup(idx);
        });

        const dropzone = col.querySelector('.group-dropzone');
        dropzone.addEventListener('dragover', (e) => {
            e.preventDefault();
            dropzone.classList.add('drag-over');
        });
        dropzone.addEventListener('dragleave', () => dropzone.classList.remove('drag-over'));
        dropzone.addEventListener('drop', (e) => {
            e.preventDefault();
            dropzone.classList.remove('drag-over');
            try {
                const names = JSON.parse(e.dataTransfer.getData('text/plain'));
                names.forEach(n => this.moveSampleToGroup(n, idx));
                // Clear selection after drop
                this.selectedSamples.clear();
                document.querySelectorAll('.sample-chip.selected').forEach(c => c.classList.remove('selected'));
                this.updateSelectionCount();
            } catch {
                const sampleName = e.dataTransfer.getData('text/plain');
                this.moveSampleToGroup(sampleName, idx);
            }
        });

        col.querySelector('.group-delete-btn').addEventListener('click', () => {
            // Move samples back to pool
            const chips = dropzone.querySelectorAll('.sample-chip');
            const pool = document.getElementById('sample-pool');
            chips.forEach(c => {
                const sName = c.dataset.sample;
                const origIdx = this.allSamples.indexOf(sName);
                pool.appendChild(this.createChip(sName, true, origIdx >= 0 ? origIdx : 0));
            });
            col.remove();
            this.updatePoolCount();
        });

        document.getElementById('groups-area').appendChild(col);
    },

    assignSelectedToGroup(groupIdx) {
        if (this.selectedSamples.size === 0) {
            App.showToast('Select samples first (click or Shift+click)', 'info');
            return;
        }
        const names = [...this.selectedSamples];
        names.forEach(n => this.moveSampleToGroup(n, groupIdx));
        this.selectedSamples.clear();
        document.querySelectorAll('.sample-chip.selected').forEach(c => c.classList.remove('selected'));
        this.updateSelectionCount();
        App.showToast(`${names.length} sample${names.length > 1 ? 's' : ''} assigned`, 'success');
    },

    moveSampleToGroup(sampleName, groupIdx) {
        // Remove from current location
        document.querySelectorAll(`.sample-chip[data-sample="${CSS.escape(sampleName)}"]`).forEach(el => el.remove());

        // Add to group
        const dropzone = document.querySelector(`.group-dropzone[data-group-idx="${groupIdx}"]`);
        if (dropzone) {
            dropzone.appendChild(this.createChip(sampleName, true));
        }
        this.updatePoolCount();
        this.updateGroupCounts();
    },

    toggleExclude(chip) {
        chip.classList.toggle('excluded');
        const name = chip.dataset.sample;
        if (chip.classList.contains('excluded')) {
            this.selectedSamples.delete(name);
            chip.classList.remove('selected');
            if (App.state.sessionId) {
                API.excludeSamples(App.state.sessionId, [name]);
            }
        } else {
            if (App.state.sessionId) {
                API.includeSamples(App.state.sessionId, [name]);
            }
        }
        this.updateSelectionCount();
    },

    autoDetectGroups() {
        const samples = App.state.sampleNames || [];
        if (samples.length === 0) return;

        // Try to find common prefix patterns
        const prefixes = {};
        samples.forEach(s => {
            let prefix = s;
            for (const sep of ['_', '-', '.', ' ']) {
                const parts = s.split(sep);
                if (parts.length >= 2) {
                    prefix = parts[0];
                    break;
                }
            }
            if (!prefixes[prefix]) prefixes[prefix] = [];
            prefixes[prefix].push(s);
        });

        const groups = Object.entries(prefixes);
        if (groups.length >= 2 && groups.length <= 10) {
            // Clear existing
            document.getElementById('groups-area').innerHTML = '';
            this.groupCount = 0;
            document.getElementById('sample-pool').innerHTML = '';
            this.selectedSamples.clear();

            groups.forEach(([prefix, members]) => {
                this.addGroup(prefix);
                const idx = this.groupCount - 1;
                members.forEach(s => this.moveSampleToGroup(s, idx));
            });

            App.showToast(`Auto-detected ${groups.length} groups`, 'success');
        } else {
            App.showToast('Could not auto-detect groups. Please assign manually.', 'info');
            this.populate(samples);
        }
    },

    regexGroupDialog() {
        const samples = App.state.sampleNames || [];
        if (samples.length === 0) return;

        const pattern = prompt(
            'Enter a regex pattern with a capture group to extract the group name.\n' +
            'Example: ^(\\w+?)_\\d+ will group "Control_1" → "Control"\n' +
            'Example: ^([A-Za-z]+) will group by leading letters',
            '^(\\w+?)[-_]'
        );

        if (!pattern) return;

        try {
            const re = new RegExp(pattern);
            const prefixes = {};
            let matched = 0;

            samples.forEach(s => {
                const m = s.match(re);
                if (m && m[1]) {
                    const key = m[1];
                    if (!prefixes[key]) prefixes[key] = [];
                    prefixes[key].push(s);
                    matched++;
                }
            });

            const groups = Object.entries(prefixes);
            if (groups.length >= 2 && groups.length <= 10 && matched > samples.length * 0.5) {
                document.getElementById('groups-area').innerHTML = '';
                this.groupCount = 0;
                document.getElementById('sample-pool').innerHTML = '';
                this.selectedSamples.clear();

                groups.forEach(([prefix, members]) => {
                    this.addGroup(prefix);
                    const idx = this.groupCount - 1;
                    members.forEach(s => this.moveSampleToGroup(s, idx));
                });

                App.showToast(`Regex grouped ${matched} samples into ${groups.length} groups`, 'success');
            } else {
                App.showToast(`Pattern matched ${matched}/${samples.length} samples into ${groups.length} groups. Needs 2-10 groups matching >50% of samples.`, 'info');
            }
        } catch (e) {
            App.showToast('Invalid regex pattern: ' + e.message, 'error');
        }
    },

    updatePoolCount() {
        const pool = document.getElementById('sample-pool');
        const count = pool.querySelectorAll('.sample-chip:not(.excluded)').length;
        document.getElementById('pool-count').textContent = count;
    },

    updateGroupCounts() {
        document.querySelectorAll('.group-column').forEach(col => {
            const idx = col.dataset.groupIdx;
            const count = col.querySelectorAll('.group-dropzone .sample-chip:not(.excluded)').length;
            const badge = col.querySelector('.group-sample-count');
            if (badge) badge.textContent = count;
        });
    },

    getGroups() {
        const groups = {};
        document.querySelectorAll('.group-column').forEach(col => {
            const name = col.querySelector('.group-name-input').value.trim();
            const samples = [];
            col.querySelectorAll('.group-dropzone .sample-chip:not(.excluded)').forEach(chip => {
                samples.push(chip.dataset.sample);
            });
            if (name && samples.length > 0) {
                groups[name] = samples;
            }
        });
        return groups;
    },

    async saveAndContinue() {
        const groups = this.getGroups();
        const groupNames = Object.keys(groups);

        if (groupNames.length < 2) {
            App.showToast('Please assign samples to at least 2 groups', 'error');
            return;
        }

        App.showLoading('Saving groups...');
        try {
            await API.setGroups(App.state.sessionId, groups);
            App.state.groups = groups;
            App.markStepCompleted('groups');
            App.goToPanel('heatmap');

            // Auto-render heatmap
            setTimeout(() => Heatmap.render(), 300);
        } catch (err) {
            App.showToast(err.message, 'error');
        } finally {
            App.hideLoading();
        }
    },
};
