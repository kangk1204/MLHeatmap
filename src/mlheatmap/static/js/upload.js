/* Upload Panel Logic */
const Upload = {
    init() {
        const zone = document.getElementById('upload-zone');
        const input = document.getElementById('file-input');

        zone.addEventListener('click', () => input.click());
        zone.addEventListener('dragover', (e) => {
            e.preventDefault();
            zone.classList.add('drag-over');
        });
        zone.addEventListener('dragleave', () => zone.classList.remove('drag-over'));
        zone.addEventListener('drop', (e) => {
            e.preventDefault();
            zone.classList.remove('drag-over');
            if (e.dataTransfer.files.length) this.handleFile(e.dataTransfer.files[0]);
        });
        input.addEventListener('change', () => {
            if (input.files.length) this.handleFile(input.files[0]);
        });

        document.getElementById('btn-to-mapping').addEventListener('click', () => {
            App.goToPanel('mapping');
        });
    },

    async handleFile(file) {
        App.showLoading('Uploading and parsing...');
        try {
            const result = await API.upload(file);
            App.state.sessionId = result.session_id;
            App.state.sampleNames = result.sample_names;
            App.state.species = result.detected_species;
            App.state.idType = result.detected_id_type;

            // Show results
            document.getElementById('badge-shape').textContent = `${result.shape[0].toLocaleString()} genes × ${result.shape[1]} samples`;
            document.getElementById('badge-species').textContent = `Species: ${result.detected_species}`;
            document.getElementById('badge-idtype').textContent = `ID type: ${result.detected_id_type}`;

            // Build preview table
            const thead = document.querySelector('#preview-table thead');
            const tbody = document.querySelector('#preview-table tbody');
            thead.innerHTML = '';
            tbody.innerHTML = '';

            if (result.preview.length > 0) {
                const cols = Object.keys(result.preview[0]);
                const headerRow = document.createElement('tr');
                cols.forEach(col => {
                    const th = document.createElement('th');
                    th.textContent = col;
                    headerRow.appendChild(th);
                });
                thead.appendChild(headerRow);

                result.preview.forEach(row => {
                    const tr = document.createElement('tr');
                    cols.forEach(col => {
                        const td = document.createElement('td');
                        const val = row[col];
                        td.textContent = typeof val === 'number' ? val.toLocaleString() : val;
                        tr.appendChild(td);
                    });
                    tbody.appendChild(tr);
                });
            }

            // Show filtering info
            if (result.filtering) {
                const fi = result.filtering;
                document.getElementById('filter-before').textContent = fi.before.toLocaleString();
                document.getElementById('filter-after').textContent = fi.after.toLocaleString();
                document.getElementById('filter-removed').textContent = fi.removed.toLocaleString();
                document.getElementById('filter-min-count').textContent = fi.min_count;
                document.getElementById('filter-min-samples').textContent = fi.min_samples;
                document.getElementById('filter-info').classList.remove('hidden');
            }

            document.getElementById('upload-result').classList.remove('hidden');
            document.getElementById('upload-zone').style.display = 'none';
            App.markStepCompleted('upload');
            App.showToast('File uploaded successfully', 'success');

            // Pre-select species radio
            if (result.detected_species !== 'unknown') {
                const radio = document.querySelector(`input[name="species"][value="${result.detected_species}"]`);
                if (radio) radio.checked = true;
            }
        } catch (err) {
            App.showToast(err.message, 'error');
        } finally {
            App.hideLoading();
        }
    },
};
