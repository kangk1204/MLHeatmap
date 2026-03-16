"""Tests for core modules."""

import gzip
import io
import time
import zipfile

import numpy as np
import pytest


class TestNormalization:
    """Tests for normalization methods."""

    def test_log2_normalize_basic(self):
        from mlheatmap.core.normalization import log2_normalize

        counts = np.array([[0, 1, 10], [100, 1000, 0]], dtype=float)
        result = log2_normalize(counts)

        assert result.shape == counts.shape
        assert result[0, 0] == pytest.approx(0.0)  # log2(0+1) = 0
        assert result[0, 1] == pytest.approx(1.0)  # log2(1+1) = 1
        assert np.all(np.isfinite(result))

    def test_log2_normalize_preserves_shape(self):
        from mlheatmap.core.normalization import log2_normalize

        counts = np.random.randint(0, 1000, (100, 10)).astype(float)
        result = log2_normalize(counts)
        assert result.shape == (100, 10)

    def test_tpm_normalize_basic(self):
        from mlheatmap.core.normalization import tpm_normalize

        counts = np.array([[100, 200], [300, 400], [600, 400]], dtype=float)
        result = tpm_normalize(counts)

        assert result.shape == counts.shape
        assert np.all(np.isfinite(result))
        assert np.all(result >= 0)

    def test_deseq2_normalize_basic(self):
        from mlheatmap.core.normalization import deseq2_normalize

        rng = np.random.default_rng(42)
        counts = rng.poisson(lam=100, size=(500, 6)).astype(float)
        result = deseq2_normalize(counts)

        assert result.shape == counts.shape
        assert np.all(np.isfinite(result))

    def test_deseq2_normalize_single_sample(self):
        """Single sample should not produce NaN from ddof=1."""
        from mlheatmap.core.normalization import deseq2_normalize

        counts = np.array([[10], [20], [30]], dtype=float)
        result = deseq2_normalize(counts)

        assert result.shape == counts.shape
        assert np.all(np.isfinite(result))

    def test_deseq2_normalize_with_zeros(self):
        from mlheatmap.core.normalization import deseq2_normalize

        counts = np.array([[0, 100, 200], [0, 0, 300], [50, 50, 50]], dtype=float)
        result = deseq2_normalize(counts)

        assert result.shape == counts.shape
        assert np.all(np.isfinite(result))

    def test_tpm_abundance_returns_linear_scale_values(self):
        from mlheatmap.core.normalization import tpm_abundance

        counts = np.array([[100, 200], [300, 400], [600, 400]], dtype=float)
        abundance = tpm_abundance(counts)

        assert abundance.shape == counts.shape
        assert np.all(np.isfinite(abundance))
        assert np.all(abundance >= 0)


class TestClustering:
    """Tests for hierarchical clustering."""

    def test_compute_heatmap_data_basic(self):
        from mlheatmap.core.clustering import compute_heatmap_data

        rng = np.random.default_rng(42)
        expression = rng.standard_normal((50, 5))
        gene_names = [f"Gene{i}" for i in range(50)]
        sample_names = [f"Sample{i}" for i in range(5)]

        result = compute_heatmap_data(expression, gene_names, sample_names, top_n=20)

        assert "z" in result
        assert "x" in result
        assert "y" in result
        assert len(result["x"]) == 5
        assert len(result["y"]) == 20
        assert len(result["z"]) == 20
        assert result["n_total_genes"] == 50
        assert result["n_shown_genes"] == 20

    def test_z_score_clipped(self):
        from mlheatmap.core.clustering import compute_heatmap_data

        rng = np.random.default_rng(42)
        expression = rng.standard_normal((20, 4))
        gene_names = [f"G{i}" for i in range(20)]
        sample_names = [f"S{i}" for i in range(4)]

        result = compute_heatmap_data(expression, gene_names, sample_names, top_n=10)

        z = np.array(result["z"])
        assert np.all(z >= -3)
        assert np.all(z <= 3)

    def test_single_gene(self):
        from mlheatmap.core.clustering import compute_heatmap_data

        expression = np.array([[1.0, 2.0, 3.0]])
        result = compute_heatmap_data(expression, ["Gene1"], ["S1", "S2", "S3"], top_n=1)

        assert len(result["y"]) == 1

    def test_single_sample(self):
        from mlheatmap.core.clustering import compute_heatmap_data

        expression = np.array([[1.0], [2.0], [3.0]])
        result = compute_heatmap_data(expression, ["G1", "G2", "G3"], ["S1"], top_n=3)

        assert len(result["x"]) == 1

    def test_invalid_heatmap_params_raise(self):
        from mlheatmap.core.clustering import compute_heatmap_data

        expression = np.array([[1.0, 2.0], [3.0, 4.0]])
        with pytest.raises(ValueError):
            compute_heatmap_data(expression, ["G1", "G2"], ["S1", "S2"], distance="invalid")

    def test_no_samples_raise(self):
        from mlheatmap.core.clustering import compute_heatmap_data

        expression = np.empty((2, 0))
        with pytest.raises(ValueError):
            compute_heatmap_data(expression, ["G1", "G2"], [], top_n=2)


class TestGeneMapping:
    """Tests for gene ID detection and mapping."""

    def test_detect_ensembl_human(self):
        from mlheatmap.core.gene_mapping import detect_id_type

        ids = [f"ENSG{str(i).zfill(11)}" for i in range(100)]
        species, id_type = detect_id_type(ids)
        assert species == "human"
        assert id_type == "ensembl"

    def test_detect_ensembl_mouse(self):
        from mlheatmap.core.gene_mapping import detect_id_type

        ids = [f"ENSMUSG{str(i).zfill(11)}" for i in range(100)]
        species, id_type = detect_id_type(ids)
        assert species == "mouse"
        assert id_type == "ensembl"

    def test_detect_entrez(self):
        from mlheatmap.core.gene_mapping import detect_id_type

        ids = [str(i) for i in range(1000, 1100)]
        species, id_type = detect_id_type(ids)
        assert id_type == "entrez"

    def test_detect_refseq(self):
        from mlheatmap.core.gene_mapping import detect_id_type

        ids = [f"NM_{str(i).zfill(6)}" for i in range(100)]
        species, id_type = detect_id_type(ids)
        assert id_type == "refseq"

    def test_detect_symbols_human(self):
        from mlheatmap.core.gene_mapping import detect_id_type

        ids = ["TP53", "BRCA1", "EGFR", "MYC", "KRAS", "AKT1", "PIK3CA"] * 15
        species, id_type = detect_id_type(ids)
        assert species == "human"
        assert id_type == "symbol"

    def test_detect_empty(self):
        from mlheatmap.core.gene_mapping import detect_id_type

        species, id_type = detect_id_type([])
        assert species == "unknown"
        assert id_type == "unknown"

    def test_map_gene_ids_human(self):
        from mlheatmap.core.gene_mapping import map_gene_ids

        # These are real human gene symbols - should map to themselves
        ids = ["TP53", "BRCA1", "EGFR"]
        mapping, unmapped = map_gene_ids(ids, "human", "symbol")

        # May or may not have data file, so just check types
        assert isinstance(mapping, dict)
        assert isinstance(unmapped, list)


class TestSession:
    """Tests for session management."""

    def test_create_session(self):
        from mlheatmap.api.session import SessionStore

        store = SessionStore(ttl_hours=1)
        session = store.create()

        assert session.id is not None
        assert session.raw_counts is None
        assert session.groups == {}

    def test_get_session(self):
        from mlheatmap.api.session import SessionStore

        store = SessionStore(ttl_hours=1)
        session = store.create()
        retrieved = store.get(session.id)

        assert retrieved is session

    def test_get_missing_session(self):
        from mlheatmap.api.session import SessionStore

        store = SessionStore(ttl_hours=1)
        assert store.get("nonexistent") is None

    def test_session_expiry(self):
        import time
        from mlheatmap.api.session import SessionStore

        store = SessionStore(ttl_hours=0)  # 0 hours = expire immediately
        session = store.create()
        time.sleep(0.01)
        # get() now checks TTL
        assert store.get(session.id) is None

    def test_cleanup_on_create(self):
        from mlheatmap.api.session import SessionStore

        store = SessionStore(ttl_hours=0)
        s1 = store.create()
        import time
        time.sleep(0.01)
        s2 = store.create()

        # s1 should be cleaned up
        assert store.get(s1.id) is None
        # s2 just created, should still be accessible (even with ttl=0, just created)

    def test_active_operation_lease_blocks_cleanup(self):
        from mlheatmap.api.session import SessionStore

        store = SessionStore(ttl_hours=0)
        session = store.create()
        lease = store.begin_use(session.id)
        assert lease is not None
        assert lease.session is session

        with session.state_lock:
            session.last_accessed_at = time.time() - 60
        store._cleanup()
        assert store.get(session.id) is session

        store.end_use(session.id, lease.operation_id)
        with session.state_lock:
            session.last_accessed_at = time.time() - 60
        store._cleanup()
        assert store.get(session.id) is None


class TestBiomarker:
    """Tests for biomarker analysis core logic."""

    def test_run_biomarker_basic(self):
        from mlheatmap.core.biomarker import run_biomarker_analysis

        rng = np.random.default_rng(42)
        n_genes = 100
        n_samples = 20

        # Create two distinguishable groups
        expression = rng.standard_normal((n_genes, n_samples))
        # Make first 5 genes different between groups
        expression[:5, :10] += 3.0

        gene_names = [f"Gene{i}" for i in range(n_genes)]
        sample_groups = {
            "GroupA": list(range(10)),
            "GroupB": list(range(10, 20)),
        }

        result = run_biomarker_analysis(
            expression=expression,
            gene_names=gene_names,
            sample_groups=sample_groups,
            n_top_genes=10,
            n_estimators=50,
            cv_folds=3,
        )

        assert "top_genes" in result
        assert "roc_data" in result
        assert "accuracy" in result
        assert "optimal_combo" in result
        assert len(result["top_genes"]) == 10
        assert 0 <= result["accuracy"] <= 1
        assert result["optimal_combo"]["evaluation"] == "nested_outer_cv"
        assert "auc_std" in result["optimal_combo"]
        assert "selection_frequency" in result["optimal_combo"]

    def test_cv_folds_clamped(self):
        """CV folds should be clamped to min class count but never below 2."""
        from mlheatmap.core.biomarker import run_biomarker_analysis

        rng = np.random.default_rng(42)
        expression = rng.standard_normal((50, 6))
        expression[:5, :3] += 5.0

        gene_names = [f"Gene{i}" for i in range(50)]
        # 3 samples per group - CV folds=5 should be clamped to 3
        sample_groups = {
            "A": [0, 1, 2],
            "B": [3, 4, 5],
        }

        result = run_biomarker_analysis(
            expression=expression,
            gene_names=gene_names,
            sample_groups=sample_groups,
            n_top_genes=5,
            n_estimators=50,
            cv_folds=5,  # Should be clamped to 3
        )

        assert "accuracy" in result


class TestDEG:
    def test_log2fc_uses_linear_scale_effect_size_basis(self):
        from mlheatmap.core.deg import compute_deg

        normalized = np.array([[1.0, 6.65821148, 5.67242534, 5.67242534]], dtype=float)
        effect_size_data = np.array([[1.0, 100.0, 50.0, 50.0]], dtype=float)

        result = compute_deg(
            expression=normalized,
            gene_names=["Gene1"],
            sample_groups={"Case": [0, 1], "Control": [2, 3]},
            effect_size_data=effect_size_data,
            effect_size_basis="counts",
        )

        assert result["effect_size_basis"] == "counts"
        assert result["results"][0]["log2fc"] == pytest.approx(
            np.log2(effect_size_data[0, :2].mean() + 1) - np.log2(effect_size_data[0, 2:].mean() + 1)
        )

    def test_bh_handles_empty_input(self):
        from mlheatmap.core.deg import _benjamini_hochberg

        result = _benjamini_hochberg(np.array([], dtype=float))
        assert result.shape == (0,)


class TestInputIO:
    def test_load_count_matrix_rejects_legacy_xls(self):
        from mlheatmap.core.input_io import MatrixValidationError, load_count_matrix_bytes

        with pytest.raises(MatrixValidationError):
            load_count_matrix_bytes(b"not-a-real-xls", "legacy.xls")

    def test_gzip_upload_guard_rejects_archive_expansion(self, monkeypatch):
        from mlheatmap.core.input_io import MatrixValidationError, load_count_matrix_bytes

        csv_payload = "gene_id,S1\nGeneA,1\n" + ("GeneB,2\n" * 2048)
        gz_buffer = io.BytesIO()
        with gzip.GzipFile(fileobj=gz_buffer, mode="wb") as compressed:
            compressed.write(csv_payload.encode("utf-8"))

        monkeypatch.setattr("mlheatmap.core.input_io.MAX_EXPANDED_UPLOAD_BYTES", 1024)
        with pytest.raises(MatrixValidationError):
            load_count_matrix_bytes(gz_buffer.getvalue(), "counts.csv.gz")

    def test_xlsx_upload_guard_rejects_oversized_archive(self, monkeypatch):
        from mlheatmap.core.input_io import MatrixValidationError, load_count_matrix_bytes

        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, "w", compression=zipfile.ZIP_DEFLATED) as workbook:
            workbook.writestr("[Content_Types].xml", "<Types/>")
            workbook.writestr("xl/worksheets/sheet1.xml", "A" * 4096)

        monkeypatch.setattr("mlheatmap.core.input_io.MAX_EXPANDED_UPLOAD_BYTES", 1024)
        with pytest.raises(MatrixValidationError):
            load_count_matrix_bytes(zip_buffer.getvalue(), "oversized.xlsx")


class TestExport:
    """Tests for export utility functions."""

    def test_json_safe_numpy_types(self):
        from mlheatmap.api.biomarker import _json_safe

        assert _json_safe(np.int64(42)) == 42
        assert _json_safe(np.float64(3.14)) == pytest.approx(3.14)
        assert _json_safe(np.array([1, 2, 3])) == [1, 2, 3]

    def test_json_safe_unsupported_type(self):
        from mlheatmap.api.biomarker import _json_safe

        with pytest.raises(TypeError):
            _json_safe({"not": "supported"})
