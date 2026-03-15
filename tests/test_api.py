"""API integration tests using FastAPI TestClient."""

import os
import json
import threading
from unittest.mock import patch

import pytest
from fastapi.testclient import TestClient

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def _prepare_deg_session(client):
    """Upload, map, normalize, and split into two groups for DEG tests."""
    path = os.path.join(DATA_DIR, "small_quick_test.csv")
    with open(path, "rb") as f:
        r = client.post("/api/v1/upload", files={"file": ("test.csv", f, "text/csv")})
    assert r.status_code == 200
    sid = r.json()["session_id"]
    samples = r.json()["sample_names"]

    r = client.post("/api/v1/gene-mapping", json={
        "session_id": sid,
        "species": "human",
        "id_type": "auto",
    })
    assert r.status_code == 200

    r = client.post("/api/v1/normalize", json={
        "session_id": sid,
        "method": "log2",
    })
    assert r.status_code == 200

    half = len(samples) // 2
    groups = {"Control": samples[:half], "Case": samples[half:]}
    r = client.post("/api/v1/groups", json={
        "session_id": sid,
        "groups": groups,
    })
    assert r.status_code == 200

    return sid, samples, groups


def _extract_sse_event_payload(text: str, event_name: str) -> dict:
    blocks = text.split("\n\n")
    for block in blocks:
        if f"event: {event_name}" not in block:
            continue
        for line in block.splitlines():
            if line.startswith("data: "):
                return json.loads(line[len("data: "):])
    raise AssertionError(f"SSE event '{event_name}' not found")


@pytest.fixture(scope="module")
def client():
    from mlheatmap.server import create_app
    app = create_app()
    return TestClient(app)


@pytest.fixture(scope="module")
def uploaded_session(client):
    """Upload the small test file and return session_id + sample_names."""
    path = os.path.join(DATA_DIR, "small_quick_test.csv")
    with open(path, "rb") as f:
        r = client.post("/api/v1/upload", files={"file": ("test.csv", f, "text/csv")})
    assert r.status_code == 200
    data = r.json()
    return {
        "session_id": data["session_id"],
        "sample_names": data["sample_names"],
        "shape": data["shape"],
        "species": data["detected_species"],
    }


# ──────────────────────────────────────────────
# Upload Tests
# ──────────────────────────────────────────────

class TestUpload:
    def test_upload_csv(self, client):
        path = os.path.join(DATA_DIR, "small_quick_test.csv")
        with open(path, "rb") as f:
            r = client.post("/api/v1/upload", files={"file": ("test.csv", f, "text/csv")})
        assert r.status_code == 200
        data = r.json()
        assert "session_id" in data
        assert data["shape"][1] == 6  # 6 samples
        assert data["detected_species"] == "human"
        assert data["detected_id_type"] == "ensembl"

    def test_upload_tsv(self, client):
        path = os.path.join(DATA_DIR, "symbol_3groups.tsv")
        with open(path, "rb") as f:
            r = client.post("/api/v1/upload", files={"file": ("test.tsv", f, "text/tab-separated-values")})
        assert r.status_code == 200
        data = r.json()
        assert data["shape"][1] == 9
        assert data["detected_id_type"] == "symbol"

    def test_upload_mouse(self, client):
        path = os.path.join(DATA_DIR, "mouse_ensembl_8samples.csv")
        with open(path, "rb") as f:
            r = client.post("/api/v1/upload", files={"file": ("mouse.csv", f, "text/csv")})
        assert r.status_code == 200
        data = r.json()
        assert data["detected_species"] == "mouse"
        assert data["detected_id_type"] == "ensembl"

    def test_upload_invalid_file(self, client):
        r = client.post("/api/v1/upload", files={"file": ("bad.csv", b"not,valid\ndata", "text/csv")})
        # Should return 400 or 200 with empty result
        assert r.status_code in (200, 400)

    def test_upload_rejects_partial_missing_or_nonnumeric_cells(self, client):
        bad_csv = b"gene_id,S1,S2\nGeneA,10,\nGeneB,5,abc\n"
        r = client.post("/api/v1/upload", files={"file": ("bad.csv", bad_csv, "text/csv")})
        assert r.status_code == 400
        data = r.json()
        assert data["invalid_cell_count"] == 2
        assert "S2" in data["invalid_columns"]

    def test_upload_with_filtering_info(self, client):
        path = os.path.join(DATA_DIR, "human_ensembl_12samples.csv")
        with open(path, "rb") as f:
            r = client.post("/api/v1/upload", files={"file": ("human.csv", f, "text/csv")})
        data = r.json()
        assert "filtering" in data
        assert data["filtering"]["before"] >= data["filtering"]["after"]

    def test_upload_sparse_data(self, client):
        path = os.path.join(DATA_DIR, "sparse_edge_case.csv")
        with open(path, "rb") as f:
            r = client.post("/api/v1/upload", files={"file": ("sparse.csv", f, "text/csv")})
        assert r.status_code == 200
        data = r.json()
        assert data["shape"][0] > 0  # Some genes should survive filtering


# ──────────────────────────────────────────────
# Gene Mapping Tests
# ──────────────────────────────────────────────

class TestGeneMapping:
    def test_map_genes_auto(self, client, uploaded_session):
        sid = uploaded_session["session_id"]
        r = client.post("/api/v1/gene-mapping", json={
            "session_id": sid,
            "species": "human",
            "id_type": "auto",
        })
        assert r.status_code == 200
        data = r.json()
        assert "mapped_count" in data
        assert "unmapped_count" in data
        assert data["mapped_count"] + data["unmapped_count"] == data["total"]

    def test_map_genes_invalid_session(self, client):
        r = client.post("/api/v1/gene-mapping", json={
            "session_id": "nonexistent",
            "species": "human",
            "id_type": "auto",
        })
        assert r.status_code == 404


# ──────────────────────────────────────────────
# Normalize Tests
# ──────────────────────────────────────────────

class TestNormalize:
    @pytest.mark.parametrize("method", ["deseq2", "tpm", "log2"])
    def test_normalize_methods(self, client, uploaded_session, method):
        sid = uploaded_session["session_id"]
        # First map genes
        client.post("/api/v1/gene-mapping", json={
            "session_id": sid, "species": "human", "id_type": "auto",
        })
        r = client.post("/api/v1/normalize", json={
            "session_id": sid,
            "method": method,
        })
        assert r.status_code == 200
        data = r.json()
        assert data["method"] == method
        assert len(data["shape"]) == 2
        assert "stats" in data
        assert data["stats"]["min"] <= data["stats"]["median"] <= data["stats"]["max"]

    def test_normalize_invalid_method(self, client, uploaded_session):
        sid = uploaded_session["session_id"]
        r = client.post("/api/v1/normalize", json={
            "session_id": sid,
            "method": "invalid_method",
        })
        assert r.status_code == 400

    def test_normalize_invalid_session(self, client):
        r = client.post("/api/v1/normalize", json={
            "session_id": "nonexistent",
            "method": "log2",
        })
        assert r.status_code == 404


# ──────────────────────────────────────────────
# Groups Tests
# ──────────────────────────────────────────────

class TestGroups:
    def test_set_and_get_groups(self, client, uploaded_session):
        sid = uploaded_session["session_id"]
        samples = uploaded_session["sample_names"]

        groups = {"Ctrl": samples[:3], "Drug": samples[3:]}
        r = client.post("/api/v1/groups", json={
            "session_id": sid, "groups": groups,
        })
        assert r.status_code == 200
        assert r.json()["n_groups"] == 2

        r = client.get(f"/api/v1/groups?session_id={sid}")
        assert r.status_code == 200
        assert len(r.json()["groups"]) == 2

    def test_exclude_samples(self, client, uploaded_session):
        sid = uploaded_session["session_id"]
        samples = uploaded_session["sample_names"]

        r = client.post("/api/v1/groups/exclude", json={
            "session_id": sid,
            "samples": [samples[0]],
        })
        assert r.status_code == 200
        assert samples[0] in r.json()["excluded"]
        assert samples[0] not in r.json()["remaining"]

    def test_include_samples(self, client, uploaded_session):
        sid = uploaded_session["session_id"]
        samples = uploaded_session["sample_names"]

        # Include back what was excluded
        r = client.post("/api/v1/groups/include", json={
            "session_id": sid,
            "samples": [samples[0]],
        })
        assert r.status_code == 200
        assert samples[0] in r.json()["remaining"]

    def test_groups_invalid_session(self, client):
        r = client.get("/api/v1/groups?session_id=nonexistent")
        assert r.status_code == 404


# ──────────────────────────────────────────────
# Heatmap Tests
# ──────────────────────────────────────────────

class TestHeatmap:
    def test_heatmap_basic(self, client, uploaded_session):
        sid = uploaded_session["session_id"]

        r = client.get(f"/api/v1/heatmap?session_id={sid}&top_n=20")
        assert r.status_code == 200
        data = r.json()
        assert "z" in data
        assert "x" in data
        assert "y" in data
        assert len(data["y"]) <= 20

    def test_heatmap_different_params(self, client, uploaded_session):
        sid = uploaded_session["session_id"]

        r = client.get(f"/api/v1/heatmap?session_id={sid}&top_n=10&distance=euclidean&linkage=ward")
        assert r.status_code == 200
        data = r.json()
        assert len(data["y"]) <= 10

    def test_heatmap_no_data(self, client):
        r = client.get("/api/v1/heatmap?session_id=nonexistent&top_n=50")
        assert r.status_code == 404

    def test_heatmap_without_column_clustering_groups_samples(self, client):
        sid, samples, _ = _prepare_deg_session(client)
        interleaved_groups = {
            "Group A": [samples[0], samples[2], samples[4]],
            "Group B": [samples[1], samples[3], samples[5]],
        }
        r = client.post("/api/v1/groups", json={
            "session_id": sid,
            "groups": interleaved_groups,
        })
        assert r.status_code == 200

        r = client.get(f"/api/v1/heatmap?session_id={sid}&top_n=20&cluster_cols=false")
        assert r.status_code == 200
        data = r.json()
        assert data["x"] == interleaved_groups["Group A"] + interleaved_groups["Group B"]


class TestDEG:
    def test_deg_reference_group_reorders_labels(self, client):
        sid, _, _ = _prepare_deg_session(client)

        r = client.get(f"/api/v1/biomarker/deg?session_id={sid}&reference_group=Control")
        assert r.status_code == 200
        data = r.json()
        assert data["reference_group"] == "Control"
        assert data["comparison_group"] == "Case"

    def test_deg_reference_group_invalid_name(self, client):
        sid, _, _ = _prepare_deg_session(client)

        r = client.get(f"/api/v1/biomarker/deg?session_id={sid}&reference_group=Unknown")
        assert r.status_code == 400
        assert "Unknown reference group" in r.json()["error"]

    def test_deg_reference_group_missing_after_exclusion_returns_400(self, client):
        sid, _, groups = _prepare_deg_session(client)

        r = client.post("/api/v1/groups/exclude", json={
            "session_id": sid,
            "samples": groups["Case"],
        })
        assert r.status_code == 200

        r = client.get(f"/api/v1/biomarker/deg?session_id={sid}&reference_group=Control")
        assert r.status_code == 400
        assert "no valid samples after exclusion" in r.json()["error"]

    def test_deg_response_includes_effect_size_basis(self, client):
        sid, _, _ = _prepare_deg_session(client)

        r = client.get(f"/api/v1/biomarker/deg?session_id={sid}&reference_group=Control")
        assert r.status_code == 200
        data = r.json()
        assert data["effect_size_basis"] == "counts"
        assert data["normalization_method"] == "log2"


class TestReentryAndConcurrency:
    def test_group_change_invalidates_existing_ml_and_deg_results(self, client):
        sid, samples, _ = _prepare_deg_session(client)

        r = client.get(f"/api/v1/biomarker/deg?session_id={sid}&reference_group=Control")
        assert r.status_code == 200

        r = client.get(
            f"/api/v1/biomarker/stream?session_id={sid}&n_top_genes=5&n_estimators=50&cv_folds=2"
        )
        assert r.status_code == 200
        assert "event: complete" in r.text
        complete = _extract_sse_event_payload(r.text, "complete")
        assert complete["optimal_combo"]["evaluation"] == "nested_outer_cv"
        assert "auc_std" in complete["optimal_combo"]
        assert isinstance(complete["optimal_combo"]["selection_frequency"], list)

        new_groups = {
            "Alt A": [samples[0], samples[2], samples[4]],
            "Alt B": [samples[1], samples[3], samples[5]],
        }
        r = client.post("/api/v1/groups", json={"session_id": sid, "groups": new_groups})
        assert r.status_code == 200

        session = client.app.state.sessions.get(sid)
        assert session.biomarker_results is None
        assert session.deg_results is None
        assert session.heatmap_data is None

        r = client.get(f"/api/v1/heatmap/shap?session_id={sid}&top_n=5")
        assert r.status_code == 400
        assert "Run biomarker analysis first" in r.json()["error"]

        r = client.get(f"/api/v1/heatmap/deg?session_id={sid}&top_n=5")
        assert r.status_code == 400
        assert "Run DEG analysis first" in r.json()["error"]

    def test_biomarker_stream_does_not_store_stale_results_after_group_change(self, client):
        sid, samples, _ = _prepare_deg_session(client)
        started = threading.Event()
        finish = threading.Event()
        response = {}

        def fake_run_biomarker_analysis(**kwargs):
            started.set()
            finish.wait(5)
            return {
                "accuracy": 0.9,
                "top_genes": [{"symbol": "G1", "importance": 1.0, "shap_mean_abs": 0.5}],
                "roc_data": [],
                "optimal_combo": None,
                "model": "Random Forest",
            }

        def run_stream():
            r = client.get(f"/api/v1/biomarker/stream?session_id={sid}&n_top_genes=5")
            response["status_code"] = r.status_code
            response["text"] = r.text

        with patch("mlheatmap.core.biomarker.run_biomarker_analysis", fake_run_biomarker_analysis):
            thread = threading.Thread(target=run_stream)
            thread.start()
            assert started.wait(5)

            r = client.post(
                "/api/v1/groups",
                json={
                    "session_id": sid,
                    "groups": {
                        "Alt A": [samples[0], samples[2], samples[4]],
                        "Alt B": [samples[1], samples[3], samples[5]],
                    },
                },
            )
            assert r.status_code == 200

            finish.set()
            thread.join(5)

        session = client.app.state.sessions.get(sid)
        assert session.biomarker_results is None
        assert response["status_code"] == 200
        assert "event: app_error" in response["text"]
        assert "inputs changed during execution" in response["text"]

    def test_heatmap_does_not_cache_stale_results_after_group_change(self, client):
        sid, samples, _ = _prepare_deg_session(client)
        started = threading.Event()
        finish = threading.Event()
        response = {}

        def fake_compute_heatmap_data(**kwargs):
            started.set()
            finish.wait(5)
            return {
                "z": [[1, 2, 3, 4, 5, 6]],
                "x": kwargs["sample_names"],
                "y": ["Gene1"],
                "row_dendrogram": {"icoord": [], "dcoord": []},
                "col_dendrogram": {"icoord": [], "dcoord": []},
                "n_total_genes": 1,
                "n_shown_genes": 1,
            }

        def run_heatmap():
            r = client.get(f"/api/v1/heatmap?session_id={sid}&top_n=20&cluster_cols=false")
            response["status_code"] = r.status_code
            response["json"] = r.json()

        with patch("mlheatmap.core.clustering.compute_heatmap_data", fake_compute_heatmap_data):
            thread = threading.Thread(target=run_heatmap)
            thread.start()
            assert started.wait(5)

            r = client.post(
                "/api/v1/groups",
                json={
                    "session_id": sid,
                    "groups": {
                        "Alt A": [samples[0], samples[2], samples[4]],
                        "Alt B": [samples[1], samples[3], samples[5]],
                    },
                },
            )
            assert r.status_code == 200

            finish.set()
            thread.join(5)

        session = client.app.state.sessions.get(sid)
        assert session.heatmap_data is None
        assert response["status_code"] == 409
        assert "inputs changed during computation" in response["json"]["error"]


# ──────────────────────────────────────────────
# Full Workflow Test (E2E)
# ──────────────────────────────────────────────

class TestFullWorkflow:
    def test_complete_pipeline(self, client):
        """Upload -> Map -> Normalize -> Groups -> Heatmap (full pipeline)."""
        # 1. Upload
        path = os.path.join(DATA_DIR, "human_ensembl_12samples.csv")
        with open(path, "rb") as f:
            r = client.post("/api/v1/upload", files={"file": ("test.csv", f, "text/csv")})
        assert r.status_code == 200
        sid = r.json()["session_id"]
        samples = r.json()["sample_names"]

        # 2. Gene mapping
        r = client.post("/api/v1/gene-mapping", json={
            "session_id": sid, "species": "human", "id_type": "auto",
        })
        assert r.status_code == 200

        # 3. Normalize
        r = client.post("/api/v1/normalize", json={
            "session_id": sid, "method": "deseq2",
        })
        assert r.status_code == 200

        # 4. Set groups
        half = len(samples) // 2
        groups = {"Control": samples[:half], "Treated": samples[half:]}
        r = client.post("/api/v1/groups", json={
            "session_id": sid, "groups": groups,
        })
        assert r.status_code == 200

        # 5. Heatmap
        r = client.get(f"/api/v1/heatmap?session_id={sid}&top_n=50")
        assert r.status_code == 200
        hm = r.json()
        assert len(hm["z"]) > 0
        assert "groups" in hm
        assert "color_scale" in hm

    def test_pipeline_3groups(self, client):
        """Test with 3-group TSV data."""
        path = os.path.join(DATA_DIR, "symbol_3groups.tsv")
        with open(path, "rb") as f:
            r = client.post("/api/v1/upload", files={"file": ("test.tsv", f, "text/tab-separated-values")})
        assert r.status_code == 200
        sid = r.json()["session_id"]
        samples = r.json()["sample_names"]

        # Gene mapping (symbols should mostly map to themselves)
        r = client.post("/api/v1/gene-mapping", json={
            "session_id": sid, "species": "human", "id_type": "symbol",
        })
        assert r.status_code == 200

        # Normalize
        r = client.post("/api/v1/normalize", json={
            "session_id": sid, "method": "log2",
        })
        assert r.status_code == 200

        # Groups (3 groups)
        groups = {
            "Normal": [s for s in samples if s.startswith("Normal")],
            "TumorA": [s for s in samples if s.startswith("TumorA")],
            "TumorB": [s for s in samples if s.startswith("TumorB")],
        }
        r = client.post("/api/v1/groups", json={
            "session_id": sid, "groups": groups,
        })
        assert r.status_code == 200
        assert r.json()["n_groups"] == 3

        # Heatmap
        r = client.get(f"/api/v1/heatmap?session_id={sid}&top_n=30")
        assert r.status_code == 200

    def test_pipeline_mouse(self, client):
        """Test with mouse Ensembl data."""
        path = os.path.join(DATA_DIR, "mouse_ensembl_8samples.csv")
        with open(path, "rb") as f:
            r = client.post("/api/v1/upload", files={"file": ("mouse.csv", f, "text/csv")})
        assert r.status_code == 200
        sid = r.json()["session_id"]
        samples = r.json()["sample_names"]
        assert r.json()["detected_species"] == "mouse"

        r = client.post("/api/v1/gene-mapping", json={
            "session_id": sid, "species": "mouse", "id_type": "auto",
        })
        assert r.status_code == 200

        r = client.post("/api/v1/normalize", json={
            "session_id": sid, "method": "tpm",
        })
        assert r.status_code == 200

        groups = {
            "WT": [s for s in samples if s.startswith("WT")],
            "KO": [s for s in samples if s.startswith("KO")],
        }
        r = client.post("/api/v1/groups", json={
            "session_id": sid, "groups": groups,
        })
        assert r.status_code == 200

        r = client.get(f"/api/v1/heatmap?session_id={sid}&top_n=50")
        assert r.status_code == 200


# ──────────────────────────────────────────────
# Static Files Tests
# ──────────────────────────────────────────────

class TestStaticFiles:
    def test_index_page(self, client):
        r = client.get("/")
        assert r.status_code == 200
        assert "MLHeatmap" in r.text

    def test_css(self, client):
        r = client.get("/static/css/main.css")
        assert r.status_code == 200

    def test_js_files(self, client):
        for js in ["app.js", "api.js", "upload.js", "groups.js",
                    "heatmap.js", "biomarker.js", "export.js"]:
            r = client.get(f"/static/js/{js}")
            assert r.status_code == 200, f"Failed to load {js}"
