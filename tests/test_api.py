"""API integration tests using FastAPI TestClient."""

import os
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
