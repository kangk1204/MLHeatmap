"""Tests for the public TCGA CRC CMS example builder."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def test_select_primary_tumor_columns_keeps_one_primary_aliquot_per_sample():
    from mlheatmap.public_crc_cms import select_primary_tumor_columns

    frame = pd.DataFrame(
        [[1, 2, 3, 4]],
        index=["ENSG000001"],
        columns=[
            "TCGA-AA-1111-01A",
            "TCGA-AA-1111-01B",
            "TCGA-BB-2222-11A",
            "TCGA-CC-3333-01A",
        ],
    )

    trimmed = select_primary_tumor_columns(frame)

    assert trimmed.columns.tolist() == ["TCGA-AA-1111-01", "TCGA-CC-3333-01"]
    assert trimmed.iloc[0].tolist() == [1, 4]


def test_xena_log2_to_counts_recovers_integer_count_scale():
    from mlheatmap.public_crc_cms import xena_log2_to_counts

    frame = pd.DataFrame(
        [[np.log2(1), np.log2(11)], [np.log2(101), np.log2(2)]],
        index=["GeneA", "GeneB"],
        columns=["Sample1", "Sample2"],
    )

    counts = xena_log2_to_counts(frame)

    assert counts.loc["GeneA", "Sample1"] == 0
    assert counts.loc["GeneA", "Sample2"] == 10
    assert counts.loc["GeneB", "Sample1"] == 100
    assert counts.loc["GeneB", "Sample2"] == 1


def test_apply_gold_labels_prefixes_columns_for_group_autodetect():
    from mlheatmap.public_crc_cms import apply_gold_labels

    frame = pd.DataFrame(
        {
            "TCGA-AA-1111-01": [10, 20],
            "TCGA-BB-2222-01": [30, 40],
            "TCGA-CC-3333-01": [50, 60],
        },
        index=["GAS1", "SPOCK1"],
    )
    labels = pd.DataFrame(
        {
            "short_barcode": ["TCGA-AA-1111", "TCGA-CC-3333"],
            "CMS_gold": ["CMS1", "CMS4"],
        }
    )

    labeled, metadata = apply_gold_labels(frame, labels)

    assert labeled.columns.tolist() == ["CMS1__AA111101", "CMS4__CC333301"]
    assert metadata.to_dict(orient="records") == [
        {
            "tcga_barcode": "TCGA-AA-1111-01",
            "short_barcode": "TCGA-AA-1111",
            "CMS_gold": "CMS1",
            "sample_id": "CMS1__AA111101",
        },
        {
            "tcga_barcode": "TCGA-CC-3333-01",
            "short_barcode": "TCGA-CC-3333",
            "CMS_gold": "CMS4",
            "sample_id": "CMS4__CC333301",
        },
    ]


def test_readme_documents_public_crc_cms_rebuild_command():
    readme = (Path(__file__).parents[1] / "README.md").read_text(encoding="utf-8")

    assert "mlheatmap-download-crc-cms" in readme
    assert "scripts/prepare_public_crc_cms_gold.py" in readme
    assert "tcga_crc_cms_gold_counts.tsv.gz" in readme
