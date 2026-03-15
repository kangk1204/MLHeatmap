#!/usr/bin/env python3
"""Reproduce paper-facing MLHeatmap results without the web UI."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from mlheatmap.api.session import Session
from mlheatmap.core.biomarker import run_biomarker_analysis
from mlheatmap.core.deg import compute_deg
from mlheatmap.core.export import build_results_metadata, export_results_excel
from mlheatmap.core.gene_mapping import detect_id_type, map_gene_ids
from mlheatmap.core.input_io import filter_low_expression, load_count_matrix_path, strict_numeric_matrix
from mlheatmap.core.normalization import deseq2_normalize, log2_normalize, tpm_abundance, tpm_normalize


def _git_commit() -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT,
            text=True,
        ).strip()
    except Exception:
        return None


def _load_groups(path: Path) -> dict[str, list[str]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("Groups JSON must be an object mapping group name to sample-name list.")
    groups = {}
    for group_name, samples in data.items():
        if not isinstance(group_name, str) or not isinstance(samples, list):
            raise ValueError("Each group entry must be {\"group\": [\"sample\", ...]}.")
        groups[group_name] = [str(sample) for sample in samples]
    return groups


def _set_group_indices(sample_names: list[str], groups: dict[str, list[str]]) -> dict[str, list[int]]:
    sample_to_idx = {sample: idx for idx, sample in enumerate(sample_names)}
    indexed = {}
    for group_name, samples in groups.items():
        missing = [sample for sample in samples if sample not in sample_to_idx]
        if missing:
            raise ValueError(f"Group '{group_name}' references unknown samples: {', '.join(missing)}")
        indexed[group_name] = [sample_to_idx[sample] for sample in samples]
    return indexed


def _map_counts(session: Session, species: str, id_type: str) -> None:
    gene_ids = session.raw_counts.index.astype(str).tolist()
    detected_species, detected_id_type = detect_id_type(gene_ids)
    if species == "auto":
        species = detected_species
    if id_type == "auto":
        id_type = detected_id_type
    if species == "unknown":
        raise ValueError("Unable to determine species automatically. Pass --species explicitly.")

    mapping, unmapped = map_gene_ids(gene_ids, species, id_type)
    if mapping:
        df = session.raw_counts.copy()
        df.index = df.index.astype(str)
        df["_symbol"] = df.index.map(lambda gene_id: mapping.get(str(gene_id)))
        df = df.dropna(subset=["_symbol"])
        df = df.groupby("_symbol").sum()
        df.index.name = None
        session.mapped_counts = df
        session.gene_names = df.index.tolist()
        session.sample_names = df.columns.tolist()
    else:
        session.mapped_counts = session.raw_counts.copy()
        session.gene_names = session.raw_counts.index.astype(str).tolist()
        session.sample_names = session.raw_counts.columns.tolist()

    session.species = species
    session.id_type = id_type
    session.metadata["mapping"] = {
        "species": species,
        "id_type": id_type,
        "mapped_count": len(mapping),
        "unmapped_count": len(unmapped),
        "total": len(gene_ids),
    }


def _normalize_session(session: Session, method: str) -> None:
    df = session.mapped_counts if session.mapped_counts is not None else session.raw_counts
    counts = df.values.astype(float)
    if method == "deseq2":
        normalized, size_factors = deseq2_normalize(counts, return_size_factors=True)
        session.size_factors = size_factors
        session.deg_abundance = counts / size_factors[np.newaxis, :]
        session.deg_effect_size_basis = "size_factor_normalized_counts"
    elif method == "tpm":
        normalized = tpm_normalize(counts)
        session.size_factors = None
        session.deg_abundance = tpm_abundance(counts)
        session.deg_effect_size_basis = "tpm"
    elif method == "log2":
        normalized = log2_normalize(counts)
        session.size_factors = None
        session.deg_abundance = counts.copy()
        session.deg_effect_size_basis = "counts"
    else:
        raise ValueError(f"Unknown normalization method: {method}")

    session.normalized = normalized
    session.norm_method = method
    session.gene_names = df.index.tolist()
    session.sample_names = df.columns.tolist()
    session.metadata["normalization"] = {
        "method": method,
        "shape": [int(normalized.shape[0]), int(normalized.shape[1])],
        "effect_size_basis": session.deg_effect_size_basis,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Reproduce paper-facing MLHeatmap outputs.")
    parser.add_argument("--input", required=True, help="Path to count matrix (.csv/.tsv/.xlsx)")
    parser.add_argument("--groups-json", required=True, help="Path to JSON mapping group name to sample-name list")
    parser.add_argument("--output-dir", required=True, help="Directory for biomarker/DEG/metadata outputs")
    parser.add_argument("--species", default="auto", help="Species for mapping (human, mouse, auto)")
    parser.add_argument("--id-type", default="auto", help="ID type for mapping (auto, ensembl, entrez, refseq, symbol)")
    parser.add_argument("--normalize", default="deseq2", choices=["deseq2", "tpm", "log2"])
    parser.add_argument("--model", default="rf")
    parser.add_argument("--panel-method", default="forward", choices=["forward", "lasso", "stability", "mrmr"])
    parser.add_argument("--n-top-genes", type=int, default=20)
    parser.add_argument("--n-estimators", type=int, default=500)
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--deg-method", default="wilcoxon", choices=["wilcoxon", "ttest"])
    parser.add_argument("--deg-log2fc-threshold", type=float, default=1.0)
    parser.add_argument("--deg-pvalue-threshold", type=float, default=0.05)
    parser.add_argument("--deg-use-raw-pvalue", action="store_true")
    parser.add_argument("--reference-group", default="", help="Optional reference group for DEG directionality")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    session = Session()
    raw_df = strict_numeric_matrix(load_count_matrix_path(args.input))
    raw_df = raw_df.loc[(raw_df != 0).any(axis=1)]
    raw_df, filtering = filter_low_expression(raw_df)
    session.raw_counts = raw_df
    session.sample_names = raw_df.columns.tolist()
    session.gene_names = raw_df.index.astype(str).tolist()
    detected_species, detected_id_type = detect_id_type(session.gene_names)
    session.species = detected_species
    session.id_type = detected_id_type
    session.metadata["upload"] = {
        "filename": Path(args.input).name,
        "shape": [int(raw_df.shape[0]), int(raw_df.shape[1])],
        "detected_species": detected_species,
        "detected_id_type": detected_id_type,
        "filtering": filtering,
    }

    _map_counts(session, args.species, args.id_type)
    _normalize_session(session, args.normalize)

    groups = _load_groups(Path(args.groups_json))
    group_indices = _set_group_indices(session.sample_names, groups)
    session.groups = groups
    session.metadata["reproducibility"] = {
        "input": str(Path(args.input).resolve()),
        "groups_json": str(Path(args.groups_json).resolve()),
        "git_commit": _git_commit(),
    }

    biomarker = run_biomarker_analysis(
        expression=session.normalized,
        gene_names=session.gene_names,
        sample_groups=group_indices,
        n_top_genes=args.n_top_genes,
        n_estimators=args.n_estimators,
        cv_folds=args.cv_folds,
        model=args.model,
        panel_method=args.panel_method,
    )
    session.biomarker_results = biomarker
    session.metadata["biomarker"] = {
        "model": args.model,
        "panel_method": args.panel_method,
        "n_top_genes": args.n_top_genes,
        "n_estimators": args.n_estimators,
        "cv_folds_requested": args.cv_folds,
        "roc_evaluation": biomarker.get("roc_evaluation"),
        "panel_evaluation": biomarker.get("optimal_combo", {}).get("evaluation"),
    }

    deg_payload: dict[str, object]
    if len(groups) == 2:
        if args.reference_group:
            if args.reference_group not in group_indices:
                raise ValueError(f"Unknown reference group: {args.reference_group}")
            comparison = next(group for group in group_indices if group != args.reference_group)
            group_indices = {
                comparison: group_indices[comparison],
                args.reference_group: group_indices[args.reference_group],
            }
        deg_payload = compute_deg(
            expression=session.normalized,
            gene_names=session.gene_names,
            sample_groups=group_indices,
            method=args.deg_method,
            log2fc_threshold=args.deg_log2fc_threshold,
            pvalue_threshold=args.deg_pvalue_threshold,
            use_raw_pvalue=args.deg_use_raw_pvalue,
            effect_size_data=session.deg_abundance,
            effect_size_basis=session.deg_effect_size_basis,
        )
        session.deg_results = deg_payload
        session.metadata["deg"] = {
            "method": args.deg_method,
            "pvalue_type": deg_payload.get("pvalue_type", "fdr"),
            "log2fc_threshold": args.deg_log2fc_threshold,
            "pvalue_threshold": args.deg_pvalue_threshold,
            "reference_group": deg_payload.get("reference_group", ""),
            "comparison_group": deg_payload.get("comparison_group", ""),
            "effect_size_basis": deg_payload.get("effect_size_basis", session.deg_effect_size_basis),
            "normalization_method": session.norm_method,
        }
    else:
        deg_payload = {
            "skipped": True,
            "reason": "DEG reproduction requires exactly 2 groups.",
        }

    metadata = build_results_metadata(session)
    (output_dir / "biomarker.json").write_text(json.dumps(biomarker, indent=2), encoding="utf-8")
    (output_dir / "deg.json").write_text(json.dumps(deg_payload, indent=2), encoding="utf-8")
    (output_dir / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    (output_dir / "results.xlsx").write_bytes(export_results_excel(session))

    print(json.dumps({"output_dir": str(output_dir.resolve()), "files": sorted(p.name for p in output_dir.iterdir())}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
