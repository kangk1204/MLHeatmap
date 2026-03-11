"""Generate test input data files for manual and automated testing."""

import csv
import os
import random

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "data")


def generate_ensembl_human_csv(path, n_genes=300, n_samples=12):
    """Ensembl Human gene IDs, Control vs Treated, CSV format."""
    random.seed(42)
    samples = [f"Control_{i+1}" for i in range(n_samples // 2)] + \
              [f"Treated_{i+1}" for i in range(n_samples // 2)]

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene_id"] + samples)
        for i in range(n_genes):
            gene_id = f"ENSG{str(100000 + i).zfill(11)}"
            # First 20 genes are differentially expressed
            if i < 20:
                ctrl = [random.randint(50, 300) for _ in range(n_samples // 2)]
                treat = [random.randint(500, 3000) for _ in range(n_samples // 2)]
            else:
                base = random.randint(10, 500)
                ctrl = [random.randint(max(1, base - 100), base + 100) for _ in range(n_samples // 2)]
                treat = [random.randint(max(1, base - 100), base + 100) for _ in range(n_samples // 2)]
            writer.writerow([gene_id] + ctrl + treat)
    print(f"  Created: {path} ({n_genes} genes x {n_samples} samples)")


def generate_symbol_tsv(path, n_genes=200, n_samples=9):
    """Gene symbol IDs, 3-group comparison, TSV format."""
    random.seed(123)
    symbols = [
        "TP53", "BRCA1", "EGFR", "MYC", "KRAS", "AKT1", "PIK3CA", "PTEN",
        "RB1", "BRAF", "NRAS", "CDH1", "CTNNB1", "APC", "SMAD4", "VHL",
        "WT1", "NF1", "NF2", "TSC1", "TSC2", "MTOR", "JAK2", "FLT3",
    ]
    # Pad with synthetic symbols
    for i in range(n_genes - len(symbols)):
        symbols.append(f"GENE{i+1}")

    samples = [f"Normal_{i+1}" for i in range(3)] + \
              [f"TumorA_{i+1}" for i in range(3)] + \
              [f"TumorB_{i+1}" for i in range(3)]

    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["gene_id"] + samples)
        for i, sym in enumerate(symbols[:n_genes]):
            if i < 8:  # DE genes between groups
                normal = [random.randint(100, 400) for _ in range(3)]
                tumorA = [random.randint(800, 2000) for _ in range(3)]
                tumorB = [random.randint(400, 900) for _ in range(3)]
            else:
                base = random.randint(50, 500)
                normal = [random.randint(max(1, base - 80), base + 80) for _ in range(3)]
                tumorA = [random.randint(max(1, base - 80), base + 80) for _ in range(3)]
                tumorB = [random.randint(max(1, base - 80), base + 80) for _ in range(3)]
            writer.writerow([sym] + normal + tumorA + tumorB)
    print(f"  Created: {path} ({n_genes} genes x {n_samples} samples)")


def generate_mouse_ensembl_csv(path, n_genes=150, n_samples=8):
    """Mouse Ensembl IDs, WT vs KO, CSV format."""
    random.seed(77)
    samples = [f"WT_{i+1}" for i in range(n_samples // 2)] + \
              [f"KO_{i+1}" for i in range(n_samples // 2)]

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene_id"] + samples)
        for i in range(n_genes):
            gene_id = f"ENSMUSG{str(50000 + i).zfill(11)}"
            if i < 10:
                wt = [random.randint(200, 600) for _ in range(n_samples // 2)]
                ko = [random.randint(20, 100) for _ in range(n_samples // 2)]
            else:
                base = random.randint(30, 400)
                wt = [random.randint(max(1, base - 80), base + 80) for _ in range(n_samples // 2)]
                ko = [random.randint(max(1, base - 80), base + 80) for _ in range(n_samples // 2)]
            writer.writerow([gene_id] + wt + ko)
    print(f"  Created: {path} ({n_genes} genes x {n_samples} samples)")


def generate_small_csv(path):
    """Minimal dataset for quick testing (few genes, few samples)."""
    random.seed(99)
    samples = ["Ctrl_1", "Ctrl_2", "Ctrl_3", "Drug_1", "Drug_2", "Drug_3"]
    genes = [f"ENSG{str(200000 + i).zfill(11)}" for i in range(30)]

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene_id"] + samples)
        for i, g in enumerate(genes):
            if i < 5:
                ctrl = [random.randint(100, 300) for _ in range(3)]
                drug = [random.randint(1000, 3000) for _ in range(3)]
            else:
                base = random.randint(50, 500)
                ctrl = [random.randint(max(1, base - 50), base + 50) for _ in range(3)]
                drug = [random.randint(max(1, base - 50), base + 50) for _ in range(3)]
            writer.writerow([g] + ctrl + drug)
    print(f"  Created: {path} (30 genes x 6 samples, minimal)")


def generate_edge_case_csv(path):
    """Edge case: sparse data with many zeros."""
    random.seed(55)
    samples = ["S1", "S2", "S3", "S4", "S5", "S6"]
    genes = [f"ENSG{str(300000 + i).zfill(11)}" for i in range(50)]

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene_id"] + samples)
        for i, g in enumerate(genes):
            row = []
            for j in range(6):
                if random.random() < 0.3:  # 30% zeros
                    row.append(0)
                else:
                    row.append(random.randint(10, 500))
            writer.writerow([g] + row)
    print(f"  Created: {path} (50 genes x 6 samples, sparse)")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print("Generating test data files...\n")

    generate_ensembl_human_csv(os.path.join(OUTPUT_DIR, "human_ensembl_12samples.csv"))
    generate_symbol_tsv(os.path.join(OUTPUT_DIR, "symbol_3groups.tsv"))
    generate_mouse_ensembl_csv(os.path.join(OUTPUT_DIR, "mouse_ensembl_8samples.csv"))
    generate_small_csv(os.path.join(OUTPUT_DIR, "small_quick_test.csv"))
    generate_edge_case_csv(os.path.join(OUTPUT_DIR, "sparse_edge_case.csv"))

    print(f"\nAll files saved to: {OUTPUT_DIR}/")
    print("\nUsage:")
    print("  1. UI test: drag any file into http://127.0.0.1:8765")
    print("  2. Unit test: .venv/bin/python -m pytest tests/ -v")
    print("  3. API test: .venv/bin/python -m pytest tests/test_api.py -v")


if __name__ == "__main__":
    main()
