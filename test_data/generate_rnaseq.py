import numpy as np
import pandas as pd

np.random.seed(42)

n_genes = 60000
samples = [
    "Control_1", "Control_2", "Control_3", "Control_4",
    "Treatment_1", "Treatment_2", "Treatment_3", "Treatment_4",
]

# Generate Ensembl-style gene IDs (ENSG + 11 digits)
gene_ids = [f"ENSG{i:011d}" for i in range(1, n_genes + 1)]

# Assign each gene to an expression tier
tier = np.random.choice(
    ["zero", "low", "medium", "high"],
    size=n_genes,
    p=[0.15, 0.70, 0.12, 0.03],
)

# Base mean for each gene
base_mean = np.zeros(n_genes)
base_mean[tier == "zero"] = 0
base_mean[tier == "low"] = np.random.exponential(15, size=(tier == "low").sum())
base_mean[tier == "medium"] = np.random.uniform(100, 5000, size=(tier == "medium").sum())
base_mean[tier == "high"] = np.random.uniform(5000, 100000, size=(tier == "high").sum())

# Generate counts using negative binomial (overdispersed Poisson)
# dispersion parameter varies by expression level
counts = np.zeros((n_genes, len(samples)), dtype=int)
for j in range(len(samples)):
    for i in range(n_genes):
        mu = base_mean[i]
        if mu < 1:
            counts[i, j] = 0
            continue
        # Add sample-level noise to the mean
        sample_mu = mu * np.random.uniform(0.8, 1.2)
        # NB parameterisation: n (size), p
        size = max(1, mu / 5)  # higher dispersion for low-count genes
        p = size / (size + sample_mu)
        counts[i, j] = np.random.negative_binomial(size, p)

# Introduce ~200 differentially expressed genes (fold change 2-5x)
de_indices = np.random.choice(
    np.where(tier != "zero")[0], size=200, replace=False
)
for idx in de_indices:
    fc = np.random.uniform(2, 5)
    if np.random.rand() < 0.5:
        # Up in treatment
        for j in range(4, 8):
            counts[idx, j] = int(counts[idx, j] * fc)
    else:
        # Down in treatment
        for j in range(4, 8):
            counts[idx, j] = max(0, int(counts[idx, j] / fc))

# Sprinkle additional zeros for sparsity (~10% of low-count entries)
low_mask = (tier == "low")
for j in range(len(samples)):
    zero_out = np.random.rand(low_mask.sum()) < 0.10
    counts[low_mask, j] = np.where(zero_out, 0, counts[low_mask, j])

df = pd.DataFrame(counts, index=gene_ids, columns=samples)
df.index.name = "gene"

out_path = "/Users/keunsoo/Projects/06_MLHeatmap/test_data/human_60k_8samples.csv"
df.to_csv(out_path)
print(f"Saved {df.shape[0]} genes x {df.shape[1]} samples to {out_path}")
print(f"File size: {pd.io.common.file_exists(out_path)}")
print(f"Non-zero fraction: {(df.values > 0).mean():.3f}")
print(f"DE genes: {len(de_indices)}")
print(df.describe())
