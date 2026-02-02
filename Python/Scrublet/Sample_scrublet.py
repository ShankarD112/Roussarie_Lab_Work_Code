"""Scrublet workflow: generate per-sample counts + doublet calls."""

import scanpy as sc
sc.__version__

import scrublet as scr
scr

import importlib.metadata as md
import pandas as pd
import scipy.sparse as sp

print(md.version("scrublet"))


import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os

# Replace with your Cell Ranger output base directory.
base_dir = "path"

samples = [
    "sample1", "sample2", "sample3", "sample4",
    "sample5", "sample6", "sample7", "sample8",
    "sample9", "sample10", "sample11", "sample12", "sample13"
]

# Where to write per-sample and merged outputs.
out_dir = "path"
os.makedirs(out_dir, exist_ok=True)



adatas = []

for s in samples:
    mtx_dir = os.path.join(base_dir, s, "outs", "filtered_feature_bc_matrix")
    if not os.path.isdir(mtx_dir):
        raise FileNotFoundError(f"Missing directory: {mtx_dir}")

    # Reads Cell Ranger mtx + barcodes + features, returns counts in .X
    ad = sc.read_10x_mtx(mtx_dir, cache=False)

    # Sample label
    ad.obs["orig.ident"] = s

    # Make cell IDs unique across samples (important for merging)
    ad.obs_names = [f"{s}_{bc}" for bc in ad.obs_names]
    ad.obs_names_make_unique()

    # Save per-sample counts AnnData
    out_path = os.path.join(out_dir, f"{s}_counts.h5ad")
    ad.write_h5ad(out_path)

    print(f"{s}: {ad.n_obs} cells x {ad.n_vars} genes -> {out_path}")
    adatas.append(ad)

adata_all = sc.concat(adatas, join="outer", index_unique=None)
merged_out = os.path.join(out_dir, "ALL_counts_only.h5ad")
adata_all.write_h5ad(merged_out)
print("Merged:", adata_all.shape, "->", merged_out)



# Where to read per-sample counts + write scrublet outputs.
in_dir = "path"
out_dir_scrub = "path"
os.makedirs(out_dir_scrub, exist_ok=True)

def expected_rate(n_cells: int) -> float:
    return (n_cells / 1000) * 0.01


for s in samples:
    fp = os.path.join(in_dir, f"{s}_counts.h5ad")
    ad = sc.read_h5ad(fp)

    X = ad.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)

    rate = expected_rate(ad.n_obs)
    print(f"{s}: n_cells={ad.n_obs}, expected_doublet_rate={rate:.4f}")

    scrub = scr.Scrublet(X, expected_doublet_rate=rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    ad.obs["scrublet_score"] = doublet_scores
    ad.obs["scrublet_doublet"] = predicted_doublets
    ad.obs["scrublet_expected_rate"] = rate

    out_fp = os.path.join(out_dir_scrub, f"{s}_counts_scrublet.h5ad")
    ad.write_h5ad(out_fp)

    print(f"Wrote: {out_fp}")




# Where to collect scrublet outputs for summary.
in_dir_scrub = "path"

dfs = []
for s in samples:
    fp = os.path.join(in_dir_scrub, f"{s}_counts_scrublet.h5ad")
    ad = sc.read_h5ad(fp)

    tmp = ad.obs[["orig.ident", "scrublet_score", "scrublet_doublet", "scrublet_expected_rate"]].copy()
    tmp["barcode"] = ad.obs_names
    dfs.append(tmp.reset_index(drop=True))

all_df = pd.concat(dfs, ignore_index=True)
all_df = all_df[["barcode", "orig.ident", "scrublet_score", "scrublet_doublet", "scrublet_expected_rate"]]

all_df.head()


summary = all_df.groupby("orig.ident").agg(
    n_cells=("scrublet_doublet", "size"),
    n_doublets=("scrublet_doublet", "sum"),
)
summary["pct_doublets"] = 100 * summary["n_doublets"] / summary["n_cells"]

summary


out_csv = "path"
all_df.to_csv(out_csv, index=False)
out_csv









