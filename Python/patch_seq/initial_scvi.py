#!/usr/bin/env python3

import scvi
import scanpy as sc
import matplotlib.pyplot as plt



scvi.settings.seed = 0
print("scvi-tools version:", scvi.__version__)



patch_seq = sc.read_h5ad("path")
snsm = sc.read_h5ad("path")

print(patch_seq)
print(snsm)




# patch_seq: batch from month
patch_seq.obs["batch"] = patch_seq.obs["month"]
patch_seq.obs["region"] = patch_seq.obs.get("region", "Unknown")

# snsm: fixed batch label
snsm.obs["batch"] = "snsm"
snsm.obs["region"] = snsm.obs.get("region", "Unknown")

print(patch_seq.obs["batch"].unique())
print(snsm.obs["batch"].unique())




adata = sc.concat(
    [patch_seq, snsm],
    label="dataset",
    keys=["patch_seq", "snsm"],
    join="outer",
    index_unique=None,
)

print(adata)
print(adata.obs["batch"].unique())
print(adata.obs["region"].unique())




scvi.model.SCVI.setup_anndata(
    adata,
    batch_key="batch",
    labels_key="region",
)


scanvi = scvi.model.SCANVI(adata)
scanvi.train()


adata.obs["predicted_region"] = scanvi.predict()
adata.obs["predicted_region"].value_counts()

adata.obsm["X_scANVI"] = scanvi.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.umap(adata)


sc.pl.umap(adata, color=["dataset", "predicted_region"])
sc.pl.umap(adata, color=["batch"])



if "cell_id" in adata.obs.columns:
    sc.pl.umap(adata, color="cell_id", frameon=False, show=False)

    umap_coords = adata.obsm["X_umap"]
    cell_ids = adata.obs["cell_id"].astype(str)

    ax = plt.gca()
    for i, txt in enumerate(cell_ids):
        if txt in ("nan", "", "NA"):
            continue
        x, y = umap_coords[i]
        ax.text(
            x, y, txt,
            fontsize=6,
            alpha=0.7,
            ha="right",
            va="bottom",
        )
    plt.show()




adata.write("path")
print("Integrated scANVI AnnData saved.")



















