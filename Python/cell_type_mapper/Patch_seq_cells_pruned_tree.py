#!/usr/bin/env python3
"""
GitHub-ready script derived from a Jupyter notebook.

This script demonstrates a MapMyCells / cell_type_mapper workflow:
1) Load ABC Atlas cache metadata
2) Select a training subset of cells (e.g., isocortex/HPF/OLF feature matrices)
3) Precompute stats (cell_type_mapper)
4) Build reference markers
5) Build query markers
6) Map a query AnnData (.h5ad) to the reference

All filesystem paths are placeholders ("path") for GitHub.
"""

from __future__ import annotations

import json
import os
import pathlib
from dataclasses import dataclass
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd

# Optional imports used by some steps:
#   pip install abc_atlas_access cell_type_mapper anndata scipy mygene h5py
try:
    from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
except ImportError as e:
    AbcProjectCache = None  # type: ignore

# cell_type_mapper CLIs (installed as a Python package)
try:
    from cell_type_mapper.cli.precompute_stats_abc import PrecomputationABCRunner
    from cell_type_mapper.cli.reference_markers import ReferenceMarkerRunner
    from cell_type_mapper.cli.query_markers import QueryMarkerRunner
    from cell_type_mapper.cli.from_specified_markers import FromSpecifiedMarkersRunner
except ImportError:
    PrecomputationABCRunner = None  # type: ignore
    ReferenceMarkerRunner = None  # type: ignore
    QueryMarkerRunner = None  # type: ignore
    FromSpecifiedMarkersRunner = None  # type: ignore


# ----------------------------
# CONFIG
# ----------------------------

@dataclass(frozen=True)
class Config:
    # ABC atlas cache directory
    abc_data_dir: pathlib.Path = pathlib.Path("path")

    # Scratch/work directories
    scratch_dir: pathlib.Path = pathlib.Path("path")
    precompute_dir: pathlib.Path = pathlib.Path("path")
    reference_dir: pathlib.Path = pathlib.Path("path")
    query_dir: pathlib.Path = pathlib.Path("path")

    # Output paths for marker + mapping
    query_marker_path: pathlib.Path = pathlib.Path("path")
    mapping_json_out: pathlib.Path = pathlib.Path("path")
    mapping_csv_out: pathlib.Path = pathlib.Path("path")

    # Input query AnnData (cells x genes)
    query_h5ad_path: pathlib.Path = pathlib.Path("path")

    # Optional: paths to pre-downloaded “mouse_markers” and “precomputed_stats”
    # (If you use the ABC cache method below, you may not need these.)
    marker_json_path: pathlib.Path = pathlib.Path("path")
    precomputed_stats_path: pathlib.Path = pathlib.Path("path")

    # Compute settings
    n_processors: int = 4
    max_gb: int = 10

    # Taxonomy hierarchy (example)
    hierarchy: Tuple[str, ...] = (
        "CCN20230722_CLAS",
        "CCN20230722_SUBC",
        "CCN20230722_SUPT",
        "CCN20230722_CLUS",
    )


# ----------------------------
# UTILITIES
# ----------------------------

def set_thread_env_vars() -> None:
    """Avoid oversubscription in BLAS/NumExpr stacks."""
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"


def ensure_dirs(*dirs: pathlib.Path) -> None:
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)


def directory_from_feature_matrix_label(feature_matrix_label: str) -> str:
    """
    Convert a feature_matrix_label into the ABC cache directory string.
    Example input might look like: 'WMB-10X-...'
    """
    params = feature_matrix_label.split("-")
    return f"{params[0]}-{params[1]}"


# ----------------------------
# CORE WORKFLOW
# ----------------------------

def load_abc_cache(cfg: Config) -> "AbcProjectCache":
    if AbcProjectCache is None:
        raise RuntimeError(
            "abc_atlas_access is not installed. Install it before running this script."
        )

    abc_cache = AbcProjectCache.from_cache_dir(cfg.abc_data_dir)
    abc_cache.load_latest_manifest()
    return abc_cache


def select_training_cells(abc_cache: "AbcProjectCache") -> pd.DataFrame:
    """
    Example: filter cell metadata to only include Isocortex/HPF/OLF feature matrices.
    """
    cell_metadata = abc_cache.get_metadata_dataframe(
        directory="WMB-10X",
        file_name="cell_metadata",
    )

    training_cells = cell_metadata[
        cell_metadata.feature_matrix_label.str.contains("Isocortex")
        | cell_metadata.feature_matrix_label.str.contains("HPF")
        | cell_metadata.feature_matrix_label.str.contains("OLF")
    ].copy()

    return training_cells


def build_h5ad_path_list(
    abc_cache: "AbcProjectCache", training_cells: pd.DataFrame
) -> List[str]:
    """
    For each unique feature matrix in training_cells, collect the ABC data path to raw h5ad.
    """
    h5ad_path_list: List[str] = []
    for feature_matrix in set(training_cells.feature_matrix_label.values):
        directory = directory_from_feature_matrix_label(feature_matrix)
        h5ad_path = abc_cache.get_data_path(
            directory=directory, file_name=f"{feature_matrix}/raw"
        )
        h5ad_path_list.append(str(h5ad_path))
    return h5ad_path_list


def run_precompute_reference_query_markers(
    cfg: Config, abc_cache: "AbcProjectCache", training_cells: pd.DataFrame, h5ad_path_list: List[str]
) -> None:
    """
    Runs:
      1) precompute_stats_abc
      2) reference_markers
      3) query_markers
    """
    if PrecomputationABCRunner is None:
        raise RuntimeError("cell_type_mapper is not installed.")

    # Save training cell metadata (csv)
    training_set_path = cfg.scratch_dir / "training_cell_metadata.csv"
    training_cells.to_csv(training_set_path, index=False)

    # Taxonomy metadata
    cluster_annotation_path = abc_cache.get_metadata_path(
        directory="WMB-taxonomy", file_name="cluster_annotation_term"
    )
    cluster_membership_path = abc_cache.get_metadata_path(
        directory="WMB-taxonomy", file_name="cluster_to_cluster_annotation_membership"
    )

    precompute_out = cfg.precompute_dir / "precomputed_stats.h5"
    precompute_config = {
        "output_path": str(precompute_out),
        "hierarchy": list(cfg.hierarchy),
        "h5ad_path_list": h5ad_path_list,
        "cell_metadata_path": str(training_set_path),
        "cluster_annotation_path": str(cluster_annotation_path),
        "cluster_membership_path": str(cluster_membership_path),
        "n_processors": cfg.n_processors,
        "split_by_dataset": True,
        "do_pruning": True,
        "tmp_dir": str(cfg.scratch_dir),
        "clobber": True,
    }

    print("Running precomputation...")
    PrecomputationABCRunner(args=[], input_data=precompute_config).run()

    # Reference markers
    from cell_type_mapper.cli.reference_markers import ReferenceMarkerRunner

    precomputed_path_list = [
        str(p)
        for p in cfg.precompute_dir.iterdir()
        if p.is_file() and not p.name.endswith("combined.h5")
    ]

    reference_config = {
        "precomputed_path_list": precomputed_path_list,
        "output_dir": str(cfg.reference_dir),
        "tmp_dir": str(cfg.scratch_dir),
        "max_gb": cfg.max_gb,
        "n_processors": cfg.n_processors,
        "clobber": True,
    }

    print("Running reference markers...")
    ReferenceMarkerRunner(args=[], input_data=reference_config).run()

    # Query markers
    from cell_type_mapper.cli.query_markers import QueryMarkerRunner

    reference_marker_path_list = [str(p) for p in cfg.reference_dir.iterdir() if p.is_file()]
    query_marker_config = {
        "output_path": str(cfg.query_marker_path),
        "reference_marker_path_list": reference_marker_path_list,
        "n_processors": cfg.n_processors,
        "tmp_dir": str(cfg.scratch_dir),
    }

    print("Running query markers...")
    QueryMarkerRunner(args=[], input_data=query_marker_config).run()


def run_mapping(cfg: Config) -> None:
    """Run FromSpecifiedMarkersRunner using precomputed stats + query markers."""
    if FromSpecifiedMarkersRunner is None:
        raise RuntimeError("cell_type_mapper is not installed.")

    # Depending on how you ran precompute, this combined file may exist:
    # (If not, replace with the correct path you generated.)
    precompute_combined_path = cfg.precompute_dir / "precomputed_stats.combined.h5"

    mapping_config = {
        "query_path": str(cfg.query_h5ad_path),
        "extended_result_path": str(cfg.mapping_json_out),
        "csv_result_path": str(cfg.mapping_csv_out),
        "tmp_dir": str(cfg.scratch_dir),
        "max_gb": cfg.max_gb,
        "cloud_safe": False,
        "verbose_stdout": True,
        "type_assignment": {
            "normalization": "raw",
            "n_processors": cfg.n_processors,
            "chunk_size": 10000,
            "bootstrap_iteration": 100,
            "bootstrap_factor": 0.5,
            "rng_seed": 12345,
        },
        "precomputed_stats": {"path": str(precompute_combined_path)},
        "query_markers": {"serialized_lookup": str(cfg.query_marker_path)},
    }

    print("Running mapping...")
    FromSpecifiedMarkersRunner(args=[], input_data=mapping_config).run()


# ----------------------------
# OPTIONAL HELPERS (kept separate)
# ----------------------------

def convert_counts_csv_to_h5ad(csv_path: str, out_h5ad_path: str) -> None:
    """
    Convert a counts matrix CSV -> AnnData .h5ad.
    Assumes genes are rows and cell IDs are columns (adjust if needed).
    """
    import anndata
    import scipy.sparse

    counts_df = pd.read_csv(csv_path, index_col=0)
    counts_df_T = counts_df.T  # cells x genes
    counts_sparse = scipy.sparse.csr_matrix(counts_df_T.values)

    adata = anndata.AnnData(X=counts_sparse)
    adata.obs_names = counts_df_T.index.tolist()
    adata.var_names = counts_df_T.columns.tolist()

    adata.write(out_h5ad_path)
    print(f"Saved AnnData to: {out_h5ad_path}")


def map_mouse_symbols_to_ensembl(
    in_csv: str, out_csv: str, species: str = "mouse"
) -> None:
    """
    Map gene symbols (index) to Ensembl IDs using mygene.
    Expects gene symbols in the row index.
    """
    from mygene import MyGeneInfo

    df = pd.read_csv(in_csv, index_col=0)
    mg = MyGeneInfo()

    gene_symbols = df.index.tolist()
    query_results = mg.querymany(
        gene_symbols, scopes="symbol", fields="ensembl.gene", species=species
    )

    mapping = {}
    for entry in query_results:
        if "ensembl" in entry:
            if isinstance(entry["ensembl"], list):
                ensembl_id = entry["ensembl"][0].get("gene")
            else:
                ensembl_id = entry["ensembl"].get("gene")
            mapping[entry["query"]] = ensembl_id
        else:
            mapping[entry.get("query")] = None

    df["ensembl_id"] = df.index.to_series().map(mapping)
    df = df[~df["ensembl_id"].isnull()].copy()
    df.index = df["ensembl_id"]
    df = df.drop(columns=["ensembl_id"])

    df.to_csv(out_csv)
    print(f"Saved Ensembl-mapped CSV to: {out_csv}")


# ----------------------------
# MAIN
# ----------------------------

def main() -> None:
    cfg = Config()

    set_thread_env_vars()
    ensure_dirs(cfg.scratch_dir, cfg.precompute_dir, cfg.reference_dir, cfg.query_dir)

    print(f"abc_data_dir: {cfg.abc_data_dir.resolve()}")
    print(f"scratch_dir:  {cfg.scratch_dir.resolve()}")

    # ABC cache driven workflow (recommended)
    abc_cache = load_abc_cache(cfg)
    training_cells = select_training_cells(abc_cache)
    h5ad_path_list = build_h5ad_path_list(abc_cache, training_cells)

    run_precompute_reference_query_markers(cfg, abc_cache, training_cells, h5ad_path_list)
    run_mapping(cfg)

    print("Done.")


if __name__ == "__main__":
    main()
