# Roussarie Lab Work Code

This repository is a grab-bag of day-to-day scripts, pipelines, and helper utilities used in the lab. The code is organized primarily by language/tooling (Bash, Python, R, Snakemake, Nextflow, and TE Tools). Many scripts are templates with **placeholder paths** (e.g., `path`) and **sample IDs** that should be updated before use.

## Repository layout

- `Bash/` — cluster job scripts (SGE) and command-line helpers for sequencing workflows.
- `Nextflow/` — Nextflow pipelines for bulk RNA-seq processing.
- `Python/` — analysis notebooks converted to scripts and utility pipelines.
- `R/` — Shiny apps, DESeq2 workflows, and Seurat helpers.
- `Snakemake/` — Snakemake pipelines for trimming, alignment, and quantification.
- `TETOOLS/` — RepeatMasker / Dfam TE Tools wrappers.

Each language folder contains its own `README.md` with more details.

## Conventions used in the scripts

- **Placeholder paths**: replace `path` with your local or cluster-specific paths.
- **Sample IDs**: replace `sample1`, `sample2`, etc. with your actual sample names.
- **Schedulers**: several Bash scripts assume Sun Grid Engine (SGE) headers (`#$` lines). Adjust for SLURM or other schedulers as needed.
- **Environment modules**: scripts that call `module load` assume a module environment is available.

## Quick start

1. Pick the folder that matches the workflow/tool you need.
2. Open the README in that folder for a brief workflow summary.
3. Update paths, sample IDs, and any environment-specific values.
4. Run the script or pipeline in your preferred environment.

## Notes

These scripts are intentionally lightweight and meant to be adapted. If you find a workflow you want to standardize, consider parameterizing it (e.g., with config files) for reproducibility.
