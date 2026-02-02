# Nextflow workflows

This folder contains Nextflow DSL2 pipelines, currently focused on bulk RNA-seq processing.

## Contents

- `Bulk_rna/` — end-to-end bulk RNA-seq (FastQC → STAR → sort/index → featureCounts).

## Usage notes

- Update `params.reads`, `params.index`, and `params.gtf` in the `.nf` files or pass them on the CLI.
- Output is written under `params.outdir` (default `./results`).
