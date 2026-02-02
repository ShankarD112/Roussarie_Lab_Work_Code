# Snakemake pipelines

This folder contains Snakemake workflows for trimming, alignment, and counting.

## Files

- `complete.snake` — full pipeline (trim → FastQC → STAR → flagstat → MultiQC → Verse → combine).
- `samtools_sort_index.snake` — minimal pipeline to sort and index BAM files.

## Usage notes

- Replace placeholder paths and sample lists before running.
- Conda environment YAMLs referenced in the rules should be created/updated locally.
