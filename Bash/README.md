# Bash scripts

This folder contains shell scripts used for day-to-day sequencing workflows and cluster job submission (mostly SGE). Most scripts are templates with placeholder paths.

## Subfolders

- `BAMs/` — merge or subset BAM files.
- `Patch_seq/` — Patch-seq pipeline helpers (Cogent + FASTQ pooling).
- `SCENIC/` — pySCENIC workflow steps run inside Singularity.
- `bcl_convert/` — BCL conversion helper.
- `cellranger_scripts/` — Cell Ranger count + aggr helpers.

## Usage notes

- Replace `path` with real paths.
- Replace sample lists with your sample names.
- Many scripts assume SGE (`#$` directives) and `module load` is available.
