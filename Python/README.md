# Python analyses

This folder collects Python scripts (often exported from notebooks) used in the lab for various analyses.

## Subfolders

- `Resonance/` — resonance frequency analysis for patch-seq traces.
- `Scrublet/` — doublet detection workflow using Scrublet.
- `cell_type_mapper/` — MapMyCells / cell_type_mapper workflow helpers.
- `patch_seq/` — scVI/scANVI integration helpers.

## Usage notes

- Many scripts use placeholder paths; update these before running.
- Some scripts assume notebook-style commands (`pip install`, `!apt-get`). These are now commented so the file can be imported as a script; run installs in your environment before execution.
