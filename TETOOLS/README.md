# TE Tools wrappers

This folder contains helpers to run Dfam TE Tools / RepeatMasker.

## Files

- `dfam-tetools.sh` — wrapper for running the Dfam TE Tools container via Docker or Singularity.
- `automated_repeatmasker.sh` — batch driver that calls `dfam-tetools.sh` across directories.

## Usage notes

- Update the placeholder `path` entries to match your data locations.
- The scripts assume SGE directives and a cluster environment for long-running jobs.
