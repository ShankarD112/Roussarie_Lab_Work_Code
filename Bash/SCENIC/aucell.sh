#!/bin/bash
#$ -o /path/pyscenic_aucell.out
#$ -j y
#$ -pe omp 16
#$ -l h_rt=48:00:00
#$ -P selneuro
#$ -N pyscenic_aucell

echo "JOB STARTED on $(hostname) at $(date)"

# Run pySCENIC AUCell scoring in Singularity.
# Replace `path` placeholders with your working directory and loom outputs.

SIF=/path/Singularity/aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif
WORKDIR=/path/

cd "$WORKDIR"

singularity exec \
  --bind /path \
  --bind "$PWD":"$PWD" \
  "$SIF" \
  pyscenic aucell \
    sample_scenic_input.loom \
    sample_reg.csv \
    --output sample_SCENIC_output.loom \
    --num_workers 16

echo "JOB FINISHED at $(date)"
