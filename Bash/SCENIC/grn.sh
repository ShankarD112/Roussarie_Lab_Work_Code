#!/bin/bash
#$ -o /path/pyscenic_grn.out
#$ -j y
#$ -pe omp 16
#$ -l h_rt=48:00:00
#$ -P selneuro
#$ -N pyscenic_grn

echo "JOB STARTED on $(hostname) at $(date)"

# Run pySCENIC GRN step inside Singularity.
# Replace `path` placeholders with your Singularity image + data directories.
# Ensure `--num_workers` matches allocated cores.
#
# Path to Singularity image
SIF=/path/Singularity/aertslab-pyscenic-scanpy-0.12.1-1.9.1.sif

# Working directory with input/output files
WORKDIR=/path
cd $WORKDIR

# Run pySCENIC inside Singularity
singularity exec \
  --bind /path \
  --bind $PWD:$PWD \
  $SIF \
  pyscenic grn \
    sample_scenic_input.loom \
    feather_cistarget/allTFs_mm.txt \
    -o sample_adj.csv \
    --num_workers 16

echo "JOB FINISHED at $(date)"
