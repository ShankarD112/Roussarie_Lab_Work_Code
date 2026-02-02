#!/bin/bash
#$ -o bcl_convert
#$ -j y
#$ -pe omp 16
#$ -l h_rt=24:00:00
#$ -P selneuro

module load bcl-convert

# Convert Illumina BCLs to FASTQs.
# Replace /home/path placeholders with your run folder + output folder.
bcl-convert \
  --bcl-input-directory /home/path/ \
  --output-directory /home/path/ \
  --no-lane-splitting true \
  --no-sample-sheet true \
  --force
