#!/bin/bash
#$ -N combined_patch
#$ -o path/analyze.log
#$ -j y
#$ -pe omp 28
#$ -l mem_per_core=8G
#$ -l h_rt=24:00:00
#$ -P selneuro

module load R/4.4.0      


cogent rna demux \
  -f path \
  -p path \
  -t smartseq_fla_umi \
  -o path \
  -b path \
  --fastqc



echo "Starting CogentAP analysis"
mkdir -p "path"

cogent rna analyze \
    -g mm39 \
    -i "path" \
    -o "path" \
    -G "path" \
    -t smartseq_fla_umi

echo "Finished"





