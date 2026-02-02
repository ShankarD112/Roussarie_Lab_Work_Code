#!/bin/bash
#$ -o path/automated_repeatmasker.txt
#$ -j y
#$ -pe omp 16
#$ -l mem_per_core=16G	
#$ -l h_rt=48:00:00
#$ -P selneuro

cd path

for dir in path/*; do
  echo "Processing directory: $dir"
  for fafile in "$dir"/*_1pct.fa; do
    if [ -f "$fafile" ]; then
      # Get base filename without extension and path
      base=$(basename "$fafile" .fa)
      outdir="path/${base}_out"
      
      echo "Running dfam-tetools.sh on $fafile, output to $outdir"
      
      ./dfam-tetools.sh --library path -- \
        RepeatMasker -species human -qq -pa 16 -html -source "$fafile" -dir "$outdir"
    fi
  done
done
