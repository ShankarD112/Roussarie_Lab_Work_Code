#!/bin/bash
# Build a Cell Ranger aggregation CSV and run `cellranger aggr`.
# Replace placeholder paths + sample IDs before submitting to the cluster.
#$ -N Sample_aggr_exon              # job name
#$ -o Sample_aggr_exon.out          # stdout log
#$ -e Sample_aggr_exon.err          # stderr log
#$ -pe omp 8                     # 8 CPU cores
#$ -l mem_per_core=8G
#$ -l h_rt=24:00:00              # wall-clock limit
#$ -P selneuro                   # SCC project

module load bcl2fastq
module load cellranger

# Base directory where all per-sample Cell Ranger runs live
OUT_BASE="path"

# List of samples that have finished `cellranger count`
samples=(
  sample1 sample2 sample3 sample4
)

echo "Building aggregation CSV..."

CSV="${OUT_BASE}/aggr_libraries.csv"
# IMPORTANT: use sample_id for Cell Ranger 6+
echo "sample_id,molecule_h5" > "$CSV"

for sample in "${samples[@]}"; do
  H5="${OUT_BASE}/${sample}/outs/molecule_info.h5"

  if [[ ! -f "$H5" ]]; then
    echo "ERROR: Missing molecule_info.h5 for sample ${sample}: ${H5}" >&2
    exit 1
  fi

  echo "${sample},${H5}" >> "$CSV"
done

echo "Aggregation CSV created at: $CSV"
echo "Running cellranger aggr..."

# Run Cell Ranger aggr from OUT_BASE so output is nicely contained there
cd "$OUT_BASE" || exit 1

cellranger aggr \
  --id=all_samples \
  --csv="$CSV" \
  --normalize=none \
  --localmem=64 \
  --localcores=8 \
  --output-dir="${OUT_BASE}/all_samples"

echo "Aggregation complete."
echo "Combined dataset is in: ${OUT_BASE}/all_samples/outs"
