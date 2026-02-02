#!/bin/bash
#$ -o bam_subset
#$ -j y
#$ -pe omp 16
#$ -l h_rt=24:00:00
#$ -P selneuro

INPUT_BAM_DIR="path"
BARCODES_FILE="path"
OUTPUT_BAM_DIR="path"

samples=(
  sample1 sample2 sample3 sample4 sample5 sample6
  sample7 sample8 sample9 sample10 sample11 sample12
)

for sample in "${samples[@]}"; do
    echo "Processing sample ${sample}"
    
    # Path to Cell Ranger BAM
    input_bam="${INPUT_BAM_DIR}/${sample}/outs/possorted_genome_bam.bam"
    if [[ ! -f "$input_bam" ]]; then
        echo "ERROR: BAM not found at $input_bam, skipping." >&2
        continue
    fi

    # Per-sample output directory + filenames
    outdir="${OUTPUT_BAM_DIR}/${sample}"
    mkdir -p "$outdir"
    output_bam="${outdir}/${sample}_subset.bam"

    # Subset reads by barcode
    subset-bam \
        --bam            "$input_bam" \
        --cell-barcodes  "$BARCODES_FILE" \
        --out-bam        "$output_bam" \
        --cores          "$NSLOTS"

    # Index the new BAM
    samtools index -@ "$NSLOTS" "$output_bam"

    echo "Finished ${sample}"
done

echo "All samples completed."
