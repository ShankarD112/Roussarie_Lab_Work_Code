# Output directory and merged BAM file
outdir="path"
outfile="${outdir}/merged_possorted.bam"

# Sample IDs
samples=(
  sample1 sample2 sample3 sample4

)

bam_list=$(mktemp)

for s in "${samples[@]}"; do
    echo "path/${s}/outs/possorted_genome_bam.bam" >> "$bam_list"
done

mkdir -p "$outdir"

samtools merge -@ "$threads" -O BAM "$outfile" -b "$bam_list"
samtools index -@ "$threads" "$outfile"

rm "$bam_list"

echo "Merged BAM written to: $outfile"
