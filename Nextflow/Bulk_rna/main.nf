#!/usr/bin/env nextflow
// Bulk RNA-seq pipeline: FastQC → STAR → sort/index → featureCounts.
nextflow.enable.dsl=2

//
// PARAMETERS (edit here or override on CLI)
//
params.reads   = 'path'
params.index   = 'path'
params.gtf     = 'path'
params.outdir  = './results'

//
// FASTQC
//
process fastqc {
  tag     { sample_id }
  cpus    2
  memory  '4 GB'
  publishDir "${params.outdir}/fastqc", mode: 'copy'

  input:
    // fromFilePairs will produce a tuple [ sample_id, [r1, r2] ]
    tuple val(sample_id), path(read1), path(read2)

  output:
    tuple val(sample_id), path("*_fastqc.zip"), path("*_fastqc.html")

  script:
  """
  fastqc -t ${task.cpus} -o ./ ${read1} ${read2}
  """
}

//
// STAR ALIGNMENT — produce an unsorted BAM
//
process star_align {
  tag     { sample_id }
  cpus    8
  memory  '32 GB'
  publishDir "${params.outdir}/bam_unsorted", mode: 'copy'

  input:
    tuple val(sample_id), path(read1), path(read2)

  output:
    tuple val(sample_id), path("${sample_id}.Aligned.out.bam")

  script:
  """
  STAR \
    --genomeDir ${params.index} \
    --readFilesIn ${read1} ${read2} \
    --readFilesCommand zcat \
    --runThreadN ${task.cpus} \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix ${sample_id}.
  """
}

//
// SORT & INDEX BAM
//
process sort_index {
  tag     { sample_id }
  cpus    2
  memory  '8 GB'
  publishDir "${params.outdir}/bam_sorted", mode: 'copy'

  input:
    tuple val(sample_id), path(bam)

  output:
    tuple val(sample_id),
          path("${sample_id}.sorted.bam"),
          path("${sample_id}.sorted.bam.bai")

  script:
  """
  samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam ${bam}
  samtools index ${sample_id}.sorted.bam
  """
}

//
// FEATURECOUNTS (exon counts, all BAMs in one go)
//
process featurecounts {
  tag     'featureCounts'
  cpus    4
  memory  '16 GB'
  publishDir "${params.outdir}/counts", mode: 'copy'

  input:
    // a single list of all sorted BAM paths
    val bam_list

  output:
    path("featurecounts.txt")

  script:
    // join the list into a space-separated string
    def bams = bam_list.join(' ')
  """
  featureCounts \
    -T ${task.cpus} \
    -a ${params.gtf} \
    -t exon \
    -g gene_id \
    -p \
    -o featurecounts.txt \
    ${bams}
  """
}

//
// WORKFLOW definition
//
workflow {

  // 1) discover and pair your reads
  //    each tuple: [ sample_id, r1_path, r2_path ]
  read_pairs = Channel
    .fromFilePairs( params.reads, flat: true )

  // 2) FastQC on both mates
  fastqc( read_pairs )

  // 3) STAR — 4) sort & index
  unsorted = star_align( read_pairs )
  sorted   = sort_index( unsorted )

  // 5) collect all sorted BAM paths into one list
  bam_list = sorted
    .map { sid, bam, bai -> bam }   // extract the BAM path
    .collect()

  // 6) run featureCounts once on that list
  featurecounts( bam_list )
}
