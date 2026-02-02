#!/bin/bash
#$ -N Name 
#$ -o results
#$ -j y
#$ -pe omp 16
#$ -l mem_per_core=8G
#$ -l h_rt=24:00:00
#$ -P selneuro

nextflow run main.nf   -with-report report.html   -with-trace trace.txt   -with-timeline timeline.html
