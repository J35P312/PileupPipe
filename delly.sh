#!/bin/bash -l
#SBATCH -A b2014152
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH -J pileup

module load bioinfo-tools
module load samtools
module load bcftools
module load vt
module load vep/87

./nextflow pileup_pipeline.nf --bam $1 --kg 1kg.tab --exac exac.tab --swefreq swefreq.tab --working_dir $2 --ref human_g1k_v37.fasta --genelist $3
