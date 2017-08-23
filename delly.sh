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



filename=$(basename $1)
extension="${filename##*.}"
sample="${filename%.*}"

mkdir $2
mkdir $2/$sample
./nextflow pileup_pipeline.nf --bam $1 --kg /proj/b2014152/private/exac/1kg_filter.tab --exac /proj/b2014152/private/exac/af_filter_data.tsv  --swefreq /proj/b2014152/private/exome_pipeline//swegen.tab --working_dir $2/$sample --ref /proj/b2016296/private/nobackup/annotation/human_g1k_v37.fasta --clinvar /proj/b2016296/private/nobackup/annotation/clinvar_20170801.vcf --dbSNP /proj/b2016296/private/nobackup/annotation/DBsnp_All_20170710.vcf.gz  --cadd_indels /proj/b2016296/private/nobackup/annotation/InDels.tsv.gz --cadd_snps  /proj/b2016296/private/nobackup/annotation/whole_genome_SNVs.tsv.gz --genelist $3 -with-trace

