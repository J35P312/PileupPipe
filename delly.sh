#!/bin/bash -l
#SBATCH -A b2014152
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 96:00:00
#SBATCH -J pileup

#these are the required variable, the path of these variables must be available
module load bioinfo-tools
module load samtools
module load bcftools
module load vt
module load vep/87
module load tabix
module load GATK/3.8-0
module load freebayes
module load picard/2.10.3

#argument1: a bam file
#argument2: an output folder
#argument3: a gene list

#exac, thousand genomes and swefreq frecuencies, these files are tab separated files on the format "chr pos alt frequency"
#cannot be skipped, but can be replaced by an empty file
exac="/proj/b2014152/private/exac/af_filter_data.tsv"
kg="/proj/b2014152/private/exac/1kg_filter.tab"
swefreq="/proj/b2014152/private/exome_pipeline//swegen.tab"
gnomad="/proj/b2014152/private/exome_pipeline/gnomad.tab"
#clinvar vcf
clinvar="/proj/b2016296/private/nobackup/annotation/clinvar_20170801.vcf"
#dbsnp vcf
dbSNP="/proj/b2016296/private/nobackup/annotation/DBsnp_All_20170710.vcf.gz"
#cadd indel and snp files, download via the cadd homepage
cadd_indels="/proj/b2016296/private/nobackup/annotation/InDels.tsv.gz"
cadd_snps="/proj/b2016296/private/nobackup/annotation/whole_genome_SNVs.tsv.gz"
#reference fasta
ref="/proj/b2016296/private/nobackup/annotation/human_g1k_v37.fasta"
#gatk
gatk="/sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar"
#picard
picard="/sw/apps/bioinfo/picard/2.10.3/milou/picard.jar"

filename=$(basename $1)
extension="${filename##*.}"
sample="${filename%.*}"

mkdir $2
echo "$sample"
mkdir $TMPDIR/$sample

cp $1 $TMPDIR/$filename
cp $1.bai $TMPDIR/$filename.bai


cp $ref $TMPDIR/human_g1k_v37.fasta
samtools faidx $TMPDIR/human_g1k_v37.fasta

java -jar $picard CreateSequenceDictionary REFERENCE=$TMPDIR/human_g1k_v37.fasta

./nextflow pileup_pipeline.nf --bam $TMPDIR/$filename --kg $kg --exac $exac  --swefreq $swefreq --working_dir $TMPDIR/$sample --ref $TMPDIR/human_g1k_v37.fasta --clinvar $clinvar  --dbSNP $dbsnp  --cadd_indels $cadd_indels --cadd_snps $cadd_snps --genelist $3 -w $TMPDIR --gatk $gatk --gnomad $gnomad
cp -rf $TMPDIR/$sample/ $2/
