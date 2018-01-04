#!/bin/bash -l
#SBATCH -A b2014152
#SBATCH -p core
#SBATCH -n 4
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
#argument3: the config file

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

./nextflow pileup_pipeline.nf --bam $TMPDIR/$filename --working_dir $TMPDIR/$sample -w $TMPDIR -c $3
cp -rf $TMPDIR/$sample/ $2/
