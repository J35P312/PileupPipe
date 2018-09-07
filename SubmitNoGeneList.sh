#!/bin/bash -l
#SBATCH -A sens2017130
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 96:00:00
#SBATCH -J pileup

#these are the required variable, the path of these variables must be available
module load bioinfo-tools
module load Nextflow
module load samtools
module load bcftools
module load vt
module load vep
module load tabix
module load freebayes
module load Nextflow
module load GATK

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

samtools index $TMPDIR/$filename
#nextflow pileup_pipeline.nf --bam $TMPDIR/$filename --working_dir $TMPDIR/$sample -w $TMPDIR -c $3
nextflow pileup_pipeline.nf --bam $TMPDIR/$filename --working_dir $TMPDIR/$sample -c $3
p -rf $TMPDIR/$sample/ $2/
