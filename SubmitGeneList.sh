#!/bin/bash -l
#SBATCH -A b2014152
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
module load tabix
module load freebayes
module load GATK

#argument1: a bam file
#argument2: an output folder
#argument3: a gene list
#argument4: a config file

filename=$(basename $1)
extension="${filename##*.}"
sample="${filename%.*}"

mkdir $2
echo "$sample"
mkdir $TMPDIR/$sample

cp $1 $TMPDIR/$filename

samtools index $TMPDIR/$filename
#./nextflow pileup_pipeline.nf --bam $TMPDIR/$filename --working_dir $TMPDIR/$sample --genelist $3 -w $TMPDIR -c $4
./nextflow pileup_pipeline.nf --bam $TMPDIR/$filename --working_dir $TMPDIR/$sample --genelist $3 -c $4
cp -rf $TMPDIR/$sample/ $2/
