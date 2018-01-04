# PileupPipe

Run SNP variant calling using GATK, bcftools pileup, and freebayes. The pipeline performs annotation using vep and multiple optional databases.

# Run Locally
Analyse a single sample, using a gene list:

    ./nextflow pileup_pipeline.nf --bam <input_bam> --working_dir <output_folder> --genelist <genelist> -w $TMPDIR -c <config_file_from_setup>

Analyse a single sample, with no gene list:

    ./nextflow pileup_pipeline.nf --bam <input_bam> --working_dir <output_folder> -w $TMPDIR -c <config_file_from_setup>

annotate a vcf file

    ./nextflow pileup_pipeline.nf --vcf <input_vcf> --working_dir <output_folder> -w $TMPDIR -c <config_file_from_setup>


# Run on slurm

The slurm scripts are wrappers around the pileup_pipeline.nf script.

Analyse a single sample, using a gene list:

    sbatch SubmitGeneList.sh <input_bam> <output_folder> <genelist> <config_file_from_setup>

Analyse a single sample, with no gene list:

    sbatch SubmitNoGeneList.sh  <input_bam> <output_folder> <config_file_from_setup>

Analyse a folder containing bam files using a gene list:

    ./snpPipe.sh <input_folder_with_bam> <output_folder> <genelist>

# Install

Download and  install nextflow

    https://www.nextflow.io/

run the setup script

    python setup.py > config.conf


open the config file using a text editor and change the path variables:
    
    nano config.conf

The Picard tools (line 48), GATK (line 47), dbSNP(line 40) and reference (line 49) paths are required:

You need to install TIDDIT and put the TIDDIT directory inside the PileupPipe directory.

Additionally, you need to put a .vep folder containing the vep database insisde the PileupPipe folder.


Per default, VEP will run the merged refseq DB. The VEP command may be edited in the nextflow script (Pileup_pipeline.nf, line 163).

the rest of the paths are optional.


If you will use the slurm wrappers, you need to edit the slurm account in the files, SubmitAnnotation.sh, SubmitGeneList.sh, and SubmitNoGeneList.sh


# Gene list

The gene list is a text file. Each line of the text file contains the HGNC symbol of a gene.
