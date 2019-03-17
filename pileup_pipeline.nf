//path of the vep executable file

print_variants="${params.pileup_pipeline_home}/internal_scripts/print_variant.py"
excel_script="${params.pileup_pipeline_home}/internal_scripts/CCCTG.py"
vep_cache="${params.pileup_pipeline_home}/.vep"
    
if(!file(params.bam).exists()) exit 1, "Missing bam:, use either --bam to provide the path to an indexed bam file"
bam_file=file(params.bam)



process gatk{
        publishDir params.working_dir , mode: 'copy', overwrite: true

        input:
		file bam_file

        output:
        file ("${bam_file.baseName}.gatk.vcf") into SNP_vcf
        """

        singularity exec ${params.pileup_pipeline_home}/PileUp.simg gatk HaplotypeCaller -R ${params.ref} -I ${params.bam} -O ${bam_file.baseName}.gatk.vcf
        singularity exec ${params.pileup_pipeline_home}/PileUp.simg vt decompose  ${bam_file.baseName}.gatk.vcf -o ${bam_file.baseName}.gatk.decomposed.vcf
        singularity exec ${params.pileup_pipeline_home}/PileUp.simg vt normalize ${bam_file.baseName}.gatk.decomposed.vcf -r ${params.ref} -o ${bam_file.baseName}.gatk.vcf

        """
}

process annotate{

    publishDir params.working_dir , mode: 'copy', overwrite: true

    input:
    file SNP_vcf
    
    output:
    file ("${SNP_vcf.baseName}.vep.filt.vcf") into annotated_filtered_SNP_vcf
    file ("${SNP_vcf.baseName}.vep.filt.xls") into annotated_filtered_SNP_xls

    """
    ${params.VEP_exec_file} -i ${SNP_vcf} -o ${SNP_vcf.baseName}.vep.vcf ${params.vep_command} --dir ${vep_cache} --fasta ${params.ref}
    mv ${SNP_vcf.baseName}.vep.vcf ${SNP_vcf.baseName}.vep.filt.vcf

    if [ "" != ${params.genelist} ]
    then
        python ${print_variants} ${SNP_vcf.baseName}.vep.filt.vcf ${params.genelist} > tmp.vcf
        mv tmp.vcf ${SNP_vcf.baseName}.vep.filt.vcf 
    fi

    python ${excel_script} --vcf ${SNP_vcf.baseName}.vep.filt.vcf --frequency 0.025
    singularity exec ${params.pileup_pipeline_home}/PileUp.simg bgzip ${SNP_vcf.baseName}.vep.filt.vcf
    singularity exec ${params.pileup_pipeline_home}/PileUp.simg tabix -p vcf ${SNP_vcf.baseName}.vep.filt.vcf 
    """
}
