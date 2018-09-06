//path of the vep executable file

tiddit="${params.pileup_pipeline_home}/TIDDIT/bin/TIDDIT"
print_variants="${params.pileup_pipeline_home}/internal_scripts/print_variant.py"
excel_script="${params.pileup_pipeline_home}/internal_scripts/CCCTG.py"
CADD_script="${params.pileup_pipeline_home}/internal_scripts/annotate_vcf_cadd.py"
vep_cache="${params.pileup_pipeline_home}/.vep"
clinvar_script="${params.pileup_pipeline_home}/internal_scripts/ClinVar_annotate.py"
sex_check="${params.pileup_pipeline_home}/internal_scripts/sex_check.py"

if (params.bam != ""){
    
    if(!file(params.bam).exists()) exit 1, "Missing bam"
    bam_file=file(params.bam)
    process tiddit{

        publishDir params.working_dir , mode: 'copy', overwrite: true

        input:
	file bam_file

        output:
        file ("${bam_file.baseName}.sex.tab") into sex_tab
        file ("${bam_file.baseName}.coverage.tab") into coverage_tab

        """
        ${tiddit} --cov -b ${params.bam} -z 10000 -o ${bam_file.baseName}.coverage
        python ${sex_check} ${bam_file.baseName}.coverage.tab ${params.ref}.fai > ${bam_file.baseName}.sex.tab

	"""
    }

    process mpileup{
        publishDir params.working_dir , mode: 'copy', overwrite: true

        input:
        file bam_file
        file sex_tab

        output:
        file ("${bam_file.baseName}.snps.samtools.vcf") into RAW_SNP_pileup_vcf
        """
        samtools mpileup -uf ${params.ref} ${params.bam} | bcftools call -cv --ploidy-file ${sex_tab} > ${bam_file.baseName}.snps.samtools.raw.vcf

        vt decompose ${bam_file.baseName}.snps.samtools.raw.vcf -o ${bam_file.baseName}.snps.samtools.decomposed.vcf 
        vt normalize ${bam_file.baseName}.snps.samtools.decomposed.vcf -r ${params.ref} -o ${bam_file.baseName}.snps.samtools.vcf
        """
    }

    process freebayes{
        publishDir params.working_dir , mode: 'copy', overwrite: true

        input:
		file bam_file

        output:
        file ("${bam_file.baseName}.snps.freebayes.vcf") into RAW_SNP_freebayes_vcf
        """
        freebayes -f ${params.ref} ${params.bam} | grep -v "	0/0:" > ${bam_file.baseName}.snps.freebayes.raw.vcf

        vt decompose  ${bam_file.baseName}.snps.freebayes.raw.vcf -o ${bam_file.baseName}.snps.freebayes.decomposed.vcf
        vt normalize ${bam_file.baseName}.snps.freebayes.decomposed.vcf -r ${params.ref} -o ${bam_file.baseName}.snps.freebayes.vcf
        """
    }

    process vt_and_merge{
        publishDir params.working_dir , mode: 'copy', overwrite: true

        input:
	file RAW_SNP_pileup_vcf
        file RAW_SNP_freebayes_vcf
        file bam_file

        output:
        file ("${bam_file.baseName}.vcf") into SNP_vcf


        """
        java -jar ${params.gatk} -T CombineVariants -R ${params.ref} --variant:samtools ${RAW_SNP_pileup_vcf.baseName}.vcf --variant:freebayes ${RAW_SNP_freebayes_vcf.baseName}.vcf -o ${bam_file.baseName}.merged.vcf -genotypeMergeOptions PRIORITIZE -priority freebayes,samtools

        vt decompose ${bam_file.baseName}.merged.vcf -o ${bam_file.baseName}.merged.decomposed.vcf
        vt normalize ${bam_file.baseName}.merged.decomposed.vcf -r ${params.ref} -o ${bam_file.baseName}.vcf
	"""
    }

}else{
    if(!file(params.vcf).exists()) exit 1, "Missing bam or vcf file"
    RAW_SNP_vcf=file(params.vcf)

    process vt_no_merge{
	publishDir params.working_dir , mode: 'copy', overwrite: true

        input:
	file RAW_SNP_vcf

        output:
        file ("${RAW_SNP_vcf.baseName}.vcf") into SNP_vcf


        """

        vt decompose ${RAW_SNP_vcf} -o ${RAW_SNP_vcf.baseName}.decomposed.vcf
        vt normalize ${RAW_SNP_vcf.baseName}.decomposed.vcf -r ${params.ref} -o ${RAW_SNP_vcf.baseName}.normalised.vcf
	
	rm ${RAW_SNP_vcf.baseName}.vcf
	mv ${RAW_SNP_vcf.baseName}.normalised.vcf ${RAW_SNP_vcf.baseName}.vcf
	
	"""
    }
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
    
    if [ "" != ${params.cadd_indels} ]
    then
        python ${CADD_script} --vcf ${SNP_vcf.baseName}.vep.filt.vcf --cadd ${params.cadd_indels} > tmp.vcf
        mv tmp.vcf ${SNP_vcf.baseName}.vep.filt.vcf 
    fi
    if [ "" != ${params.cadd_snps} ]
    then
        python ${CADD_script} --vcf ${SNP_vcf.baseName}.vep.filt.vcf --cadd ${params.cadd_snps} > tmp.vcf
        mv tmp.vcf ${SNP_vcf.baseName}.vep.filt.vcf 
    fi
    if [ "" != ${params.clinvar} ]
    then
        python ${clinvar_script} --vcf ${SNP_vcf.baseName}.vep.filt.vcf --db ${params.clinvar} > tmp.vcf
        mv tmp.vcf ${SNP_vcf.baseName}.vep.filt.vcf 
    fi

    python ${excel_script} --vcf ${SNP_vcf.baseName}.vep.filt.vcf --frequency 0.025
    """
}
