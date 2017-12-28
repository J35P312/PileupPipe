params.bam=""
params.vcf=""

params.ref=""
params.working_dir=""
params.exac='""'
params.kg='""'
params.swefreq='""'
params.genelist='""'
params.cadd_indels='""'
params.cadd_snps='""'
params.dbSNP='""'
params.clinvar='""'
params.gnomad='""'
params.gatk='""'


//path of the vep executable file
VEP_exec_file="variant_effect_predictor.pl"
//path of the tiddit executable file
tiddit="/proj/b2016296/private/nobackup/annotation/TIDDIT/bin/TIDDIT"

//This variable needs to be set to the path of the pipeline folder
pileup_pipeline_home="/home/jesperei/PileupPipe"

frequency_script="${pileup_pipeline_home}/exac_annotation_sqlite.py"
print_variants="${pileup_pipeline_home}/print_variant.py"
excel_script="${pileup_pipeline_home}/CCCTG.py"
CADD_script="${pileup_pipeline_home}/annotate_vcf_cadd.py"
vep_cache="${pileup_pipeline_home}/.vep"
clinvar_script="${pileup_pipeline_home}/ClinVar_annotate.py"
sex_check="${pileup_pipeline_home}/sex_check.py"

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
        errorStrategy 'retry'
        maxRetries 20

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

    process GATK{
        publishDir params.working_dir , mode: 'copy', overwrite: true

        input:
        file bam_file

        output:
        file ("${bam_file.baseName}.snps.GATK.vcf") into RAW_SNP_GATK_vcf
        """
        
        java -Xmx8G -jar ${params.gatk} -T HaplotypeCaller -R ${params.ref} -I ${params.bam} --genotyping_mode DISCOVERY -stand_call_conf 10 -o ${bam_file.baseName}.snps.GATK.raw.vcf  -mmq 10

        vt decompose ${bam_file.baseName}.snps.GATK.raw.vcf -o ${bam_file.baseName}.snps.GATK.decomposed.vcf
        vt normalize ${bam_file.baseName}.snps.GATK.decomposed.vcf -r ${params.ref} -o ${bam_file.baseName}.snps.GATK.vcf
        """
    }




    process vt_and_merge{
        publishDir params.working_dir , mode: 'copy', overwrite: true
        errorStrategy 'retry'
        maxRetries 20

        input:
		file RAW_SNP_pileup_vcf
        file RAW_SNP_freebayes_vcf
        file RAW_SNP_GATK_vcf
        file bam_file

        output:
        file ("${bam_file.baseName}.vcf") into SNP_vcf


        """
        java -jar ${params.gatk} -T CombineVariants -R ${params.ref} --variant:GATK ${RAW_SNP_GATK_vcf} --variant:samtools ${RAW_SNP_pileup_vcf.baseName}.vcf --variant:freebayes ${RAW_SNP_freebayes_vcf.baseName}.vcf -o ${bam_file.baseName}.merged.vcf -genotypeMergeOptions PRIORITIZE -priority freebayes,samtools,GATK

        vt decompose ${bam_file.baseName}.merged.vcf -o ${bam_file.baseName}.merged.decomposed.vcf
        vt normalize ${bam_file.baseName}.merged.decomposed.vcf -r ${params.ref} -o ${bam_file.baseName}.vcf

        if [ "" != ${params.dbSNP} ]
        then
            bgzip ${bam_file.baseName}.vcf
            bcftools index ${bam_file.baseName}.vcf.gz
            bcftools  annotate -a ${params.dbSNP} -c ID ${bam_file.baseName}.vcf.gz > ${bam_file.baseName}.vcf
        fi
	"""
    }

}else{
    if(!file(params.vcf).exists()) exit 1, "Missing bam or vcf file"
    RAW_SNP_vcf=file(params.vcf)

    process vt_no_merge{
	publishDir params.working_dir , mode: 'copy', overwrite: true
        errorStrategy 'retry'
        maxRetries 20

        input:
	file RAW_SNP_vcf

        output:
        file ("${RAW_SNP_vcf.baseName}.vcf") into SNP_vcf


        """

        vt decompose ${RAW_SNP_vcf} -o ${RAW_SNP_vcf.baseName}.decomposed.vcf
        vt normalize ${RAW_SNP_vcf.baseName}.decomposed.vcf -r $TMPDIR/ref.fa -o ${RAW_SNP_vcf.baseName}.vcf
        if [ "" != ${params.dbSNP} ]
        then
            bgzip ${RAW_SNP_vcf.baseName}.vcf
            bcftools index ${RAW_SNP_vcf.baseName}.vcf.gz
            bcftools  annotate -a ${params.dbSNP} -c ID ${RAW_SNP_vcf.baseName}.vcf.gz > ${RAW_SNP_vcf.baseName}.vcf
        fi
	"""
    }
}

process annotate{
    publishDir params.working_dir , mode: 'copy', overwrite: true
    errorStrategy 'retry'
    maxRetries 20


    input:
    file SNP_vcf
    
    output:
    file ("${SNP_vcf.baseName}.vep.filt.vcf") into annotated_filtered_SNP_vcf
    file ("${SNP_vcf.baseName}.vep.filt.xls") into annotated_filtered_SNP_xls

    """
    ${VEP_exec_file} -i ${SNP_vcf} -o ${SNP_vcf.baseName}.vep.vcf --force_overwrite --hgvs --symbol --fasta ${params.ref}  --sift b --polyphen b  --vcf --offline --per_gene --no_intergenic --cache --dir ${vep_cache} --symbol --check_existing
    mv ${SNP_vcf.baseName}.vep.vcf ${SNP_vcf.baseName}.vep.filt.vcf

    if [ "" != ${params.genelist} ]
    then
        python ${print_variants} ${SNP_vcf.baseName}.vep.filt.vcf ${params.genelist} > tmp.vcf
        mv tmp.vcf ${SNP_vcf.baseName}.vep.filt.vcf 
    fi
    
    if [ "" != ${params.exac} ]
    then
        python ${frequency_script} --vcf ${SNP_vcf.baseName}.vep.filt.vcf --db ${params.exac} --tag EXACAF > tmp.vcf
        mv tmp.vcf ${SNP_vcf.baseName}.vep.filt.vcf 
    fi

    if [ "" != ${params.kg} ]
    then
        python ${frequency_script} --vcf ${SNP_vcf.baseName}.vep.filt.vcf --db ${params.kg} --tag 1000GAF  > tmp.vcf
        mv tmp.vcf ${SNP_vcf.baseName}.vep.filt.vcf 
    fi

    if [ "" != ${params.swefreq} ]
    then
        python ${frequency_script} --vcf ${SNP_vcf.baseName}.vep.filt.vcf --db ${params.swefreq} --tag SWEAF   > tmp.vcf
        mv tmp.vcf ${SNP_vcf.baseName}.vep.filt.vcf 
    fi

    if [ "" != ${params.gnomad} ]
    then
        python ${frequency_script} --vcf ${SNP_vcf.baseName}.vep.filt.vcf --db ${params.gnomad} --tag GNOMADAF   > tmp.vcf
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
