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

VEP_exec_file="variant_effect_predictor.pl"
frequency_script="/home/jesperei/new_pipeline/exac_annotation_sqlite.py"
print_variants="/home/jesperei/new_pipeline/print_variant.py"
excel_script="/home/jesperei/new_pipeline/CCCTG.py"
CADD_script="/home/jesperei/new_pipeline/annotate_vcf_cadd.py"
vep_cache="/home/jesperei/.vep"
clinvar_script="/home/jesperei/new_pipeline/ClinVar_annotate.py"

if (params.bam != ""){
    if(!file(params.bam).exists()) exit 1, "Missing bam"
    bam_file=file(params.bam)

    process mpileup{
        publishDir params.working_dir , mode: 'copy', overwrite: true

        input:
        file bam_file
    
        output:
        file ("${bam_file.baseName}.snps.vcf") into SNP_vcf

        """

        samtools mpileup -uf ${params.ref} ${params.bam} | bcftools call -cv > ${bam_file.baseName}.snps.raw.vcf
        vt decompose ${bam_file.baseName}.snps.raw.vcf -o ${bam_file.baseName}.snps.decomposed.vcf
        vt normalize ${bam_file.baseName}.snps.decomposed.vcf -r ${params.ref} -o ${bam_file.baseName}.snps.vcf
        if [ "" != ${params.dbSNP} ]
        then
            bgzip ${bam_file.baseName}.snps.vcf
            bcftools index ${bam_file.baseName}.snps.vcf.gz
            bcftools  annotate -a ${params.dbSNP} -c ID ${bam_file.baseName}.snps.vcf.gz > ${bam_file.baseName}.snps.vcf
        fi
        """
    }   

}else{
    if(!file(params.vcf).exists()) exit 1, "Missing bam or vcf file"
    SNP_vcf=file(params.vcf)

}


process annotate{
    publishDir params.working_dir , mode: 'copy', overwrite: true

    input:
    file SNP_vcf
    
    output:
    file ("${SNP_vcf.baseName}.vep.vcf") into annotated_SNP_vcf
    file ("${SNP_vcf.baseName}.vep.filt.vcf") into annotated_filtered_SNP_vcf
    file ("${SNP_vcf.baseName}.vep.filt.xls") into annotated_filtered_SNP_xls

    """
    ${VEP_exec_file} -i ${SNP_vcf} -o ${SNP_vcf.baseName}.vep.vcf --force_overwrite --hgvs --fasta ${params.ref}  --sift b --polyphen b  --vcf --offline --per_gene --no_intergenic --cache ${vep_cache}
    grep -E "MODERATE|HIGH|#" ${SNP_vcf.baseName}.vep.vcf > ${SNP_vcf.baseName}.vep.filt.vcf

    if [ "" != ${params.genelist} ]
    then
        python ${print_variants} ${SNP_vcf.baseName}.vep.filt.vcf > tmp.vcf
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

    python ${excel_script} --vcf ${SNP_vcf.baseName}.vep.filt.vcf
    """
}

