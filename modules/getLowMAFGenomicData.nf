process getLowMAFGenomicData {
    storeDir 'outputs/cohort/final_cohort_data/genomic_data/'

    input:
        tuple(
            val(chr_num),
            path('ancestrally_matched_cases_and_controls_platekeys.txt'),
            path('ancestrally_matched_cases_platekeys.txt')
        )

    output:
        path "final_cohort_genomic_data_to_be_scored_chr${chr_num}.*", emit: bfiles

    script:
        """      
        module load plink/1.9

        BFILE_NAME=\$(ls ${params.source_lowMAF_bfiles_dir} | grep .\\*${params.bfile_chr_num_prefix}${chr_num}${params.bfile_chr_num_suffix}.\\*.bim)
        BFILE_NAME=\${BFILE_NAME::-4} # removes the .bim from the filename

        plink -bfile ${params.source_lowMAF_bfiles_dir}/\$BFILE_NAME \
            -keep ancestrally_matched_cases_and_controls_platekeys.txt \
            --allow-no-sex \
            --make-pheno ancestrally_matched_cases_platekeys.txt '*' \
            --make-bed \
            --out final_cohort_genomic_data_to_be_scored_chr${chr_num}
        """
}
