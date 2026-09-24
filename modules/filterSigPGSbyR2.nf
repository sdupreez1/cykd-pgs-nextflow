process process filterSigPGSbyR2 {
    publishDir 'outputs/comparison_results/', mode: 'copy'

    input:
        path(significant_pgs_list)

    output:
        path("combined_sig_pgs.tsv")
        path("final_significant_pgs.tsv")

    script:
        """
        # FNR - Current line number of file (FNR==1 means header of current file)
        # NR - Total number of lines seen so far (NR==1 means header of first file)
        awk 'FNR==1 && NR!=1 { next } { print }' ${significant_pgs_list} > combined_sig_pgs.tsv
        module load R/4.5.2
        filter_sig_pgs_by_r2.R ${projectDir}
        """
}
