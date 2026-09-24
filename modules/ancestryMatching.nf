process getGenomicPrincipalComponenets {
    // Only using storeDir when all that is emitted from a process is paths that are stored in the same storeDir 
    // (prevents vals or publishDir paths from being emitted when process is skipped)
    storeDir 'outputs/cohort/genomic_PCs'

    input:
        tuple(
            path('completely_unrelated_cohort.bim'),
            path('completely_unrelated_cohort.bed'),
            path('completely_unrelated_cohort.fam')
        )
        path 'completely_unrelated_cohort_platekeys.txt'

    output:
        path 'plink.eigenvec'
        path 'plink.eigenval'
        path 'completely_unrelated_cohort.fam'

    script:
        """
        module load plink/1.9
        plink --bfile completely_unrelated_cohort --pca 10 --allow-no-sex
        """
}

process ancestryMatching {
    storeDir 'outputs/cohort/final_cohort_data'

    input:
        path 'plink.eigenvec'
        path 'plink.eigenval'
        path 'completely_unrelated_cohort.fam'

    output:
        path 'PCA.pdf'
        path 'plot_ancestry_matched_t0.002_x2.pdf'
        tuple(
            path('ancestrally_matched_cases_and_controls_platekeys.txt'),
            path('ancestrally_matched_cases_platekeys.txt'), 
            emit: platekeys
        )

    script:
        """
        module load R/4.5.2
        pca_case_control_matching.R plink.eigenvec plink.eigenval completely_unrelated_cohort.fam
        """
}
