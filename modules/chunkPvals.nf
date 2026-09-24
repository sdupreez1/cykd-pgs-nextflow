process chunkPvals {
    publishDir { "outputs/comparison_results/pvals/${platform}/" }, mode: 'copy'

    input:
        tuple(
            val(chunk_id),
            val(platform),
            path(chunk_traits),
            path(reference_panel_platekeys),
            path(case_and_control_platekeys),
            path(case_platekeys),
            path(final_cohort_platekeys_sex_age)
        )

    output:
        tuple val("${platform}"), path("${chunk_id}_pvals.tsv")

    script:
        """
        module load R/4.5.2
        compare_case_and_control_pgs_scores.R ${platform} ${chunk_id} ${projectDir}
        """
}
