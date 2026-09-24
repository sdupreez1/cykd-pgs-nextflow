process finalCohortStats {
    storeDir 'outputs/cohort/final_cohort'

    input:
        tuple(
            path('ancestrally_matched_cases_and_controls_platekeys.txt'),
            path('ancestrally_matched_cases_platekeys.txt')
        )

    output:
        path 'full_final_cohort_platekey_sex_age.tsv', emit: platekey_sex_age
        path 'final_cohort_summary_stats.tsv'

    script:
        """
        module load R/4.5.2
        calculate_cohort_stats.R
        """
}
