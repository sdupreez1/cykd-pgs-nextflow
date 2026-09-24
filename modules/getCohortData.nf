process getCohortData {
    publishDir 'outputs/cohort/non_final_cohort_data'
    
    input:
        val  which_cohort
        path filter_data_files 

    output:
        path "related_${which_cohort}_platekeys.txt"        , emit:     related_cohort_platekeys
        path "related_${which_cohort}_participant_ids.txt"  , optional: true
        val  "${which_cohort}"                              , emit:     which_cohort

    script:
        """
        module load R/4.5.2
        get_cohort_data.R ${which_cohort} ${params.provide_participant_ids} 
        """
}
