process removeRelatedWithinCohorts {
    storeDir 'outputs/cohort/non_final_cohort_data'

    input:
        path(cohort_platekeys)

    output:
        path "un${cohort_platekeys}"

    script:
        """
        module load king/2.3.2
        module load plink/1.9

        which_cohort=\$(grep -Po "(?<=related_).*(?=_platekeys.txt)" <<< ${cohort_platekeys})

        plink -bfile ${params.source_ancestry_bfiles} \
        --keep-fam ${cohort_platekeys} \
        --allow-no-sex \
        --make-bed --out \${which_cohort}_data

        king -b \${which_cohort}_data.bed --related --degree 2

        plink --bfile \${which_cohort}_data \
        --remove king.kin0 \
        --allow-no-sex \
        --make-bed --out unrelated_\${which_cohort}_data

        awk '{print \$1, \$1}' unrelated_\${which_cohort}_data.fam > unrelated_\${which_cohort}_platekeys.txt
        """
}

process removeRelatedBetweenCohorts {
    storeDir 'outputs/cohort/non_final_cohort_data'

    input:
        path internally_unrelated_cohorts // This is a list

    output:
        tuple(
            path('completely_unrelated_cohort.bim'),
            path('completely_unrelated_cohort.bed'),
            path('completely_unrelated_cohort.fam')
        )
        path('completely_unrelated_cohort_platekeys.txt')

    script:
        """
        module load king/2.3.2
        module load plink/1.9
        module load python/3.11

        cat ${internally_unrelated_cohorts} > internally_unrelated_whole_cohort_platekeys.txt

        plink --bfile ${params.source_ancestry_bfiles} \
        --keep-fam internally_unrelated_whole_cohort_platekeys.txt \
        --allow-no-sex \
        --make-bed \
        --out internally_unrelated_whole_cohort

        king -b internally_unrelated_whole_cohort.bed --related --degree 2

        remove-related-controls.py ${internally_unrelated_cohorts} king.kin0 to_keep.txt

        plink --bfile internally_unrelated_whole_cohort \
        --keep to_keep.txt \
        --allow-no-sex \
        --make-pheno unrelated_case_platekeys.txt '*' \
        --make-bed \
        --out completely_unrelated_cohort

        awk '{print \$1, \$1}' completely_unrelated_cohort.fam > completely_unrelated_cohort_platekeys.txt
        """
}

