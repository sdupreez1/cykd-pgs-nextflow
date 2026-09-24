process process calculatePGS {
    // No publish/storeDir since pgsc_calc itself can output to an accesible location (check --outidr), easier to handle chunk_id than passing a directory between processes
    input:
        tuple(
            val(chunk_id),
            val(platform),
            path(formatted_scorefiles, stageAs: "formatted_scorefiles/formatted_scorefile_?.txt"), 
            path(low_maf_data,         stageAs: "cohort_genomic_data/*")
        )

    output:
        val "${chunk_id}"

    script:
        """
        module load singularity/4.1.1
        module load java

        echo 'sampleset,path_prefix,chrom,format' > samplesheet.csv
        for chr_num in \$(seq 1 22) X; do
            echo "${chunk_id},\$PWD/cohort_genomic_data/final_cohort_genomic_data_to_be_scored_chr\${chr_num},\${chr_num},bfile" >> samplesheet.csv
        done

        to_bind=(${projectDir} ${params.temp_dir})
        BIND_DIRS=""
        for path in \${to_bind[@]}; do
            BIND_DIRS+="\${path}:\${path},"
        done
        BIND_DIRS=\${BIND_DIRS::-1}

        SIF_IMAGE="${projectDir}/misc/pgsc_calc/pgsc_calc_container.sif"

        echo 'Loading singularity...'
        singularity exec --bind \$BIND_DIRS \$SIF_IMAGE bash <<END
            echo 'Initiating pgsc_calc nextflow pipeline...'
            nextflow run /opt/pgsc_calc/main.nf \\
                -profile conda \\
                --outdir "${projectDir}/outputs/pgsc_calc_results/${platform}" \\
                --input "samplesheet.csv" \\
                --scorefile "${projectDir}/outputs/OPGS_files/${platform}/${chunk_id}/*" \\
                --target_build GRCh38 \\
                --run_ancestry "${projectDir}/misc/pgsc_calc/pgsc_HGDP+1kGP_v1.tar.zst"
        END
        """
}

process getRefPanelPlatekeys {
    storeDir 'outputs/'

    input:
        val  first_pgs_chunk_id
        tuple path(final_cohort_platekeys), path(case_platekeys)

    output:
        path "reference_panel_platekeys.txt"

    script:
        """
        set -eux

        chunk=\$(grep -Po "^[a-zA-Z]+(?=_)" <<< ${first_pgs_chunk_id})
        pgs_file="${first_pgs_chunk_id}/score/\${chunk}_pgs.txt.gz"

        count-pgs() {
            zcat \$1 | tail -n +2 | awk '{print \$4}' | uniq | wc -l
        }

        first_OPGS=\$(zgrep -m 1 "OPGS.*" \${pgs_file} | awk '{print \$4}')
        n_pgs=\$(count-pgs \${pgs_file})
        if [[ \$n_pgs == 1 ]]; then
            n_individuals_incl_ref_panel=\$(( \$(zcat \${pgs_file} | wc -l) - 1 ))
        elif [[ \$n_pgs > 1 ]]; then
            n_individuals_incl_ref_panel=\$(zgrep -v -m 2 -n \${first_OPGS} \${pgs_file} | awk -F ':' '{print \$1}' | tail -n 1)
        else
            echo "ERROR: no individuals found in PGS file \${pgs_file}"
            exit 1
        fi

        awk '{print \$1}' ${final_cohort_platekeys} > final_cohort_platekeys_single_col.txt

        zcat \${pgs_file} | \\
        head -n \${n_individuals_incl_ref_panel} | \\
        grep -v -f final_cohort_platekeys_single_col.txt | \\
        awk '{print \$2}' | \\
        uniq > reference_panel_platekeys.txt
        """
}
