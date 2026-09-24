process process formatHarmonisedOPGS {
    publishDir { "outputs/OPGS_traits/${platform}/" }, mode: 'copy', pattern: "*to*_traits.tsv"

    input:
        tuple(
            val(chunk_id),
            val(platform),
            path(unformatted_files, stageAs: 'unformatted_files/*')
        )

    // Don't emit files since pgsc_calc struggles to deal with symlinked files in calculatePGS
    output:
        tuple val("${chunk_id}"), val("${platform}"), emit: chunk_id
        path "${platform}*to*_traits.tsv"           , emit: traits

    script:
        """
        echo 'platform = ${platform}'
        echo 'unformatted_files'
        echo "${unformatted_files}"

        formatted_dir="${projectDir}/outputs/OPGS_files/${platform}/${chunk_id}"
        mkdir -p \${formatted_dir}

        echo -e "rsID\\tchr_name\\tchr_position\\teffect_allele\\tother_allele\\teffect_weight\\thm_source\\thm_rsID\\thm_chr\\thm_pos\\thm_inferOtherAllele" > col_titles_hm.txt

        # Formatting
        for OPGS in ${unformatted_files};
        do
            gz_file=0
            grep -q ".gz" <<< \$OPGS && \\
            gz_file=1 && filename=\$(basename \$OPGS .gz) || filename=\$(basename \$OPGS)

            if [[ \$gz_file == 1 ]]; then
                gunzip -c \$OPGS > unformatted_files/\$filename
            fi

            grep "^#" unformatted_files/\$filename > formating_temp.txt
            awk -v OFS=\$'\\t' '{print \$1, \$2, \$3, \$4, \$5, \$6}' col_titles_hm.txt >> formating_temp.txt

            grep -e "ENSEMBL" -e "liftover" unformatted_files/\$filename > body_temp.txt

            head_lc=\$(wc -l < formating_temp.txt)
            orig_lc=\$((\$(wc -l < unformatted_files/\$filename)-\$head_lc))  
            filtered_lc=\$(wc -l < body_temp.txt)

            echo "Removed" \$((\${orig_lc}-\${filtered_lc})) "variants from" \$filename

            if [ \$(wc -l < body_temp.txt) == 0 ]; then
                echo "NO VALID BODY FOUND FOR" \$filename
            else
                > formatted_body_temp.txt
                        
                while IFS= read -r line; do
                    if [ \$(echo \$line | grep -c liftover) == 1 ]; then
                        awk -v OFS=\$'\\t' '{print \$1, \$2, \$9, \$4, \$5, \$6}' <<< \$line >> formatted_body_temp.txt
                    else 
                        awk -v OFS=\$'\\t' '{print \$1, \$2, \$10, \$4, \$5, \$6}' <<< \$line >> formatted_body_temp.txt
                    fi
                done < body_temp.txt

                cat formatted_body_temp.txt >> formating_temp.txt

                sed -i -e 's/#genome_build=GRCh37/#genome_build=GRCh38/' -e 's/#omicspred_id/#pgs_id/' formating_temp.txt       
                cat formating_temp.txt > \${formatted_dir}/formatted_${platform}_\${filename}

            fi
        done

        # Get traits
        if [[ ${platform} != 'RNAseq' ]]; then
            
            # awk program structure:
            #    (for each file in the platform OPGS directory)
            #    - Print column titles; 
            #    - Store pgs_id and trait_reported of the current scorefile; 
            #    - Print these values in the same row

            awk -v OFS='\\t' --field-separator '=' \\
            'BEGIN {print "PGS", "trait"}; \\
            /pgs_id/ {id=\$2}; \\
            /trait_reported/ {trait=\$2}; \\
            ENDFILE {print id, trait}' \\
            \${formatted_dir}/* > ${chunk_id}_traits.tsv
        else
            awk -v OFS='\\t' --field-separator '=' \\
            'BEGIN {print "PGS", "trait"}; \\
            /pgs_id/ {id=\$2}; \\
            /trait_reported/ {trait=\$2}; \\
            ENDFILE {print id, trait}' \\
            \${formatted_dir}/* > ${chunk_id}_traits_without_functional_name.tsv

            # Get all transcript names (of the scorefiles which have names listed) to feed into zgrep, and make them match only the exact name
            tail -n +2 ${chunk_id}_traits_without_functional_name.tsv | \\
            awk -F '\\t' '{print \$2}'| \\
            grep '^\\w' | \\
            sed -e 's/^/;Name=/' -e 's/\$/;/' > ${chunk_id}_names_only.txt 

            echo -e 'name\\ttrait' > ${chunk_id}_named_descriptions.tsv

            zgrep -f ${chunk_id}_names_only.txt ${params.ensembl_GRCh38_gff3} | \\
            awk -F '\\t' '{print \$9}' | \\
            awk -F ';' '{for(i=1;i<=NF;i++)if(\$i~/Name=|description=/)printf("%s;", \$i)}{print ""}' | \\
            sed -e 's/=/\\t/g' -e 's/;/\\t/' -e 's/\\[Source:/\\t/' | \\
            awk -v OFS='\\t' -F '\\t' '{print \$2,\$4}'>> ${chunk_id}_named_descriptions.tsv

            module load R/4.5.2
            join_RNAseq_pgs_to_descriptions.R "${chunk_id}"
        fi
        """
}
