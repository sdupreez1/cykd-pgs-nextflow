process platformFDRandPlotting {
    publishDir { "outputs/comparison_results/stats_and_figures/${platform}" }, mode: 'copy'

    input:
        tuple val(platform), path(platform_pvals, stageAs: "platform_pvals/*")
    
    output:
        path "${platform}_significant_pgs.tsv", emit: sig_pgs

    script:
        """
        module load R/4.5.2
        platform_fdrs.R ${platform} ${projectDir}
        """
}
