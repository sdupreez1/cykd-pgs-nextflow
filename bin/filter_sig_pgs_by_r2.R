#!/bin/env Rscript
job_args = commandArgs(trailingOnly=TRUE)
projectDir = job_args[1]

.libPaths(c(.libPaths(), "/tools/aws-workspace-ubuntu-apps/ce/R/4.5.2"))
library(readxl)
library(dplyr)

sig_scorefiles = read.table("combined_sig_pgs.tsv",
                            sep='\t',
                            header=TRUE)
print('sig_scorefiles:')
print(sig_scorefiles)

supp_table_path = paste0(projectDir, "/misc/OmicsPred_supp_tables.xlsx")
metabolon_summary = read_excel(supp_table_path,
                                sheet=1,
                                skip=2,
                                .name_repair="universal")
metabolon_summary = metabolon_summary[,c(1:15)]
metabolon_r2 = metabolon_summary[,c("OMICSPRED.ID", "Internal_R2")] %>%
    cbind(platform=rep("Metabolon", nrow(metabolon_summary)))

nightingale_summary = read_excel(supp_table_path,
                                sheet=2,
                                skip=2,
                                .name_repair="universal")
nightingale_summary = nightingale_summary[,c(1:10)]
nightingale_r2 = nightingale_summary[,c("OMICSPRED.ID", "Internal_R2")] %>%
    cbind(platform=rep("Nightingale", nrow(nightingale_summary)))

olink_summary = read_excel(supp_table_path,
                                sheet=3,
                                skip=2,
                                .name_repair="universal")
olink_summary = olink_summary[,c(1:9)]
olink_r2 = olink_summary[,c("OMICSPRED.ID", "Internal_R2")] %>%
    cbind(platform=rep("Olink", nrow(olink_summary)))

somascan_summary = read_excel(supp_table_path,
                                sheet=4,
                                skip=2,
                                .name_repair="universal")
somascan_summary = somascan_summary[,c(1:11)]
somascan_r2 = somascan_summary[,c("OMICSPRED.ID", "Internal_R2")] %>%
    cbind(platform=rep("SomaScan", nrow(somascan_summary)))

rnaseq_summary = read_excel(supp_table_path,
                                sheet=5,
                                skip=2,
                                .name_repair="universal")
rnaseq_summary = rnaseq_summary[,c(1:8)]
rnaseq_r2 = rnaseq_summary[,c("OMICSPRED.ID", "Internal_R2")] %>%
    cbind(platform=rep("RNAseq", nrow(rnaseq_summary)))


# Add UKB sig scorefile data below from OmicsPred.org
ukb_r2 = data.frame(
    OMICSPRED.ID = c(),
    Internal_R2  = c(),
    platform     = c()
)

all_r2 = rbind(olink_r2, metabolon_r2, rnaseq_r2, somascan_r2, nightingale_r2, ukb_r2)
print('all_r2:')
print(all_r2)

sig_scorefiles_r2 = left_join(
    sig_scorefiles, all_r2,
    by=join_by(`PGS` == `OMICSPRED.ID`)
)
valid_sig_scores = sig_scorefiles_r2 %>%
    filter(Internal_R2 > 0.01)

associated_scorefiles = valid_sig_scores %>%
    filter(hazard_ratio != 1.00) %>%
    unique()

write.table(associated_scorefiles,
            "final_significant_pgs.tsv",
            quote=FALSE,
            row.names=FALSE,
            sep='\t'
)
