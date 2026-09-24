#!/bin/env Rscript
job_args = commandArgs(trailingOnly=TRUE)
platform   = job_args[1]
projectDir = job_args[2]

.libPaths(c(.libPaths(), "/tools/aws-workspace-ubuntu-apps/ce/R/4.5.2"))
library(tidyr)
library(dplyr)
library(withr)     
library(labeling) 
library(farver)   
library(ggplot2)
library(survival)

pgs_pvals = data.frame()
pval_file_list = list.files("platform_pvals/", full.names=TRUE)
print("pval_file_list:")
print(pval_file_list)

for (file in pval_file_list) {
    current_pvals = read.table(file, sep='\t', header=TRUE)
    pgs_pvals = rbind(pgs_pvals, current_pvals)
}
print("pgs_pvals:")
head(pgs_pvals)

pgs_pvals = pgs_pvals %>%
  mutate(FDR=p.adjust(`wald_pval`, method="fdr")) 

print("pgs_pvals:")
head(pgs_pvals[,c("PGS","wald_pval","FDR","trait","hazard_ratio","HR_CI")])

sig_lvl = 0.05
sig_pgs = pgs_pvals[,c("PGS","wald_pval","FDR","trait","hazard_ratio","HR_CI")] %>%
  filter(FDR < sig_lvl)

print("sig_pgs:")
head(sig_pgs)

print("Saving significant scorefiles")
write.table(sig_pgs,
            paste0(platform, "_significant_pgs.tsv"),
            quote=FALSE,
            row.names=FALSE,
            sep='\t')

