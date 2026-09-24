#!/bin/env Rscript
job_args = commandArgs(trailingOnly=TRUE)
platform       = job_args[1]
platform_chunk = job_args[2]
projectDir     = job_args[3]

.libPaths(c(.libPaths(), '/tools/aws-workspace-ubuntu-apps/ce/R/4.5.2'))
library(tidyr)
library(dplyr)
library(survival)
library(stringr)

print('Packages loaded successfully')

#=============================
#       Loading data
#=============================
print('')
print(paste('### LOADING', platform_chunk, 'DATA ###'))

print('Loading platekeys')
case.platekeys          = read.table('ancestrally_matched_cases_platekeys.txt')$V1
full_cohort.platekeys   = read.table('ancestrally_matched_cases_and_controls_platekeys.txt')$V1
ref_panel.platekeys     = read.table('reference_panel_platekeys.txt', header=TRUE)$FID

control.platekeys       = setdiff(full_cohort.platekeys, ref_panel.platekeys)
all.platekeys           = c(full_cohort.platekeys, ref_panel.platekeys)
n.platekeys             = length(all.platekeys)

full_cohort = read.table('full_final_cohort_platekey_sex_age.tsv', 
                         sep='\t', header=TRUE) %>%
    rename(IID=plate_key,
            sex=participant_phenotyped_sex,
            age=age_to_use) %>%
    mutate(sex=as.factor(sex))

# Disabled quoting (quote="") since some traits mentions 5' which gives "EOF within quoted string" error when parsed
scorefile_traits  = read.table(
    paste0(platform_chunk, '_traits.tsv'),
    sep='\t', header=TRUE, quote=""
    )
print('scorefile_traits:')
head(scorefile_traits)
    
chunk.path = paste0(projectDir, '/outputs/pgsc_calc_results/', platform, '/', platform_chunk, '/score/')
pca_results = read.table(
  paste0(chunk.path, platform_chunk, '_popsimilarity.txt.gz'), 
  sep='\t', header=TRUE)

ancestry_norm_scores_cols = read.table(
  paste0(chunk.path, platform_chunk, '_pgs.txt.gz'), 
  header=TRUE, 
  nrows=1)

ancestry_norm_scores = read.table(
  paste0(chunk.path, platform_chunk, '_pgs.txt.gz'), 
  sep='\t', 
  header=TRUE
)
colnames(ancestry_norm_scores) <- colnames(ancestry_norm_scores_cols)

# PGS file also has reference panel scores, so we filter for only CyKD cohort platekeys
ancestry_norm_scores = ancestry_norm_scores %>% 
  left_join(pca_results[,c(3:13)], by='IID') %>%
  filter(IID %in% full_cohort.platekeys) %>% 
  mutate(CASE = IID %in% case.platekeys) %>%
  left_join(full_cohort, by='IID')

print('ancestry_norm_scores:')
head(ancestry_norm_scores)

scorefiles_used = unique(ancestry_norm_scores$PGS)
print(paste('Number of scorefiles (from score data):', length(scorefiles_used)))

mean_sd = ancestry_norm_scores %>% 
  filter(PGS == ancestry_norm_scores$PGS[1]) %>%
  select(SUM, Z_MostSimilarPop, Z_norm1, Z_norm2) %>% 
  summarise_all( list(mean = mean, sd= sd) ) 

print(paste('Finished loading', platform_chunk, 'scoring data'))

#======================================================================
#   Testing differences in distributions of case and control PGSs
#======================================================================
print('')
print('### FITTING COX PH MODEL TO SCORING DATA ###')

pgs_wald_pvals = data.frame() # Empty data frame to add data to via below
PCs = c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
for (OPGS in scorefiles_used) {
  single_pgs_scores = filter(ancestry_norm_scores, PGS==OPGS)
  single_pgs_scores[,PCs] = as.data.frame(scale(single_pgs_scores[,PCs]))
  
  fitted.cox = coxph(
    formula = Surv(time=age, event=CASE) ~ Z_norm2 + strata(sex) + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
    data = single_pgs_scores
    )
  pgs_coef_pval = as.data.frame(summary(fitted.cox)$coefficients)[1,5]
  pgs_hazard_ratio = round(summary(fitted.cox)$conf.int[1,1], 2)
  pgs_hazard_ration_CI_lower = round(summary(fitted.cox)$conf.int[1,3], 2)
  pgs_hazard_ration_CI_upper = round(summary(fitted.cox)$conf.int[1,4], 2)
  
  cox_results = data.frame(PGS=OPGS, 
                           wald_pval= pgs_coef_pval,
                           hazard_ratio = pgs_hazard_ratio,
                           HR_CI = paste0('[', pgs_hazard_ration_CI_lower, ', ', pgs_hazard_ration_CI_upper, ']')
                           )
  
  pgs_wald_pvals = pgs_wald_pvals %>%
    rbind(cox_results)
}
pgs_wald_pvals = left_join(pgs_wald_pvals, scorefile_traits[,c('PGS','trait')], by='PGS')

print('pgs_wald_pvals:')
head(pgs_wald_pvals)

write.table(pgs_wald_pvals,
            paste0(platform_chunk, '_pvals.tsv'),
            quote=FALSE,
            row.names=FALSE,
            sep='\t')
