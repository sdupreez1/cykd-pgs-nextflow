#!/usr/bin/env Rscript

.libPaths(c(.libPaths(), "/tools/aws-workspace-ubuntu-apps/ce/R/4.5.2"))
library(httr)
library(jsonlite)
library(curl)
library(Rlabkey)
library(data.table)
library(dplyr)
library(stringr)

labkey.setWafEncoding(FALSE) # Needed for sql query to work
labkey.setDefaults(baseUrl="https://labkey-embassy.gel.zone/labkey/")

all_unrelated_platekeys = read.table(
  'ancestrally_matched_cases_and_controls_platekeys.txt',
  header=FALSE,
  sep='\t'
)$V1

cases_platekeys = read.table(
  'ancestrally_matched_cases_platekeys.txt',
  header=FALSE,
  sep='\t'
)$V1

#==================
#     Controls
#==================

query_controls <- "SELECT rd.participant_id, rd.participant_type, par.participant_phenotyped_sex, par.rare_disease_diagnosis_age, rd.plate_key, par.death_date, par.yob, part.latest_case_activity_datetime
FROM rare_disease_analysis AS rd
LEFT JOIN participant_summary AS par
ON rd.participant_id = par.participant_id
LEFT JOIN rare_diseases_pedigree_member AS rdp
ON rd.participant_id = rdp.participant_id
LEFT JOIN genome_file_paths_and_types AS sq
ON rd.participant_id = sq.participant_id
LEFT JOIN participant AS part
ON rd.participant_id = part.participant_id
WHERE rd.participant_type = 'Relative'
AND rdp.affection_status = 'Unaffected'
AND sq.genome_build = 'GRCh38'"

mysql_controls <- labkey.executeSql(
  schemaName="lists",                                                 
  colNameOpt = "rname",                                              
  maxRows = 100000000,                                              
  folderPath="/main-programme/main-programme_v19_2024-10-31",      
  sql = query_controls                                            
)

# All participants admitted as unaffected relatives of 100KGP probands
# We then filter for people admitted for non-kidney problems
mysql_controls <-data.frame(mysql_controls)
y <- unique(mysql_controls)

controls <- subset(y, y$plate_key %in% all_unrelated_platekeys) %>%
  mutate(rare_disease_diagnosis_age = NA,
         last_updated_year = as.integer(str_sub(latest_case_activity_datetime, start=1, end=4)),
         death_year = as.integer(str_sub(death_date, start=1, end=4)),
         year_last_seen = pmin(death_year, last_updated_year, na.rm=TRUE), # pmin() does for each element, min() does for whole field
         age_to_use = year_last_seen - (yob + 1) # assume dob is 31st Dec (lower bound for age)
  )

#==============
#    Cases
#==============

query_cases <- "SELECT rd.participant_id, rd.participant_type, par.participant_phenotyped_sex, par.rare_disease_diagnosis_age, rd.plate_key, par.death_date, par.yob, part.latest_case_activity_datetime
FROM rare_disease_analysis AS rd
LEFT JOIN participant_summary AS par
ON rd.participant_id = par.participant_id
LEFT JOIN participant AS part
ON rd.participant_id = part.participant_id
WHERE rd.normalised_specific_disease = 'Cystic kidney disease'
AND rd.participant_type = 'Proband'"


mysql_cases <- labkey.executeSql(
  schemaName="lists",                                                 
  colNameOpt = "rname",                                              
  maxRows = 100000000,                                              
  folderPath="/main-programme/main-programme_v19_2024-10-31",      
  sql = query_cases                                               
)

x <-data.frame(mysql_cases)
x <- x[!duplicated(x$participant_id),]

cases <-subset(x, x$plate_key %in% cases_platekeys) %>%
  mutate(last_updated_year = as.integer(str_sub(latest_case_activity_datetime, start=1, end=4)),
         death_year = as.integer(str_sub(death_date, start=1, end=4)),
         year_last_seen = pmin(death_year, last_updated_year, na.rm=TRUE), 
         age_to_use = rare_disease_diagnosis_age
  )

#================
#     Stats
#================

whole_cohort = rbind(cases, controls) %>%
  filter(participant_phenotyped_sex != "<NA>") 

stats = whole_cohort %>% 
  summarise(n_cohort = n(),
            n_cases = length(cases_platekeys),
            n_controls = n() - length(cases_platekeys),
            mean_age_used = mean(age_to_use, na.rm=TRUE),
            age_used_sd = sd(age_to_use, na.rm=TRUE),
            diagnosis_mean_age = mean(rare_disease_diagnosis_age, na.rm=TRUE),
            diagnosis_age_sd = sd(rare_disease_diagnosis_age, na.rm=TRUE),
            percent_female = 100*sum(participant_phenotyped_sex=="Female", na.rm=TRUE)/n()
            )

write.table(whole_cohort[,c("plate_key", "participant_phenotyped_sex", "age_to_use")],
            "full_final_cohort_platekey_sex_age.tsv",
            sep='\t',
            quote=FALSE,
            row.names=FALSE)
write.table(stats,
            "final_cohort_summary_stats.tsv",
            sep='\t',
            quote=FALSE,
            row.names=FALSE)
