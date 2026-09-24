#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
which_cohort = args[1]
provide_participant_ids = args[2]

.libPaths(c(.libPaths(), "/tools/aws-workspace-ubuntu-apps/ce/R/4.5.2"))
library(httr)
library(jsonlite)
library(curl)
library(Rlabkey)
library(data.table)

if (which_cohort == "control") {
  labkey.setWafEncoding(FALSE)
  labkey.setDefaults(baseUrl="https://labkey-embassy.gel.zone/labkey/")

  query <- "SELECT rd.plate_key, hp.hpo_term, hp.hpo_present, rd.participant_id
  FROM rare_disease_analysis AS rd
  LEFT JOIN genome_file_paths_and_types AS sq
  ON rd.participant_id = sq.participant_id
  LEFT JOIN rare_diseases_pedigree_member AS rdp
  ON rd.participant_id = rdp.participant_id
  LEFT JOIN participant AS par
  ON rd.participant_id = par.participant_id
  LEFT JOIN rare_diseases_participant_phenotype AS hp
  ON rd.participant_id = hp.participant_id
  WHERE rd.participant_type = 'Relative'
  AND rdp.affection_status = 'Unaffected'
  AND sq.genome_build = 'GRCh38'"

  mysql <- labkey.executeSql(
    schemaName = "lists",                                               
    colNameOpt = "rname",                                               
    maxRows = 100000000,                                                
    folderPath="/main-programme/main-programme_v19_2024-10-31",         
    sql = query                                                         
  )

  mysql <-data.frame(mysql)
  y <- unique(mysql)

  a <- y[which(y$hpo_present=="Yes"),]
  b <- a[grep("renal", a$hpo_term,ignore.case = T),]
  c <- a[grep("kidney", a$hpo_term,ignore.case = T),]
  d <- a[grep("hematuria", a$hpo_term,ignore.case = T),]
  e <- a[grep("proteinuria", a$hpo_term,ignore.case = T),]
  g <- a[grep("glom", a$hpo_term,ignore.case = T),]
  f <- rbind(b,c,d,e,g)
  f <- as.data.frame(f$plate_key)

  no_hpo_full <- subset(y, !y$plate_key %in% f$plate_key)
  no_hpo <- no_hpo_full$plate_key
  no_hpo <- as.data.frame(unique(no_hpo))

  ckd_to_remove <- fread('participant_explorer_all_mention_of_cystic_platekeys_july2022.csv', header = F)

  no_hpo_no_ckd <- subset(no_hpo, !no_hpo$`unique(no_hpo)` %in% ckd_to_remove$V1)
  cystic_to_remove <- fread('controls_with_renal_disease_or_ckd.tsv', header = F)
  no_hpo_no_ckd_no_cystic <- subset(no_hpo_no_ckd , !no_hpo_no_ckd $`unique(no_hpo)` %in% cystic_to_remove$V1)

  control_part_ids <- subset(no_hpo_full, no_hpo_full$plate_key %in% no_hpo_no_ckd_no_cystic$`unique(no_hpo)`)$participant_id

  no_hpo_no_ckd_no_cystic.for_plink = cbind(no_hpo_no_ckd_no_cystic$plate_key, no_hpo_no_ckd_no_cystic$plate_key)
  write.table(no_hpo_no_ckd_no_cystic, "related_control_platekeys.txt", col.names = F, row.names = F, quote = F)

  if (provide_participant_ids == 1) {
    write.table(control_part_ids, "related_control_participant_ids.txt", col.names = F, row.names = F, quote = F)
  }
  
} else if (which_cohort == "case") {
  labkey.setWafEncoding(FALSE)
  labkey.setDefaults(baseUrl = "https://labkey-embassy.gel.zone/labkey/")

  query <- "SELECT rd.participant_id, rd.participant_type, par.consanguinity, par.participant_phenotypic_sex, par.year_of_birth, par.mother_affected, par.father_affected, par.full_brothers_affected, par.full_sisters_affected, rd.normalised_specific_disease, rd.plate_key, gc.acmg_classification, gc.case_solved_family, gc.additional_comments, gc.gene_name, gc.assembly, gc.chromosome, gc.position, gc.reference, gc.alternate
  FROM rare_disease_analysis AS rd
  LEFT JOIN gmc_exit_questionnaire gc
  ON rd.participant_id = gc.participant_id
  LEFT JOIN participant AS par
  ON gc.participant_id = par.participant_id
  WHERE rd.normalised_specific_disease = 'Cystic kidney disease'
  AND rd.participant_type = 'Proband'"


  mysql <- labkey.executeSql(
    schemaName="lists",                                                 
    colNameOpt = "rname",                                               
    maxRows = 100000000,                                                
    folderPath="/main-programme/main-programme_v19_2024-10-31",         
    sql = query                                                        
  )

  y <-data.frame(mysql)
  y<- y[!duplicated(y$participant_id),]

  cykd <-subset(y, participant_type =='Proband' & normalised_specific_disease == "Cystic kidney disease")
  cykd[is.na(cykd)] <- "Unknown"

  cykd.for_plink = cbind(cykd$plate_key, cykd$plate_key)
  write.table(cykd.for_plink, "related_case_platekeys.txt", col.names = F, row.names = F, quote = F)

  if (provide_participant_ids == 1) {
    write.table(cykd$participant_id, "related_case_participant_ids.txt", col.names = F, row.names = F, quote = F)
  }

} else {
  stop("Unrecognised cohort name. Should only be either 'case' or 'control'")
}
