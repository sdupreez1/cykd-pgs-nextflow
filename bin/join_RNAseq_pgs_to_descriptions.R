#!/bin/env Rscript
job_args = commandArgs(trailingOnly=TRUE)
platchunk = job_args[1]

.libPaths(c(.libPaths(), "/tools/aws-workspace-ubuntu-apps/ce/R/4.5.2"))
library(tidyr)
library(dplyr)
library(withr)

named_pgs = read.table(paste0(job_args, "_traits_without_functional_name.tsv"),
                       sep='\t', 
                       header=TRUE, 
                       quote="") %>%
  rename(name=trait)
descriptions = read.table(paste0(job_args, "_named_descriptions.tsv"),
                          sep='\t',
                          header=TRUE,
                          quote="")

traits = left_join(named_pgs, descriptions, by="name") %>% mutate(trait=paste0(name, ": ", trait)) %>% select(PGS, trait)

write.table(traits, paste0(job_args, "_traits.tsv"), sep='\t', row.names=FALSE, quote=FALSE)
print(paste("Successfully saved", job_args, "descriptions"))
