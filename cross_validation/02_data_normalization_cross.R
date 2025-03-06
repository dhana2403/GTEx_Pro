############################ Data normalization by TMM ############################

library(edgeR)
library(dplyr)

# Initialize an empty list to collect sample counts
sample_counts_list <- list()

# Quality check function for read counts with edgeR TMM normalization
qc_check_readcounts <- function(tis, outliers = NA) {
  
  # Load attributes and read counts data
  att <- readRDS('./data/processed/attphe_all_tis.rds')
  readcounts <- readRDS(paste0('./data/processed/expression/readcounts_all_cross/', tis))
  
  # Filter sample names and read counts excluding outliers
  sampnames <- setdiff(intersect(colnames(readcounts), att$sample_id), outliers)
  att_filtered <- att %>% filter(sample_id %in% sampnames)
  readcounts_filtered <- readcounts[, sampnames]
  
  # Perform TMM normalization using edgeR
  dge <- DGEList(counts = readcounts_filtered) %>%
    calcNormFactors()
  readcounts_norm <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
  
  output_dir <- './data/processed/expression/readcounts_tmm_all_cross/'
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(readcounts_norm, paste0(output_dir, tis))
  
  # Append the sample count to the list
  sample_counts_list[[tis]] <<- data.frame(
    Tissue = tis,
    SampleCount = length(sampnames),
    stringsAsFactors = FALSE
  )
  
  print(paste("Normalized data saved for:", tis))
}

alltis <- list.files('data/processed/expression/readcounts_all_cross/')

# Apply the QC function to all tissues
sapply(alltis, qc_check_readcounts)

# Combine the results and save as a CSV
global_sample_counts_df <- do.call(rbind, sample_counts_list)
write.csv(global_sample_counts_df, './data/processed/sample_counts.csv', row.names = FALSE)

rm(list = ls())