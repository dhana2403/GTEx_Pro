############################
# Data Normalization by TMM
############################


library(edgeR)
library(dplyr)

# Get the output directory from the environment variable
processed_dir <- Sys.getenv("PROCESSED_DIR")


# Check if `readcounts_all/` exists
readcounts_dir <- file.path(processed_dir, 'expression/readcounts_all/')
if (!dir.exists(readcounts_dir)) {
  stop("Error: 'expression/readcounts_all/' directory does not exist!")
}

# Initialize an empty list to collect sample counts
sample_counts_list <- list()

# Quality check function for read counts with edgeR TMM normalization
qc_check_readcounts <- function(tis, outliers = NA) {
  print(paste("Processing tissue:", tis))  # Debugging
  
  # Load attributes and read counts data
  att_path <- file.path(processed_dir, 'attphe_all.rds')
  readcounts_path <- file.path(readcounts_dir, tis)

  if (!file.exists(att_path)) {
    stop(paste("Error: Metadata file not found at", att_path))
  }
  if (!file.exists(readcounts_path)) {
    stop(paste("Error: Read counts file not found for", tis))
  }

  att <- readRDS(att_path)
  readcounts <- readRDS(readcounts_path)
  
  # Filter sample names and read counts excluding outliers
  sampnames <- setdiff(intersect(colnames(readcounts), att$sample_id), outliers)
  att_filtered <- att %>% filter(sample_id %in% sampnames)
  
  if (length(sampnames) == 0) {
    print(paste("Skipping", tis, "because no valid samples were found."))
    return(NULL)
  }

  readcounts_filtered <- readcounts[, sampnames]
  
  # Perform TMM normalization using edgeR
  dge <- DGEList(counts = readcounts_filtered) %>%
    calcNormFactors()
  readcounts_norm <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
  
  output_dir <- file.path(processed_dir, 'expression/readcounts_tmm_all')
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  output_path <- file.path(output_dir, tis)
  saveRDS(readcounts_norm, output_path)
  
  # Append the sample count to the list
  sample_counts_list[[tis]] <<- data.frame(
    Tissue = tis,
    SampleCount = length(sampnames),
    stringsAsFactors = FALSE
  )
  
  print(paste("Normalized data saved for:", tis, "at", output_path))
}

# Get list of tissue files
alltis <- list.files(readcounts_dir)

if (length(alltis) == 0) {
  stop("Error: No tissue files found in 'expression/readcounts_all/'")
}

# Apply the QC function to all tissues
sapply(alltis, qc_check_readcounts)

# Check if normalization was successful
if (length(sample_counts_list) == 0) {
  stop("Error: No tissues were processed. Check the input files!")
}

# Combine the results and save as a CSV
global_sample_counts_df <- do.call(rbind, sample_counts_list)
write.csv(global_sample_counts_df, file.path(processed_dir, 'sample_counts.csv'), row.names = FALSE)

rm(list = ls())

