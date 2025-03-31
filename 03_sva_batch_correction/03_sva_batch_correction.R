# sva_batch_correction.R

# Load required libraries
library(limma)    # For removeBatchEffect
library(sva)      # For sva
library(dplyr)    # For data manipulation

# Get the output directory from the environment variable
processed_dir <- Sys.getenv("PROCESSED_DIR")

metadata_path <- file.path(processed_dir, "attphe_all.rds")  # Corrected metadata path

# Load metadata
metadata <- readRDS(metadata_path)

# List all .rds files under the processed expression data path
tissue_files <- list.files(path = file.path(processed_dir, 'expression/readcounts_tmm_all/'), pattern = '*.rds', full.names = TRUE)

# Create directory for adjusted data (only if needed)
dir.create(file.path(processed_dir, 'expression/adjusted_sva_all'), recursive = TRUE, showWarnings = FALSE)

# Define sex-specific tissues (SVA should be skipped)
sex_specific_tissues <- c("Cervix-Ectocervix", "Cervix-Endocervix", "FallopianTube",
                          "Testis", "Uterus", "Vagina", "Ovary", "Prostate", "Breast-MammaryTissue")

# Function to process each tissue file
process_tissue <- function(tissue_file, metadata) {
  tissue_name <- gsub('.rds$', '', basename(tissue_file))
  
  # Load the normalized read counts data
  normalized_counts <- readRDS(tissue_file)
  
  # Ensure matrix format
  normalized_counts <- as.matrix(normalized_counts)
  
  # Extract sample IDs
  sample_ids <- colnames(normalized_counts)
  
  # Filter metadata to include only matching samples
  attr_filtered <- metadata %>% filter(sample_id %in% sample_ids)
  
  # Check the number of samples
  num_samples <- ncol(normalized_counts)
  
  if (num_samples < 20) {
    cat('Skipping tissue:', tissue_name, '(Insufficient samples:', num_samples, ')\n')
    return(NULL)  # Skip this file
  }
  
  # Define model matrices
  mod <- model.matrix(~ sex, data = attr_filtered)

  rownames(mod) <- sample_ids

  # Check if SVA should be skipped
  if (tissue_name %in% sex_specific_tissues) {
    cat("Skipping SVA for tissue:", tissue_name, "\n")
    # Apply batch effect removal using only sex as a covariate
    adjusted_expression_data <- removeBatchEffect(normalized_counts, batch = attr_filtered$batch1)
  } else {
    mod0 <- model.matrix(~ 1, data = attr_filtered)
    rownames(mod0) <- sample_ids
  
  
  # Estimate number of surrogate variables
  num_svs <- num.sv(normalized_counts, mod, method = 'be')
  
  # Perform SVA with error handling
  sva_results <- sva(normalized_counts, mod, mod0, n.sv = num_svs)
  
  # Extract surrogate variables and remove batch effects
  sv <- sva_results$sv
  adjusted_expression_data <- removeBatchEffect(normalized_counts, covariates = sv)
  
}
  # Convert back to data frame
  adjusted_expression_data <- as.data.frame(adjusted_expression_data)
  
  # Define the file path to save the adjusted data
  result_file <- file.path(processed_dir, 'expression/adjusted_sva_all', paste0(tissue_name, '.rds'))
  
  # Save only successfully processed tissues
  saveRDS(adjusted_expression_data, file = result_file)
  
  # Print confirmation
  cat('Processed tissue:', tissue_name, '\n')
  cat('Dimensions of adjusted data:', dim(adjusted_expression_data), '\n')
  }


# Process each tissue file
for (tissue_file in tissue_files) {
  process_tissue(tissue_file, metadata)
}



