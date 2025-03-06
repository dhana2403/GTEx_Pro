# Load required packages
library(limma)    # For removeBatchEffect
library(sva)      # For sva
library(dplyr)    # For data manipulation

# Create directory for adjusted data (only if needed)
dir.create("./data/processed/expression/adjusted_sva_all_cross", recursive = TRUE, showWarnings = FALSE)

# Define paths
data_path <- "./data/processed/expression/readcounts_tmm_all_cross/"
metadata_path <- "./data/processed/attphe_all_cross.rds"

# Load metadata
metadata <- readRDS(metadata_path)

# List all .rds files under the data path
tissue_files <- list.files(path = data_path, pattern = "*.rds", full.names = TRUE)

# Function to process each tissue file
process_tissue <- function(tissue_file, metadata) {
  tissue_name <- gsub(".rds$", "", basename(tissue_file))
  
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
    cat("Skipping tissue:", tissue_name, "(Insufficient samples:", num_samples, ")\n")
    return(NULL)  # Skip this file
  }
  
  # Define model matrices
  mod <- model.matrix(~ sex, data = attr_filtered)
  mod0 <- model.matrix(~ 1, data = attr_filtered)
  
  rownames(mod) <- sample_ids
  rownames(mod0) <- sample_ids
  
  # Estimate number of surrogate variables
  num_svs <- num.sv(normalized_counts, mod, method = "be")
  

  sva_results <- sva(normalized_counts, mod, mod0, n.sv = num_svs)
 
  if (is.null(sva_results)) {
    cat("Skipping SVA for tissue:", tissue_name, "\n")
    return(NULL)  # Do not proceed to save the file
  }
  
  # Extract surrogate variables and remove batch effects
  sv <- sva_results$sv
  adjusted_expression_data <- removeBatchEffect(normalized_counts, covariates = sv)
  
  # Convert back to data frame
  adjusted_expression_data <- as.data.frame(adjusted_expression_data)
  
  # Define the file path to save the adjusted data
  result_file <- file.path("./data/processed/expression/adjusted_sva_all_cross", paste0(tissue_name, ".rds"))
  
  # Save only successfully processed tissues
  saveRDS(adjusted_expression_data, file = result_file)
  
  # Print confirmation
  cat("Processed tissue:", tissue_name, "\n")
  cat("Dimensions of adjusted data:", dim(adjusted_expression_data), "\n")
}

# Process each tissue file
for (tissue_file in tissue_files) {
  process_tissue(tissue_file, metadata)
}
