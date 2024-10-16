################################## SVA Batch Correction ####################################

library(limma)  # For removeBatchEffect
library(sva)    # For sva
library(dplyr)  # For data manipulation

# Create directory for adjusted data
dir.create("./data/processed/expression/adjusted_sva", recursive = TRUE)

# Define the path to the folder containing normalized count values
data_path <- "./data/processed/expression/readcounts_tmm/"

# Define the path to your metadata file
metadata_path <- "./data/processed/attphe.rds"

# Load metadata
metadata <- readRDS(metadata_path)

# List all .rds files under the data path
tissue_files <- list.files(path = data_path, pattern = "*.rds", full.names = TRUE)

# Define a function to process each tissue file
process_tissue <- function(tissue_file, metadata) {
  
  # Extract tissue name
  tissue_name <- gsub(".rds$", "", basename(tissue_file))
  
  # Load the normalized read counts data from the .rds file
  normalized_counts <- readRDS(tissue_file)
  
  # Extract sample IDs from expression data column names
  sample_ids <- colnames(normalized_counts)
  
  # Filter metadata to include only the matching samples
  attr_filtered <- metadata %>% filter(sample_id %in% sample_ids)
  
  # Check the number of samples
  if (ncol(normalized_counts) != nrow(attr_filtered)) {
    stop(paste("Number of samples in", tissue_name, "does not match number of rows in filtered metadata."))
  }
  
  # Rebuild model matrices with filtered metadata
  mod <- model.matrix(~ sex, data = attr_filtered)  # Model including sex
  mod0 <- model.matrix(~ 1, data = attr_filtered)   # Null model
  
  # Ensure that the rows of mod and mod0 correspond to the columns of normalized_counts
  rownames(mod) <- sample_ids
  rownames(mod0) <- sample_ids
  
  # Perform SVA to estimate surrogate variables
  sva_results <- sva(as.matrix(normalized_counts), mod, mod0)
  
  # Extract surrogate variables
  sv <- sva_results$sv
  
  print(sv)
  
  # Remove batch effects using the surrogate variables but keep sex effects
  adjusted_expression_data <- removeBatchEffect(
    as.matrix(normalized_counts),
    covariates = sv
  )
  
  # Convert back to data frame if needed
  adjusted_expression_data <- as.data.frame(adjusted_expression_data)
  
  # Define the file path to save the adjusted data
  result_file <- file.path("./data/processed/expression/adjusted_sva", paste0(tissue_name, ".rds"))
  
  saveRDS(adjusted_expression_data, file = result_file)
  cat("Processed tissue:", tissue_name, "\n")
  cat("Dimensions of adjusted data:", dim(adjusted_expression_data), "\n")
}

# Loop through each tissue file and process it
for (tissue_file in tissue_files) {
  process_tissue(tissue_file, metadata)
}
