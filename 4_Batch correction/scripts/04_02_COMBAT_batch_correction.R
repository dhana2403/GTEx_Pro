
################################## COMBAT BATCH CORRECTION ####################################

library(limma)  # For ComBat
library(sva)
library(dplyr)  # For data manipulation

dir.create("./data/processed/expression/adjusted_combat", recursive = TRUE)

data_path <- "./data/processed/expression/readcounts_tmm/"

metadata_path <- "./data/processed/attphe.rds"

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
  
  # Fit a linear model to remove the effect of sex
  mod_sex <- model.matrix(~ sex, data = attr_filtered)
  residuals <- apply(normalized_counts, 1, function(gene) {
    lm(gene ~ mod_sex - 1)$residuals
  })
  
  # Convert residuals to matrix
  residuals <- t(residuals)
  
  # Extract batch information from metadata
  batch <- attr_filtered$batch1  
  
  # Ensure that batch is a factor
  batch <- as.factor(batch)
  
  # Apply ComBat for batch correction on residuals
  adjusted_expression_data <- ComBat(
    dat = residuals,
    batch = batch
  )
  
  # Convert back to data frame if needed
  adjusted_expression_data <- as.data.frame(adjusted_expression_data)
  
  # Define the file path to save the adjusted data
  result_file <- file.path("./data/processed/expression/adjusted_combat", paste0(tissue_name, '.rds'))
  
  saveRDS(adjusted_expression_data, file = result_file)
  
  cat("Processed tissue:", tissue_name, "\n")
  cat("Dimensions of adjusted data:", dim(adjusted_expression_data), "\n")
}

# Loop through each tissue file and process it
for (tissue_file in tissue_files) {
  process_tissue(tissue_file, metadata)
}
