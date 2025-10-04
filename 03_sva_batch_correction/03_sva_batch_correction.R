# sva_batch_correction_academic.R
# Academic/demo version of GTEx_Pro
# Simplified batch correction for academic/research use

# Load required libraries
library(limma)    # For removeBatchEffect
library(dplyr)    # For data manipulation

# Get the output directory from the environment variable
processed_dir <- Sys.getenv("PROCESSED_DIR")

metadata_path <- file.path(processed_dir, "attphe_all.rds")  

# Load metadata
metadata <- readRDS(metadata_path)

# List all .rds files under the processed expression data path
tissue_files <- list.files(
  path = file.path(processed_dir, 'expression/readcounts_tmm_all/'),
  pattern = '*.rds', 
  full.names = TRUE
)

# Create directory for adjusted data
dir.create(file.path(processed_dir, 'expression/adjusted_sva_all'),
           recursive = TRUE, showWarnings = FALSE)

# Define sex-specific tissues (for special handling if needed)
sex_specific_tissues <- c("Cervix-Ectocervix", "Cervix-Endocervix", "FallopianTube",
                          "Testis", "Uterus", "Vagina", "Ovary", "Prostate", "Breast-MammaryTissue")

# -----------------------------
# Academic/demo batch correction
# -----------------------------
academic_batch_correction <- function(normalized_counts, metadata, tissue_name) {
  cat("Academic demo: running simplified correction for:", tissue_name, "\n")
  
  # Simple covariate-based adjustment for demo purposes
  # (Full SVA-based batch correction is proprietary and available under commercial license)
  if("sex" %in% colnames(metadata)) {
    adjusted <- limma::removeBatchEffect(normalized_counts, batch = metadata$sex)
  } else {
    adjusted <- normalized_counts
  }
  
  return(adjusted)
}

# -----------------------------
# Process each tissue file
# -----------------------------
process_tissue <- function(tissue_file, metadata) {
  tissue_name <- gsub('.rds$', '', basename(tissue_file))
  
  # Load normalized read counts
  normalized_counts <- readRDS(tissue_file)
  normalized_counts <- as.matrix(normalized_counts)
  
  sample_ids <- colnames(normalized_counts)
  attr_filtered <- metadata %>% filter(sample_id %in% sample_ids)
  
  num_samples <- ncol(normalized_counts)
  if (num_samples < 20) {
    cat('Skipping tissue:', tissue_name, '(Insufficient samples:', num_samples, ')\n')
    return(NULL)
  }
  
  # Call simplified academic batch correction
  adjusted_expression_data <- academic_batch_correction(normalized_counts, attr_filtered, tissue_name)
  
  # Convert back to data frame
  adjusted_expression_data <- as.data.frame(adjusted_expression_data)
  
  # Save adjusted data
  result_file <- file.path(processed_dir, 'expression/adjusted_sva_all',
                           paste0(tissue_name, '.rds'))
  saveRDS(adjusted_expression_data, file = result_file)
  
  cat('Processed tissue:', tissue_name, '\n')
  cat('Dimensions of adjusted data:', dim(adjusted_expression_data), '\n')
}

# -----------------------------
# Run pipeline for all tissues
# -----------------------------
for (tissue_file in tissue_files) {
  process_tissue(tissue_file, metadata)
}



