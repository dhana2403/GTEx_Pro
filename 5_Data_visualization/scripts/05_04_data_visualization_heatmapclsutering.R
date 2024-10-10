
################################################ VISUALIZATION BY HEATMAPCLUSTERING #####################################

# Load necessary libraries
library(ComplexHeatmap)
library(circlize)  # For colorRamp2 function

# Define the path to the folder containing normalized count values
data_path <- "./adjusted_sva"

# List all .rds files under the data path
tissue_files <- list.files(path = data_path, pattern = "*.rds", full.names = TRUE)

# Extract tissue names from the file names
tissue_names <- gsub(".rds$", "", basename(tissue_files))

# Define gene names manually
 gene_names <- # Input your geneid of interest

# Loop through each tissue file
for (tissue_file in tissue_files) {
  
  # Extract the tissue name
  tissue_name <- gsub(".rds$", "", basename(tissue_file))
  
  # Load the normalized read counts data from the .rds file
  normalized_counts <- readRDS(tissue_file)
  
  # Ensure the number of gene names matches the number of rows
  if (nrow(normalized_counts) != length(gene_names)) {
    stop(paste("Error: The number of gene names does not match the number of rows in", tissue_name))
  }
  
  # Replace rownames with the manually defined gene names
  rownames(normalized_counts) <- gene_names
  
  # Convert to data frame and log-transform
  normalized_counts <- as.data.frame(normalized_counts)
  data_log10 <- log10(normalized_counts + 1e-9)  # Add a small constant to avoid log(0)
  
  # Z-score normalization
  data_zscore <- as.data.frame(t(scale(t(data_log10))))
  
  # Check for any NA values and remove rows/columns that are fully NA
  data_zscore <- na.omit(data_zscore)  # Removes any rows with NA values
  data_zscore <- data_zscore[, colSums(is.na(data_zscore)) == 0]  # Removes columns with NA values
  
  data_zscore_reordered <- as.matrix(data_zscore)

  # Create a directory to save the plots for the current tissue
  result_directory <- file.path("results_after_sva_heatmap_cluster", tissue_name)
  dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)
  
  # Generate the vertical heatmap with adjusted parameters
  p <- Heatmap(
    data_zscore,
    name = "Expression",
    cluster_rows = TRUE,           # Cluster rows (genes)
    cluster_columns = FALSE,       # Do not cluster columns to preserve order
    show_column_names = TRUE,      # Show column names
    show_row_names = TRUE,         # Show row names
    color = colorRampPalette(c("blue", "white", "red"))(50)  # Color gradient
  )
  
  # Save the heatmap as PDF and PNG with larger dimensions
  pdf(file.path(result_directory, paste0('heatmap_cluster_', tissue_name, '_vertical.pdf')), width = 20, height = 16)
  draw(p)  # Use draw() to render the heatmap
  dev.off()
  
  png(file.path(result_directory, paste0('heatmap_cluster_', tissue_name, '_vertical.png')), width = 6000, height = 8000, res = 300)
  draw(p)  # Use draw() to render the heatmap
  dev.off()
}

# Clean up
rm(list = ls())
