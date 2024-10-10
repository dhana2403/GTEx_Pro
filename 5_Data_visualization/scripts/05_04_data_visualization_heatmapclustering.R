
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

# Define group labels for each tissue
group_labels_list <- list(
  "Brain-Cortex" = c(rep("20-49", 9), rep("50-59", 18), rep("60-79", 55)),
  "Heart-AtrialAppendage" = c(rep("20-49", 9), rep("50-59", 20), rep("60-79", 49)), 
  "Heart-LeftVentricle" = c(rep("20-49", 10), rep("50-59", 14), rep("60-79", 38)),
  "Lung" = c(rep("20-49", 12), rep("50-59", 24), rep("60-79", 59)),
  "Muscle-Skeletal" = c(rep("20-49", 16), rep("50-59", 28), rep("60-79", 89)),
  "Skin-NotSunExposed" = c(rep("20-49", 11), rep("50-59", 20), rep("60-79", 76)),
  "Skin-SunExposed" = c(rep("20-49", 15), rep("50-59", 26), rep("60-79", 79))
)

# Define gene names manually
gene_names <- # Input your geneid of interest

# Loop through each tissue file
for (tissue_file in tissue_files) {
  
  # Extract the tissue name
  tissue_name <- gsub(".rds$", "", basename(tissue_file))
  
  # Check if the tissue is in the group labels list
  if (!tissue_name %in% names(group_labels_list)) {
    next  # Skip this tissue if group labels are not defined
  }
  
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
  
  # Define the group labels for the current tissue
  group_labels <- group_labels_list[[tissue_name]]
  
  # Ensure the length of group_labels matches the number of columns in data_zscore
  if (length(group_labels) != ncol(data_zscore)) {
    stop(paste("Error: The length of group_labels does not match the number of columns in", tissue_name))
  }
  
  # Convert data_zscore to a data frame
  data_zscore_transposed <- as.data.frame(t(data_zscore))
  
  # Convert 'group_labels' to a factor with explicit levels
  group_labels_factor <- factor(group_labels, levels = c("20-49", "50-59", "60-79"))
  
  # Add group labels as a column to the data matrix
  data_with_labels <- cbind(data_zscore_transposed, Group = group_labels_factor)
  
  # Ensure data_with_labels is a data frame
  if (!is.data.frame(data_with_labels)) {
    stop("Error: data_with_labels is not a data frame.")
  }
  
  # Reorder rows by group labels
  ordered_data <- data_with_labels[order(data_with_labels$Group), ]
  data_zscore_reordered <- ordered_data[, -ncol(ordered_data)]  # Remove the Group column for heatmap
  group_labels_reordered <- ordered_data$Group
  
  # Create the annotation data frame for columns
  annotation_col <- data.frame(AgeGroup = group_labels_reordered)
  
  # Define colors for annotations
  annotation_colors <- list(
    AgeGroup = c(
      "20-49" = "#87CEEB", 
      "50-59" = "#388E3C", 
      "60-79" = "#FFC107"
    )
  )
  
  # Create the heatmap annotation for columns
  ha <- HeatmapAnnotation(
    AgeGroup = annotation_col$AgeGroup,
    col = annotation_colors
  )
  
  data_zscore_reordered <- as.matrix(data_zscore_reordered)
  data_zscore_reordered <- t(data_zscore_reordered)

  # Create a directory to save the plots for the current tissue
  result_directory <- file.path("results_after_sva_cluster", tissue_name)
  dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)
  
  # Generate the vertical heatmap with adjusted parameters
  p <- Heatmap(
    data_zscore_reordered,
    name = "Expression",
    top_annotation = ha,           # Apply column annotation on the top
    cluster_rows = TRUE,           # Cluster rows (genes)
    cluster_columns = FALSE,       # Do not cluster columns to preserve order
    show_column_names = TRUE,      # Show column names
    show_row_names = TRUE,         # Show row names
    color = colorRampPalette(c("blue", "white", "red"))(50)  # Color gradient
  )
  
  # Save the heatmap as PDF and PNG with larger dimensions
  pdf(file.path(result_directory, paste0('heatmap_cluster_', tissue_name, '_vertical.pdf')), width = 20, height = 16)
  draw(p)  
  dev.off()
  
  png(file.path(result_directory, paste0('heatmap_cluster_', tissue_name, '_vertical.png')), width = 6000, height = 8000, res = 300)
  draw(p)  
  dev.off()
}

# Clean up
rm(list = ls())
