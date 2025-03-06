# Load necessary libraries
library(dplyr)
library(ggplot2)
library(edgeR)
library(plotly)
library(htmlwidgets)  # Load htmlwidgets for saveWidget function

# Define the path to the folder containing normalized count values
data_path <- "./data/processed/expression/readcounts_tmm_all_cross"

# List all .rds files under the data path
tissue_files <- list.files(path = data_path, pattern = "*.rds", full.names = TRUE)

# Initialize a list to store combined data
combined_data <- list()

# Loop through each tissue file
for (tissue_file in tissue_files) {
  
  # Extract the tissue name
  tissue_name <- gsub(".rds$", "", basename(tissue_file))
  
  # Load the normalized read counts data from the .rds file
  normalized_counts <- readRDS(tissue_file)
  
  # Check for empty normalized counts
  if (nrow(normalized_counts) == 0) {
    stop(paste("Error: The normalized counts for", tissue_name, "are empty."))
  }
  
  # Replace negative counts with a small positive value
  normalized_counts[normalized_counts < 0] <- 1e-9
  
  # Convert to data frame and transpose for PCA
  data_for_pca <- as.data.frame(t(normalized_counts))
  data_for_pca$tissue <- tissue_name
  
  # Append to the combined data list
  combined_data[[tissue_name]] <- data_for_pca
}

# Combine all data into one data frame
combined_data_df <- do.call(rbind, combined_data)

# Perform PCA on combined data
pca_result <- prcomp(combined_data_df[, -ncol(combined_data_df)], center = TRUE, scale. = TRUE)

# Extract PCA results for samples
pca_data <- as.data.frame(pca_result$x)
pca_data$tissue <- combined_data_df$tissue

# Calculate percentage of variance explained
pca_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
pc1_var <- round(pca_variance[1], 2)
pc2_var <- round(pca_variance[2], 2)
pc3_var <- round(pca_variance[3], 2)  # For 3D, we will also use PC3

# Define specific colors for certain tissues
specific_tissue_colors <- c(
  "Brain-Cortex" = "black",
  "Heart-AtrialAppendage" = "red",
  "Heart-LeftVentricle" = "blue",
  "Liver" = "cyan",
  "Skin-NotSunExposed" = "pink",
  "Skin-SunExposed" = "maroon"
)

# Identify the remaining tissues that are not in the predefined list
remaining_tissues <- setdiff(unique(combined_data_df$tissue), names(specific_tissue_colors))

# Generate random colors for the remaining tissues
set.seed(42)  # Set seed for reproducibility
random_colors <- colors()[sample(length(colors()), length(remaining_tissues))]

# Assign the random colors to the remaining tissues
random_tissue_colors <- setNames(random_colors, remaining_tissues)

# Merge specific and random tissue colors
tissue_colors <- c(specific_tissue_colors, random_tissue_colors)

# Create 3D PCA plot using plotly
pca_3d_plot <- plot_ly(
  data = pca_data, 
  x = ~PC1, 
  y = ~PC2, 
  z = ~PC3, 
  color = ~tissue, 
  colors = tissue_colors, 
  type = 'scatter3d', 
  mode = 'markers', 
  marker = list(size = 5, opacity = 0.7)
) %>%
  layout(
    title = paste("3D PCA of Gene Expression Across Tissues\nPC1:", pc1_var, "%, PC2:", pc2_var, "%, PC3:", pc3_var, "%"),
    scene = list(
      xaxis = list(title = paste("Principal Component 1 (", pc1_var, "%)", sep = "")),
      yaxis = list(title = paste("Principal Component 2 (", pc2_var, "%)", sep = "")),
      zaxis = list(title = paste("Principal Component 3 (", pc3_var, "%)", sep = ""))
    ),
    legend = list(title = list(text = 'Tissue'))
  )

# Print the 3D PCA plot
pca_3d_plot

# Save the plot as an HTML file for interactive viewing
result_directory <- "results_after_readcounts_pca_all_tissues"
dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)

# Save as an interactive HTML file
html_file <- file.path(result_directory, "3d_pca_plot_all_tissues_tmm_cross.html")
saveWidget(pca_3d_plot, html_file)
