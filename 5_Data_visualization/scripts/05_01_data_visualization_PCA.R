
########################################### PCA VISUALIZATION BY TISSUES ########################################################

# Load required packages
library(dplyr)
library(ggplot2)
library(edgeR)

# Define the path to the folder containing normalized count values
data_path <- "./data/processed/expression/adjusted_sva"

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
  
  # Add tissue name to the data frame
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

# Define colors for tissues
tissue_colors <- c(
  "Brain-Cortex" = "black",
  "Heart-AtrialAppendage" = "red",
  "Heart-LeftVentricle" = "blue",
  "Lung" = "orange",
  "Muscle-Skeletal" = "cyan",
  "Skin-NotSunExposed" = "magenta",
  "Skin-SunExposed" = "brown"
)

# Create PCA plot for all tissues
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = tissue)) +
  geom_point(alpha = 0.7, size = 3) +
  scale_color_manual(values = tissue_colors) +
  labs(
    title = paste("PCA of Gene Expression Across Tissues\nPC1:", pc1_var, "%, PC2:", pc2_var, "%"),
    x = paste("Principal Component 1 (", pc1_var, "%)", sep = ""),
    y = paste("Principal Component 2 (", pc2_var, "%)", sep = ""),
    color = "Tissue"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Print PCA plot to check if it renders
print(pca_plot)

# Create a directory to save the plots
result_directory <- "results_after_sva_pca_all_tissues"
dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)

# Save the PCA plot for all tissues
ggsave(file.path(result_directory, "pca_plot_all_tissues_sva.pdf"), pca_plot, units = 'cm', width = 18, height = 18, useDingbats = FALSE)
ggsave(filename = file.path(result_directory, "pca_plot_all_tissues_sva.png"), plot = pca_plot, units = "cm", width = 18, height = 18, dpi = 300)
