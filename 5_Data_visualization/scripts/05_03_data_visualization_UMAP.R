########################################### UMAP VISUALIZATION BY TISSUES ########################################################

# Load required packages
library(dplyr)
library(ggplot2)
library(edgeR)
library(umap)  # Load the umap package for UMAP

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
  
  # Convert to data frame and transpose for UMAP
  data_for_umap <- as.data.frame(t(normalized_counts))
  
  # Add tissue name to the data frame
  data_for_umap$tissue <- tissue_name
  
  # Append to the combined data list
  combined_data[[tissue_name]] <- data_for_umap
}

# Combine all data into one data frame
combined_data_df <- do.call(rbind, combined_data)

# Run UMAP
set.seed(42)  # Set seed for reproducibility
umap_result <- umap(as.matrix(combined_data_df[, -ncol(combined_data_df)]), 
                    n_neighbors = 15, 
                    min_dist = 0.1, 
                    metric = "euclidean")

# Create a data frame with the UMAP results
umap_data <- as.data.frame(umap_result$layout)  # Use the $layout attribute
umap_data$tissue <- combined_data_df$tissue

# Define colors for tissues
tissue_colors <- c(
  "Brain-Cortex" = "black",
  "Heart-AtrialAppendage" = "red",
  "Heart-LeftVentricle" = "blue",
  "Lung" = "orange",
  "Muscle-Skeletal" = "cyan",
  "Skin-NotSunExposed(Suprapubic)" = "magenta",
  "Skin-SunExposed(Lowerleg)" = "brown"
)

# Create UMAP plot for all tissues
umap_plot <- ggplot(umap_data, aes(x = V1, y = V2, color = tissue)) +
  geom_point(alpha = 0.7, size = 3) +
  scale_color_manual(values = tissue_colors) +
  labs(
    title = "UMAP of Gene Expression Across Tissues",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2",
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

# Print UMAP plot to check if it renders
print(umap_plot)

# Create a directory to save the plots
result_directory <- "results_after_sva_umap_all_tissues"
dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)

# Save the UMAP plot for all tissues
ggsave(file.path(result_directory, "umap_plot_all_tissues_sva.pdf"), umap_plot, units = 'cm', width = 18, height = 18, useDingbats = FALSE)
ggsave(filename = file.path(result_directory, "umap_plot_all_tissues_sva.png"), plot = umap_plot, units = "cm", width = 18, height = 18, dpi = 300)
