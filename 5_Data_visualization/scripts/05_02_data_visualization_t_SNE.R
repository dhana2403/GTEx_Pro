########################################### t-SNE VISUALIZATION BY TISSUES ########################################################

library(dplyr)
library(ggplot2)
library(edgeR)
library(Rtsne) 

data_path <- "./data/processed/expression/adjusted_sva"

# List all .rds files under the data path
tissue_files <- list.files(path = data_path, pattern = "*.rds", full.names = TRUE)

# Initialize a list to store combined data
combined_data <- list()

# Loop through each tissue file
for (tissue_file in tissue_files) {
  
  # Extract the tissue name
  tissue_name <- gsub(".rds$", "", basename(tissue_file))
  
  normalized_counts <- readRDS(tissue_file)
  
  # Check for empty normalized counts
  if (nrow(normalized_counts) == 0) {
    stop(paste("Error: The normalized counts for", tissue_name, "are empty."))
  }
  
  # Replace negative counts with a small positive value
  normalized_counts[normalized_counts < 0] <- 1e-9
  
  # Convert to data frame and transpose for t-SNE
  data_for_tsne <- as.data.frame(t(normalized_counts))
  
  data_for_tsne$tissue <- tissue_name
  
  # Append to the combined data list
  combined_data[[tissue_name]] <- data_for_tsne
}

# Combine all data into one data frame
combined_data_df <- do.call(rbind, combined_data)

# Run t-SNE
set.seed(42)  # Set seed for reproducibility
tsne_result <- Rtsne(as.matrix(combined_data_df[, -ncol(combined_data_df)]), 
                     dims = 2, 
                     perplexity = 30, 
                     verbose = TRUE, 
                     check_duplicates = FALSE)

# Create a data frame with the t-SNE results
tsne_data <- as.data.frame(tsne_result$Y)
tsne_data$tissue <- combined_data_df$tissue

tissue_colors <- c(
  "Brain-Cortex" = "black",
  "Heart-AtrialAppendage" = "red",
  "Heart-LeftVentricle" = "blue",
  "Lung" = "orange",
  "Muscle-Skeletal" = "cyan",
  "Skin-NotSunExposed(Suprapubic)" = "magenta",
  "Skin-SunExposed(Lowerleg)" = "brown"
)

# Create t-SNE plot for all tissues
tsne_plot <- ggplot(tsne_data, aes(x = V1, y = V2, color = tissue)) +
  geom_point(alpha = 0.7, size = 3) +
  scale_color_manual(values = tissue_colors) +
  labs(
    title = "t-SNE of Gene Expression Across Tissues",
    x = "t-SNE Dimension 1",
    y = "t-SNE Dimension 2",
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

print(tsne_plot)

result_directory <- "results_after_sva_tsne_all_tissues"
dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(result_directory, "tsne_plot_all_tissues_sva.pdf"), tsne_plot, units = 'cm', width = 18, height = 18, useDingbats = FALSE)
ggsave(filename = file.path(result_directory, "tsne_plot_all_tissues_sva.png"), plot = tsne_plot, units = "cm", width = 18, height = 18, dpi = 300)
