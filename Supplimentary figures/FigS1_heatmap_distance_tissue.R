# Load necessary libraries
library(dplyr)
library(ggplot2)
library(edgeR)
library(cluster)
library(tibble)
library(reshape2)

# Define the data paths
tmm_data_path <- "./data/processed/expression/readcounts_tmm_all"
sva_data_path <- "./data/processed/expression/adjusted_sva_all"

# Function to process data, perform PCA, and calculate cluster distances
process_data <- function(data_path) {
  tissue_files <- list.files(path = data_path, pattern = "*.rds", full.names = TRUE)
  combined_data <- list()
  
  # Loop through each tissue file
  for (tissue_file in tissue_files) {
    tissue_name <- gsub(".rds$", "", basename(tissue_file))
    normalized_counts <- readRDS(tissue_file)
    
    # Check for empty normalized counts
    if (nrow(normalized_counts) == 0) {
      stop(paste("Error: The normalized counts for", tissue_name, "are empty."))
    }
    
    # Replace negative counts with a small positive value
    normalized_counts[normalized_counts < 0] <- 1e-9
    
    # Prepare data for PCA (transpose the data for samples to rows)
    data_for_pca <- as.data.frame(t(normalized_counts))
    data_for_pca$tissue <- tissue_name
    
    # Append to the combined data list
    combined_data[[tissue_name]] <- data_for_pca
  }
  
  # Combine all data into one data frame
  combined_data_df <- do.call(rbind, combined_data)
  
  # Perform PCA on the combined data
  pca_result <- prcomp(combined_data_df[, -ncol(combined_data_df)], center = TRUE, scale. = TRUE)
  
  # Extract PCA results for samples
  pca_data <- as.data.frame(pca_result$x)
  pca_data$tissue <- combined_data_df$tissue  # Adding tissue information
  
  # Compute cluster centroids
  cluster_centroids <- pca_data %>%
    group_by(tissue) %>%
    summarise(PC1 = mean(PC1), PC2 = mean(PC2), PC3 = mean(PC3))
  
  # Convert to a data frame (instead of tibble) to set rownames properly
  cluster_centroids <- as.data.frame(cluster_centroids)
  
  # Convert the tissue column to rownames
  rownames(cluster_centroids) <- cluster_centroids$tissue
  cluster_centroids$tissue <- NULL  # Remove the tissue column from the data
  
  # Compute pairwise Euclidean distances between clusters
  dist_matrix <- as.matrix(dist(cluster_centroids))
  
  # Set row and column names of the distance matrix to tissue names
  rownames(dist_matrix) <- colnames(dist_matrix) <- rownames(cluster_centroids)
  
  return(dist_matrix)
}

# Function to combine the distance matrices for plotting
combine_distances_for_plotting <- function(tmm_distances, sva_distances) {
  # Convert distance matrices into long format using melt
  tmm_melt <- melt(tmm_distances)
  tmm_melt$dataset <- "TMM"
  
  sva_melt <- melt(sva_distances)
  sva_melt$dataset <- "SVA"
  
  # Combine all datasets into one data frame
  combined_distances <- bind_rows(tmm_melt, sva_melt)
  
  # Rename columns for clarity
  colnames(combined_distances) <- c("Tissue1", "Tissue2", "Distance", "Dataset")
  
  return(combined_distances)
}

# Process the datasets (TMM+CPM, TMM+CPM+SVA)
tmm_distances <- process_data(tmm_data_path)
sva_distances <- process_data(sva_data_path)

# Combine the distances for the two datasets
combined_distances <- combine_distances_for_plotting(tmm_distances, sva_distances)

# Calculate global min and max distance values for consistent color scale
global_min_distance <- min(combined_distances$Distance)
global_max_distance <- max(combined_distances$Distance)

# Create a function to save plots for each dataset
save_heatmap_for_each_dataset <- function(dataset_name, data) {
  heatmap_plot <- ggplot(data, aes(x = Tissue1, y = Tissue2, fill = Distance)) + 
    geom_tile() +  # Create the heatmap
    geom_text(aes(label = round(Distance, 2)), color = "black", size = 3) +  # Add data values inside the tiles
    scale_fill_gradientn(
      colours = c("white", "lightblue", "blue", "darkblue"),  # Use a fine gradient with more shades
      limits = c(global_min_distance, global_max_distance),  # Set global limits to ensure consistent scale
      oob = scales::squish  # Ensure extreme values are mapped within the color scale
    ) +  # Color gradient for distances
    theme_minimal() +  # Minimal theme
    theme(axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate axis labels for clarity
          axis.text.y = element_text(angle = 0, hjust = 1),  # Rotate y-axis labels
          strip.text = element_text(size = 12),  # Adjust facet title size
          axis.title = element_text(size = 14),  # Axis title size
          axis.text = element_text(size = 10)) +  # Axis text size
    labs(title = paste("Tissue Distance Comparison for", dataset_name), x = "Tissue", y = "Tissue")
  
  # Save the heatmap plot for the dataset
  ggsave(paste0(dataset_name, "_tissue_distance_comparison_heatmap.png"), plot = heatmap_plot, width = 25, height = 23, dpi = 300)
  ggsave(paste0(dataset_name, "_tissue_distance_comparison_heatmap.pdf"), plot = heatmap_plot, width = 25, height = 23)
}

# Filter data for each dataset and save separate plots
tmm_data <- filter(combined_distances, Dataset == "TMM")
sva_data <- filter(combined_distances, Dataset == "SVA")

# Save plots for each dataset
save_heatmap_for_each_dataset("TMM", tmm_data)
save_heatmap_for_each_dataset("SVA", sva_data)
