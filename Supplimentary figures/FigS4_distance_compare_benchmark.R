# Load necessary libraries
library(dplyr)
library(ggplot2)
library(edgeR)
library(cluster)
library(tibble)
library(reshape2)
library(tidyr)

# Define the data paths
tpm_data_path <- "./data/processed/expression/tpm_all"
tpm_sva_data_path <- "./data/processed/expression/adjusted_sva_tpm_all"
tmm_data_path <- "./data/processed/expression/readcounts_tmm_all"
tmm_sva_data_path <- "./data/processed/expression/adjusted_sva_all"

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

# Process the datasets (TPM, TPM+SVA, TMM, TMM+SVA)
tpm_distances <- process_data(tpm_data_path)
tpm_sva_distances <- process_data(tpm_sva_data_path)
tmm_distances <- process_data(tmm_data_path)
tmm_sva_data_distances <- process_data(tmm_sva_data_path)

# Function to combine the distance matrices for plotting
combine_distances_for_plotting <- function(tmm_distances, tmm_sva_distances, tpm_distances, tpm_sva_distances) {
  # Convert distance matrices into long format using melt
  tmm_melt <- melt(tmm_distances)
  tmm_melt$dataset <- "TMM+CPM"
  
  tmm_sva_melt <- melt(tmm_sva_distances)
  tmm_sva_melt$dataset <- "TMM+CPM+SVA"
  
  tpm_melt <- melt(tpm_distances)
  tpm_melt$dataset <- "TPM"
  
  tpm_sva_melt <- melt(tpm_sva_distances)
  tpm_sva_melt$dataset <- "TPM+SVA"
  
  # Combine all datasets into one data frame
  combined_distances <- bind_rows(tmm_melt, tmm_sva_melt, tpm_melt, tpm_sva_melt)
  
  # Rename columns for clarity
  colnames(combined_distances) <- c("Tissue1", "Tissue2", "Distance", "Dataset")
  
  return(combined_distances)
}

# Combine the distances for the datasets
combined_distances <- combine_distances_for_plotting(tmm_distances, tmm_sva_data_distances, tpm_distances, tpm_sva_distances)

# Calculate the average distance for each tissue across all datasets
average_distances <- combined_distances %>%
  group_by(Tissue1, Dataset) %>%
  summarise(avg_distance = mean(Distance, na.rm = TRUE)) %>%
  spread(Dataset, avg_distance)  # Spread to have each dataset in its own column

print(colnames(average_distances))

# Remove the 'Raw' dataset (as per your request)
average_distances <- average_distances %>%
  select(Tissue1, `TMM+CPM`, `TMM+CPM+SVA`, `TPM`, `TPM+SVA`)

# Convert data into long format for ggplot
long_distances <- average_distances %>%
  pivot_longer(cols = c(`TMM+CPM`, `TMM+CPM+SVA`, `TPM`, `TPM+SVA`), names_to = "Groups", values_to = "Average_Distance")

# Rename dataset labels for legend clarity
long_distances$Groups <- factor(long_distances$Groups,
                                levels = c("TMM+CPM", "TMM+CPM+SVA", "TPM", "TPM+SVA"),
                                labels = c("TMM+CPM", "TMM+CPM+SVA", "TPM", "TPM+SVA"))

# Visualize the average distance between tissues for each dataset using a bar plot
bar_plot <- ggplot(long_distances, aes(x = Tissue1, y = Average_Distance, fill = Groups)) +
  geom_bar(stat = "identity", position = "dodge") +  # Dodge to separate bars
  scale_fill_manual(values = c("lightblue", "pink", "lightgreen","purple")) +  # Custom colors
  labs(x = "Tissue", y = "Average Distance Between Tissues", fill = "Dataset") +  # Legend title
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.background = element_rect(fill = "white", color = NA),  # White background, no borders
    plot.background = element_rect(fill = "white", color = NA),  # White plot background
    panel.grid = element_blank(),  # Remove gridlines
    axis.line = element_blank(),  # Remove axis lines
    legend.position = "top"
  )

# Save the plot as a PNG image
ggsave("average_distance_between_tissues_benchmark.png", plot = bar_plot, width = 12, height = 8, dpi = 300)

# Save the plot as a PDF
ggsave("average_distance_between_tissues_benchmark.pdf", plot = bar_plot, width = 12, height = 8)
