# Load necessary libraries
library(dplyr)
library(ggplot2)
library(edgeR)
library(cluster)
library(clusterSim)  # For Davies-Bouldin Index
library(ggpubr)  # Use ggpubr for better significance testing

# Define the data paths
tmm_data_path <- "./data/processed/expression/readcounts_tmm_all"
sva_data_path <- "./data/processed/expression/adjusted_sva_all"

# Function to process data, perform PCA, and calculate DBI
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
  
  # Assign clusters based on tissue names
  cluster_labels <- as.numeric(as.factor(pca_data$tissue))
  
  # Compute Davies-Bouldin Index
  dbi_value <- index.DB(pca_data[, 1:3], cluster_labels, centrotypes = "centroids")$DB
  
  # Return a data frame with DBI values for the dataset
  return(data.frame(Groups = basename(data_path), DBI = dbi_value))
}

# Process the datasets and compute DBI
tmm_dbi <- process_data(tmm_data_path)
sva_dbi <- process_data(sva_data_path)

# Combine the results
dbi_results <- bind_rows(
  mutate(tmm_dbi, Groups = "TMM+CPM"),
  mutate(sva_dbi, Groups = "TMM+CPM+SVA")
)

# Convert Groups to factor for plotting
dbi_results$Groups <- factor(dbi_results$Groups, levels = c("TMM+CPM", "TMM+CPM+SVA"))

# Generate the DBI comparison plot with significance testing
ggplot(dbi_results, aes(x = Groups, y = DBI, fill = Groups)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
  labs(
    x = "Processing Method",
    y = "Davies-Bouldin Index (DBI)"
  ) +
  scale_fill_manual(values = c("TMM+CPM" = "lightblue", "TMM+CPM+SVA" = "lightpink")) +
  theme_minimal() +
  theme(
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  ) 