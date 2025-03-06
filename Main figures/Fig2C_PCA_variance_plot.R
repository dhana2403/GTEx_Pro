library(dplyr)
library(ggplot2)
library(edgeR)

process_pca_data <- function(data_path, label) {
  # List all .rds files under the data path
  tissue_files <- list.files(path = data_path, pattern = "*.rds", full.names = TRUE)
  
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
    
    data_for_pca <- as.data.frame(t(normalized_counts))
    data_for_pca$tissue <- tissue_name
    
    # Append to the combined data list
    combined_data[[tissue_name]] <- data_for_pca
  }
  
  # Combine all data into one data frame
  combined_data_df <- do.call(rbind, combined_data)
  
  # Perform PCA on combined data
  pca_result <- prcomp(combined_data_df[, -ncol(combined_data_df)], center = TRUE, scale. = TRUE)
  
  # Calculate percentage of variance explained
  pca_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
  
  # Create a data frame for variance explained by each PC
  variance_data <- data.frame(
    PC = factor(paste0("PC", seq_along(pca_variance)), levels = paste0("PC", seq_along(pca_variance))),
    Variance = pca_variance,
    Groups = label
  )
  
  return(variance_data)
}

# Process both datasets
readcounts_tmm_path <- "./data/processed/expression/readcounts_tmm_all"
adjusted_sva_path <- "./data/processed/expression/adjusted_sva_all"


readcounts_tmm_data <- process_pca_data(readcounts_tmm_path, "TMM+CPM")
adjusted_sva_data <- process_pca_data(adjusted_sva_path, "TMM+CPM+SVA")


# Combine data for plotting
combined_variance_data <- rbind(adjusted_sva_data, readcounts_tmm_data)

# Filter the combined_variance_data to only include PC1 to PC10
combined_variance_data_filtered <- combined_variance_data %>%
  filter(PC %in% paste0("PC", 1:10))

# Grouped bar plot of variance explained with increased spacing between PCs, showing only PC1 to PC10
variance_plot <- ggplot(combined_variance_data_filtered, aes(x = PC, y = Variance, fill = Groups)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +  # Increased width for more space
  geom_text(aes(label = paste0(round(Variance, 2), "%")),
            position = position_dodge(width = 1.2), vjust = -1, size = 3) + 
  labs(
    x = "Principal Component",
    y = "Percentage of Variance Explained"
  ) + 
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) + 
  scale_x_discrete(expand = c(0, 0)) +  # Removes extra space before PC1
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.01)),  # Keep a little space at the bottom
    limits = c(0, 40)  # Adjust the y-axis to go from 0 to 50 (or any value you prefer)
  ) + 
  scale_fill_manual(values = c("TMM+CPM" = "grey", "TMM+CPM+SVA" = "lightblue"))

# Save the filtered variance plot
result_directory <- "results_after_readcounts_pca_comparison"
dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(result_directory, "combined_variance_explained_comparison_PC1_to_PC10_barplot.pdf"), variance_plot, units = 'cm', width = 25, height = 15, useDingbats = FALSE)  # Reduced height here
ggsave(filename = file.path(result_directory, "combined_variance_explained_comparison_PC1_to_PC10_barplot.png"), plot = variance_plot, units = "cm", width = 25, height = 15, dpi = 300)  # Reduced height here
