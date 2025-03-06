library(dplyr)
library(ggplot2)
library(tidyr)

# Define the paths to the two datasets
data_path_1 <- "./data/processed/expression/adjusted_sva_tpm_all_tis"  # TMM+CPM
data_path_2 <- "./data/processed/expression/adjusted_sva_all_tis"   # TMM+CPM+SVA

# Define the gene ID of interest
gene_id_of_interest <- "ENSG00000132639"

# Function to extract gene expression data for a given dataset
extract_gene_expression <- function(data_path, gene_id_of_interest, dataset_name) {
  # List all .rds files in the data path
  tissue_files <- list.files(path = data_path, pattern = "*.rds", full.names = TRUE)
  
  combined_data <- list()  # List to store data for all tissues
  
  for (tissue_file in tissue_files) {
    # Extract the tissue name
    tissue_name <- gsub(".rds$", "", basename(tissue_file))
    
    # Load the normalized read counts data from the .rds file
    normalized_counts <- readRDS(tissue_file)
    
    # Check if the gene of interest is in the data
    if (gene_id_of_interest %in% rownames(normalized_counts)) {
      
      # Extract expression values for the gene of interest
      gene_expression <- normalized_counts[gene_id_of_interest, , drop = FALSE]
      
      # Transform the data (log10 transformation with a small constant)
      small_constant <- 1e-9
      gene_expression <- log10(gene_expression + small_constant)
      
      # Convert to long format
      df_long <- data.frame(
        Sample = colnames(gene_expression),
        Expression = as.numeric(gene_expression),
        Dataset = dataset_name,  # Assign the dataset name
        Tissue = tissue_name
      )
      
      # Append the data to the combined list
      combined_data[[tissue_name]] <- df_long
    }
  }
  
  # Combine all tissues' data into a single data frame
  combined_df <- bind_rows(combined_data)
  return(combined_df)
}

# Extract gene expression data from both datasets with specific dataset names
expression_data_1 <- extract_gene_expression(data_path_1, gene_id_of_interest, "TMM+CPM")
expression_data_2 <- extract_gene_expression(data_path_2, gene_id_of_interest, "TMM+CPM+SVA")

# Combine the data from both datasets
expression_data <- bind_rows(expression_data_1, expression_data_2)

# Plot the expression trajectory of the gene across all tissues and datasets
p <- ggplot(expression_data, aes(x = Tissue, y = Expression, group = Dataset, color = Dataset)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  labs(
    title = paste("Expression Trajectory of", gene_id_of_interest),
    x = "Tissue", 
    y = "Log10(Expression + Small Constant)"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.line = element_line(color = "black"),  # Add x and y-axis lines
    plot.title = element_text(size = 10, hjust = 0.5)
  )

# Save the plot output
output_dir <- "results_gene_trajectory_tis"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(output_dir, paste0("trajectory_plot_", gene_id_of_interest, "_all_tissues.pdf")), p, units = "cm", width = 24, height = 12)
ggsave(file.path(output_dir, paste0("trajectory_plot_", gene_id_of_interest, "_all_tissues.png")), p, units = "cm", width = 24, height = 12, dpi = 300)
