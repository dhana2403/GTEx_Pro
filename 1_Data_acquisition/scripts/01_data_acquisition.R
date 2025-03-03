################################Data acquisition#################################################

# Ensure the output directory paths are passed to the script
output_dir <- Sys.getenv("OUTPUT_DIR")  # Get the output directory path from the environment variable

library(tidyverse)

# Create directories within the output directory
raw_dir <- file.path(output_dir, "data", "raw")
processed_dir <- file.path(output_dir, "data", "processed")
metadata_dir <- file.path(output_dir, "data", "metadata")

dir.create(raw_dir, recursive = TRUE)
dir.create(processed_dir, recursive = TRUE)
dir.create(metadata_dir, recursive = TRUE)

# Check the directory structure
print(paste("Created directories: ", raw_dir, processed_dir, metadata_dir))

# Use wget to download files into the specified directories
system(paste('wget -nv "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx" -P', shQuote(metadata_dir)))
system(paste('wget -nv "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx" -P', shQuote(metadata_dir)))
system(paste('wget -nv "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt" -P', shQuote(metadata_dir)))
system(paste('wget -nv "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt" -P', shQuote(metadata_dir)))
system(paste('wget -nv "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz" -P', shQuote(raw_dir)))

# Unzip the downloaded files in the 'raw' directory
system(paste('cd', shQuote(raw_dir), '&& gunzip *'))

# Output to confirm successful execution
print("Download and extraction complete.")
