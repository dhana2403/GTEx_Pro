##################################################
# Data Download and Directory Creation (First part)
##################################################

# Ensure the output directory paths are passed to the script
output_dir <- Sys.getenv("OUTPUT_DIR")  # Get the output directory path from the environment variable

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

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
system(paste('wget -nv "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt" -P', shQuote(metadata_dir)))
system(paste('wget -nv "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt" -P', shQuote(metadata_dir)))
system(paste('wget -nv "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx" -P', shQuote(metadata_dir)))
system(paste('wget -nv "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx" -P', shQuote(metadata_dir)))
system(paste('wget -nv "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz" -P', shQuote(raw_dir)))

# Unzip the downloaded files in the 'raw' directory
system(paste('cd', shQuote(raw_dir), '&& gunzip *'))

# Output to confirm successful execution
print("Download and extraction complete.")

##################################################
# Data Preprocessing (Second part)
##################################################

args <- commandArgs(trailingOnly = TRUE)
print(paste("Raw arguments from Nextflow:", paste(args, collapse = " ")))

if (length(args) == 0) {
  stop("Error: No genes_of_interest were provided. Check Nextflow input.")
}

genes_of_interest <- as.character(args)
print(paste("Genes of interest:", paste(genes_of_interest, collapse = ", ")))

# Load required libraries
library(data.table)
library(tidyverse)
library(dplyr)

# Load and process metadata files
att <- read_tsv(file.path(metadata_dir, 'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')) %>%
  mutate(SUBJID = sapply(strsplit(SAMPID, '-'), function(x) paste(x[1], x[2], sep = '-'))) %>%
  select(SAMPID, SMTS, SMTSD, SMNABTCH, SMGEBTCH, SUBJID, SMRIN, SMTSISCH) %>%
  set_names(c('sample_id', 'major_tissue', 'minor_tissue', 'batch1', 'batch2', 'subj_id', 'rin', 'ischemic_time')) %>%
  unique()

phe <- read_tsv(file.path(metadata_dir, 'GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')) %>%
  set_names(c('subj_id', 'sex', 'age', 'death')) %>%
  mutate(
    sex = as.factor(c('male', 'female')[sex]),
    age = as.factor(age),
    death = factor(death, levels = 0:4,
                   labels = c('ventilator', 'fastdeath_violent', 'fastdeath_naturalcause', 'intermediatedeath', 'slowdeath'))
  )

# Merge and filter metadata
attphe <- att %>%
  full_join(phe, by = "subj_id") %>%
  unique() %>%
  drop_na()  # Omit missing values

# Filtered metadata
attphe_filtered <- attphe %>%
  filter(!minor_tissue %in% c("Cervix - Ectocervix", "Cervix - Endocervix", "Fallopian Tube",
                              "Testis", "Uterus", "Vagina", "Ovary", "Prostate","Breast - Mammary Tissue"))

# Load and process gene expression data
dat <- fread(file.path(raw_dir, 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct'), select = c('Name', attphe_filtered$sample_id)) %>%
  as.data.frame()

# Remove version numbers from gene names
rownames(dat) <- gsub("\\.\\d+$", "", dat$Name)
dat$Name <- NULL

# Convert to matrix and filter based on genes of interest
dat_matrix <- as.matrix(dat)
dat_matrix[!is.finite(dat_matrix)] <- NA # Convert infinite values to NA
dat_matrix <- na.omit(dat_matrix) # Remove rows with NA values

genes_of_interest <- gsub("\\.\\d+$", "", genes_of_interest)

# Filter the data for genes of interest
dat_filtered <- dat_matrix[rownames(dat_matrix) %in% genes_of_interest, ]

# Ensure samples match metadata and filter accordingly
samplesx <- intersect(attphe_filtered$sample_id, colnames(dat_filtered))
dat_filtered <- dat_filtered[, samplesx]

# Remove columns with zero variance
dat_filtered <- dat_filtered[, apply(dat_filtered, 2, var) != 0]

# Group samples by tissue based on filtered metadata
samples_by_tissues <- tapply(attphe_filtered$sample_id, INDEX = attphe_filtered$minor_tissue, FUN = unique)
dat_filtered_tissues <- lapply(samples_by_tissues, function(samps) {
  samps <- intersect(samps, samplesx)
  if (length(samps) > 1) { # Ensure each tissue group has at least 2 samples
    return(dat_filtered[, samps, drop = FALSE])
  } else {
    return(NULL) # Skip tissue groups with insufficient samples
  }
})

# Remove NULL elements (tissues with too few samples)
dat_filtered_tissues <- Filter(Negate(is.null), dat_filtered_tissues)

# Clean Row names
names(dat_filtered_tissues) <- sapply(strsplit(gsub(' ', '', names(dat_filtered_tissues)), '[()]'), function(x) x[[1]])

# Create directory if it does not exist
dir.create(file.path(processed_dir, 'expression', 'readcounts_all'), recursive = TRUE, showWarnings = FALSE)

# Save the filtered data and metadata
sapply(names(dat_filtered_tissues), function(nm) {
  saveRDS(dat_filtered_tissues[[nm]], file.path(processed_dir, 'expression', 'readcounts_all', paste0(nm, '.rds')))
})

# Save the metadata
saveRDS(attphe_filtered, file.path(processed_dir, 'attphe_all.rds'))

print("Data preprocessing and saving complete.")
##################################################
# Data Download and Directory Creation (First part)
##################################################
#' @export

# Ensure the output directory paths are passed to the script
output_dir <- Sys.getenv("OUTPUT_DIR")  # Get the output directory path from the environment variable

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

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
system(paste('wget -nv "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt" -P', shQuote(metadata_dir)))
system(paste('wget -nv "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt" -P', shQuote(metadata_dir)))
system(paste('wget -nv "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx" -P', shQuote(metadata_dir)))
system(paste('wget -nv "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx" -P', shQuote(metadata_dir)))
system(paste('wget -nv "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz" -P', shQuote(raw_dir)))

# Unzip the downloaded files in the 'raw' directory
system(paste('cd', shQuote(raw_dir), '&& gunzip *'))

# Output to confirm successful execution
print("Download and extraction complete.")


##################################################
# Data Preprocessing (Second part)
##################################################

#' @export
args <- commandArgs(trailingOnly = TRUE)
print(paste("Raw arguments from Nextflow:", paste(args, collapse = " ")))

if (length(args) == 0) {
  stop("Error: No genes_of_interest were provided. Check Nextflow input.")
}

genes_of_interest <- as.character(args)
print(paste("Genes of interest:", paste(genes_of_interest, collapse = ", ")))

# Load required libraries
library(data.table)
library(tidyverse)
library(dplyr)

# Load and process metadata files
att <- read_tsv(file.path(metadata_dir, 'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')) %>%
  mutate(SUBJID = sapply(strsplit(SAMPID, '-'), function(x) paste(x[1], x[2], sep = '-'))) %>%
  select(SAMPID, SMTS, SMTSD, SMNABTCH, SMGEBTCH, SUBJID, SMRIN, SMTSISCH) %>%
  set_names(c('sample_id', 'major_tissue', 'minor_tissue', 'batch1', 'batch2', 'subj_id', 'rin', 'ischemic_time')) %>%
  unique()

phe <- read_tsv(file.path(metadata_dir, 'GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')) %>%
  set_names(c('subj_id', 'sex', 'age', 'death')) %>%
  mutate(
    sex = as.factor(c('male', 'female')[sex]),
    age = as.factor(age),
    death = factor(death, levels = 0:4,
                   labels = c('ventilator', 'fastdeath_violent', 'fastdeath_naturalcause', 'intermediatedeath', 'slowdeath'))
  )

# Merge and filter metadata
attphe <- att %>%
  full_join(phe, by = "subj_id") %>%
  unique() %>%
  drop_na()  # Omit missing values
  
# Load and process gene expression data
dat <- fread(file.path(raw_dir, 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct'), select = c('Name', attphe$sample_id)) %>%
  as.data.frame()

# Remove version numbers from gene names
rownames(dat) <- gsub("\\.\\d+$", "", dat$Name)
dat$Name <- NULL

# Convert to matrix and filter based on genes of interest
dat_matrix <- as.matrix(dat)
dat_matrix[!is.finite(dat_matrix)] <- NA # Convert infinite values to NA
dat_matrix <- na.omit(dat_matrix) # Remove rows with NA values

genes_of_interest <- gsub("\\.\\d+$", "", genes_of_interest)

# Filter the data for genes of interest
dat_filtered <- dat_matrix[rownames(dat_matrix) %in% genes_of_interest, ]

# Ensure samples match metadata and filter accordingly
samplesx <- intersect(attphe$sample_id, colnames(dat_filtered))
dat_filtered <- dat_filtered[, samplesx]

# Remove columns with zero variance
dat_filtered <- dat_filtered[, apply(dat_filtered, 2, var) != 0]

# Group samples by tissue based on filtered metadata
samples_by_tissues <- tapply(attphe$sample_id, INDEX = attphe$minor_tissue, FUN = unique)
dat_filtered_tissues <- lapply(samples_by_tissues, function(samps) {
  samps <- intersect(samps, samplesx)
  if (length(samps) > 1) { # Ensure each tissue group has at least 2 samples
    return(dat_filtered[, samps, drop = FALSE])
  } else {
    return(NULL) # Skip tissue groups with insufficient samples
  }
})

# Remove NULL elements (tissues with too few samples)
dat_filtered_tissues <- Filter(Negate(is.null), dat_filtered_tissues)

# Clean Row names
names(dat_filtered_tissues) <- sapply(strsplit(gsub(' ', '', names(dat_filtered_tissues)), '[()]'), function(x) x[[1]])

# Create directory if it does not exist
dir.create(file.path(processed_dir, 'expression', 'readcounts_all'), recursive = TRUE, showWarnings = FALSE)

# Save the filtered data and metadata
sapply(names(dat_filtered_tissues), function(nm) {
  saveRDS(dat_filtered_tissues[[nm]], file.path(processed_dir, 'expression', 'readcounts_all', paste0(nm, '.rds')))
})

# Save the metadata
saveRDS(attphe, file.path(processed_dir, 'attphe_all.rds'))

print("Data preprocessing and saving complete.")


