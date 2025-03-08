# Load required libraries
library(data.table)
library(tidyverse)
library(dplyr)
library(argparse)  # Argument parsing library

# Set up argument parser
parser <- ArgumentParser(description = "Preprocess Gene Expression Data")
parser$add_argument("--output", type = "character", help = "Directory to store output")
parser$add_argument("--gene", type = "character", help = "Comma-separated list of gene IDs")

# Parse the command line arguments
args <- parser$parse_args()

# Ensure output directory exists
dir.create(args$output, recursive = TRUE, showWarnings = FALSE)

# Parse gene IDs from the string passed in command line
genes_of_interest <- strsplit(args$gene, ",")[[1]]

################################################## Load and process metadata files #################################################

att <- read_tsv('./data/metadata/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt') %>%
  mutate(SUBJID = sapply(strsplit(SAMPID, '-'), function(x) paste(x[1], x[2], sep = '-'))) %>%
  select(SAMPID, SMTS, SMTSD, SMNABTCH, SMGEBTCH, SUBJID, SMRIN, SMTSISCH) %>%
  set_names(c('sample_id', 'major_tissue', 'minor_tissue', 'batch1', 'batch2', 'subj_id', 'rin', 'ischemic_time')) %>%
  unique()

phe <- read_tsv('./data/metadata/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt') %>%
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
  unique()

################################################## All attributes and phenotypes are in attphe #################################################

attphe_filtered <- attphe %>%
  # filter(rin >= 7.0) %>%  # include it as a quality control if there are many samples
  filter(minor_tissue %in% c('Brain - Cortex')) %>%  # Example of filtering based on a specific tissue
  filter(death %in% c('intermediatedeath', 'slowdeath'))  # Example of filtering based on experimental conditions

################################################## Modify count data #################################################

# Load and process gene expression data
dat <- fread('./data/raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct', select = c('Name', attphe_filtered$sample_id)) %>%
  as.data.frame()

# Remove version numbers from gene names
rownames(dat) <- gsub("\\.\\d+$", "", dat$Name)
dat$Name <- NULL

# Convert to matrix and filter based on genes of interest (using the input gene list)
dat_matrix <- as.matrix(dat)
genes_of_interest <- gsub("\\.\\d+$", "", genes_of_interest)  # Clean up gene IDs if necessary
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
  if (length(samps) > 1) {  # Ensure each tissue group has at least 2 samples
    return(dat_filtered[, samps, drop = FALSE])
  } else {
    return(NULL)  # Skip tissue groups with insufficient samples
  }
})

# Remove NULL elements (tissues with too few samples)
dat_filtered_tissues <- Filter(Negate(is.null), dat_filtered_tissues)

# Clean column names
names(dat_filtered_tissues) <- sapply(strsplit(gsub(' ', '', names(dat_filtered_tissues)), '[()]'), function(x) x[[1]])

# Save the filtered data and metadata
dir.create(file.path(args$output, 'data/processed/expression/readcounts/'), recursive = TRUE, showWarnings = FALSE)
sapply(names(dat_filtered_tissues), function(nm) {
  saveRDS(dat_filtered_tissues[[nm]], file.path(args$output, 'data/processed/expression/readcounts/', paste0(nm, '.rds')))
})
saveRDS(attphe_filtered, file.path(args$output, 'data/processed/attphe.rds'))



