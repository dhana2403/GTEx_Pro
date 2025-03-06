################################################## DATA PREPROCESSING #################################################

# Load required libraries
library(data.table)
library(tidyverse)
library(dplyr)

# Load and process metadata files
att <- read_tsv('./data/metadata/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt') %>%
  mutate(SUBJID = sapply(strsplit(SAMPID, '-'), function(x) paste(x[1], x[2], sep = '-'))) %>%
  select(SAMPID, SMTS, SMTSD, SMNABTCH, SMGEBTCH, SUBJID, SMRIN, SMTSISCH) %>%
  set_names(c('sample_id', 'major_tissue', 'minor_tissue', 'batch1', 'batch2', 'subj_id', 'rin', 'ischemic_time')) %>%
  unique()

print(att$minor_tissue)


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
  unique() %>%
  drop_na() # Omit missing values

################################################## All attributes and phenotypes are in attphe #################################################

attphe_filtered <- attphe %>%
  filter(!minor_tissue %in% c("Cervix - Ectocervix", "Cervix - Endocervix", "Fallopian Tube",
                              "Testis", "Uterus", "Vagina", "Ovary", "Prostate","Breast - Mammary Tissue"))
################################################## Modify count data #################################################

# Load and process gene expression data
dat <- fread('./data/raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct', select = c('Name', attphe_filtered$sample_id)) %>%
  as.data.frame()

# Remove version numbers from gene names
rownames(dat) <- gsub("\\.\\d+$", "", dat$Name)
dat$Name <- NULL

# Convert to matrix and filter based on genes of interest
dat_matrix <- as.matrix(dat)
dat_matrix[!is.finite(dat_matrix)] <- NA # Convert infinite values to NA
dat_matrix <- na.omit(dat_matrix) # Remove rows with NA values

genes_of_interest <- c('ENSG00000142208', 'ENSG00000082701', 'ENSG00000135414', 'ENSG00000150907', 'ENSG00000130766', 'ENSG00000177169', 'ENSG00000096088', 'ENSG00000158828', 'ENSG00000140992', 'ENSG00000171791', 'ENSG00000100292', 'ENSG00000214253', 'ENSG00000204490','ENSG00000143799')

genes_of_interest <- gsub("\\.\\d+$", "", genes_of_interest)

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

# Clean column names
names(dat_filtered_tissues) <- sapply(strsplit(gsub(' ', '', names(dat_filtered_tissues)), '[()]'), function(x) x[[1]])

# Create directory if it does not exist
dir.create('data/processed/expression/readcounts_all_cross/', recursive = TRUE, showWarnings = FALSE)

# Save the filtered data and metadata
sapply(names(dat_filtered_tissues), function(nm) {
  saveRDS(dat_filtered_tissues[[nm]], file.path('data/processed/expression/readcounts_all_cross/', paste0(nm, '.rds')))
})
saveRDS(attphe_filtered, './data/processed/attphe_all_cross.rds')