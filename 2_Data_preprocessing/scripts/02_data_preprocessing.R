################################################## DATA PREPROCESSING #################################################

library(data.table)
library(tidyverse)

# Preprocess metadata files - attributes and phenotypes
att <- read_tsv('./data/metadata/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt') %>%
  mutate(SUBJID = sapply(strsplit(SAMPID, '-'), function(x) paste(x[1], x[2], sep = '-'))) %>%
  select(SAMPID, SMTS, SMTSD, SMNABTCH, SMGEBTCH, SUBJID, SMRIN, SMTSISCH) %>%
  set_names(c('sample_id', 'major_tissue', 'minor_tissue', 'batch1', 'batch2', 'subj_id', 'rin', 'ischemic_time'))

phe <- read_tsv('./data/metadata/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt') %>%
  set_names(c('subj_id', 'sex', 'age', 'death')) %>%
  mutate(
    sex = as.factor(c('male', 'female')[sex]),
    age = as.factor(age),
    death = factor(death, levels = 0:4,
                   labels = c('ventilator', 'fastdeath_violent', 'fastdeath_naturalcause', 'intermediatedeath', 'slowdeath'))
  )

attphe <- full_join(att, phe, by = "subj_id")

############################################### ALL ATTRIBUTES AND PHENOTYPES ARE IN ATTPHE ##########################

##Quality control##
##Check for low RIN values 
rin_threshold <- 7.0
     attphe <- attphe %>%
     filter(rin >= rin_threshold) 

# Filter specific tissues of your choice
attphe_filtered = attphe %>%
  filter(minor_tissue %in% c('Brain - Cortex', 'Lung'.......... #####PLEASE INSERT TISSUE OF YOUR CHOICE)) 

############################################ MODIFY COUNT DATA ########################################################
# Load and filter count data
dat <- fread('./data/raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct', select = c('Name', attphe_filtered$sample_id)) %>%
  as.data.frame()

# Strip version numbers from gene names
rownames(dat) <- gsub("\\.\\d+$", "", dat$Name)
dat$Name <- NULL

# Define genes of interest and match them in the dataset
genes_of_interest <- c('ENSG00000198793', .....#PLEASE INSERT GENE ID OF YOUR CHOICE)
genes_of_interest <- gsub("\\.\\d+$", "", genes_of_interest)
dat_filtered <- dat[rownames(dat) %in% genes_of_interest, , drop = FALSE]
samplesx <- intersect(attphe_filtered$sample_id, colnames(dat_filtered))
dat_filtered <- dat_filtered[, samplesx, drop = FALSE]

# Separate data by tissues
dat_tissues <- split(samplesx, attphe_filtered$minor_tissue)
dat_tissues <- lapply(dat_tissues, function(samps) dat_filtered[, samps, drop = FALSE])
names(dat_tissues) <- gsub(' ', '', names(dat_tissues))

# save the processed count and metadata 
dir.create('data/processed/expression/readcounts/', recursive = TRUE, showWarnings = FALSE)
sapply(names(dat_tissues), function(nm) saveRDS(dat_tissues[[nm]], file.path('data/processed/expression/readcounts/', paste0(nm, '.rds'))))
saveRDS(attphe, './data/processed/attphe.rds')

