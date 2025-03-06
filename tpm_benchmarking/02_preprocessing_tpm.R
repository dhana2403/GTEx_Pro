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

################################################## Modify count data #################################################  

# Define genes of interest  
genes_of_interest <- c('ENSG00000198793', 'ENSG00000118689', 'ENSG00000096717', 'ENSG00000142082',  
                       'ENSG00000133818', 'ENSG00000121691', 'ENSG00000017427', 'ENSG00000140443',  
                       'ENSG00000141510', 'ENSG00000077463', 'ENSG00000130203', 'ENSG00000126458',  
                       'ENSG00000142168', 'ENSG00000133116')  

dat_filtered <- fread('./data/raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct',  
                      select = c('Name', attphe$sample_id)) %>%  
  mutate(Name = gsub("\\.\\d+$", "", Name)) %>%
  filter(Name %in% genes_of_interest) %>%  
  as.data.frame()  

rownames(dat_filtered) <- dat_filtered$Name  
dat_filtered$Name <- NULL  

# Ensure samples match metadata and filter accordingly  
samplesx <- intersect(attphe$sample_id, colnames(dat_filtered))  
dat_filtered <- dat_filtered[, samplesx, drop = FALSE]  


# Group samples by tissue based on filtered metadata  
samples_by_tissues <- split(attphe$sample_id, attphe$minor_tissue)  
dat_filtered_tissues <- lapply(samples_by_tissues, function(samps) {  
  samps <- intersect(samps, samplesx)  
  if (length(samps) > 1) {  
    return(dat_filtered[, samps, drop = FALSE])  
  } else {  
    return(NULL)  
  }  
})  

# Remove NULL elements (tissues with too few samples)  
dat_filtered_tissues <- Filter(Negate(is.null), dat_filtered_tissues)  

# Clean column names  
names(dat_filtered_tissues) <- sapply(strsplit(gsub(' ', '', names(dat_filtered_tissues)), '[()]'), function(x) x[[1]])  

# Create directory if it does not exist  
dir.create('data/processed/expression/tpm_all/', recursive = TRUE, showWarnings = FALSE)  

# Save the filtered data and metadata  
sapply(names(dat_filtered_tissues), function(nm) {  
  saveRDS(dat_filtered_tissues[[nm]], file.path('data/processed/expression/tpm_all/', paste0(nm, '.rds')))  
})  

saveRDS(attphe, './data/processed/attphe_tpm_all.rds')  
