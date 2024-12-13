
################################Data acquisition#################################################

ystem('mkdir -p data/raw')
system('mkdir -p data/processed')
system('mkdir -p data/metadata')

library(tidyverse)

system('wget -nv https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx -P data/metadata')
system('wget -nv https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx -P data/metadata')
system('wget -nv https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt -P data/metadata')
system('wget -nv https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt -P data/metadata')
system('wget -nv https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz -P data/raw')
system('cd data/raw; gunzip *')


