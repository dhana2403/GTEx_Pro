
################################Data acquisition#################################################

#create directories 
system('mkdir -p data/raw')
system('mkdir -p data/processed')
system('mkdir -p data/metadata')
system('mkdir -p scripts_gtex')

#download the metadata and gene expression files from https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression and move the metadata to data/metadata/ and gene expression data to data/raw

