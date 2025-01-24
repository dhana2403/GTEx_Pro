#######################GTEx data retrieval#####################


#Load library for acquiring metadata and gene expression data
import os      #creates directories
import requests  #downloads files
import subprocess #unzips command

#create directories for saving the data
os.makedirs("data/raw", exist_ok=True)
os.makedirs("data/processed", exist_ok=True)
os.makedirs("data/metadata", exist_ok=True)

#assign variables to the urls to be downloaded
urls = [
    "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx",
    "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx",
    "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
    "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
]

#Function for getting the data from the website

for url in urls:

 response = requests.get(url)

 filename = os.path.basename(url)

 if "metadata" in url:
     save_dir = "data/metadata"
 else:
     save_dir = "data/raw"

 file_path = os.path.join(save_dir, filename)

 with open(file_path, 'wb') as f:
     f.write(response.content)
     

#####unzip files 

if filename.endswith(".gz"):

 unzip_command = f"gunzip {file_path}"

 subprocess.run(unzip_command, shell=True, check=True)

 print(f"Unzipped file: {filename}")



