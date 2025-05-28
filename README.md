<table>
  <tr>
    <td><h1>GTEx Pro Pipeline </h1></td>
    <td><img src="https://raw.githubusercontent.com/dhana2403/GTEx_sample/main/2.png" alt="GTEx Pro Logo" width="200" /></td>
  </tr>
</table>

### GTEx_pro pipeline - Academic Version
This repository contains the academic version of the GTEx_pro pipeline for RNA-seq data analysis. It provides a comprehensive workflow for analyzing gene expression data from various tissues, facilitating downstream analysis and interpretation. This workflow implements a robust and accurate preprocessing pipeline for analyzing gene expression data of GTEx effectively correcting for batch effects while using sex as a covariate

**Features:**
- Integrates basic normalization using TMM followed by CPM scaling and batch correction using SVA.
- Provides reproducible preprocessing pipeline for tissue-specific gene expression analysis.


**Key Terms:**
- The pipeline is provided free of charge for academic use. If used for academic purposes, please cite https://www.biorxiv.org/content/10.1101/2025.04.26.650748v1 
- Redistribution or modification for commercial purposes is strictly prohibited.

**Prerequisites**

1. Install Nextflow:
   For Linux or macOS:
       curl -s https://get.nextflow.io | bash

2. Install R:
   Download R from https://www.r-project.org/ and install necessary packages (data.table, tidyverse, dplyr, edgeR, limma, sva)
   
3. Clone this repository:
   git clone <repository_url>

4. Configure the workflow:
   cd <repository_path>
   In the workflow.nf file, you'll need to specify the output_dir and the genes of interest:

6. Customize your genes_of_interest in workflow.nf:
   Replace the example gene IDs in the genes_of_interest variable with your own gene IDs.
   Make sure the gene IDs are separated by commas and enclosed in single quotes, e.g., 'ENSG00000000001', 'ENSG00000000002'.

7. Run the workflow:
   nextflow run workflow.nf
   
8. Results in the output directory:
   The results will be saved in the directories defined by output_dir and processed_dir within the workflow.nf file.
   Specifically, you will find:
         Raw data and metadata in {output_dir}/data/raw/ and ${output_dir}/data/metadata/.
         Normalized expression data in {processed_dir}/expression/readcounts_tmm_all/.
         Adjusted (SVA-corrected) expression data in {processed_dir}/expression/adjusted_sva_all/.
   
**Results**

   <div style="display: flex; justify-content: space-between; gap: 30px;">
   <img src = "https://github.com/dhana2403/GTEx_sample/blob/main/images/TMM%2BCPM.png?raw=true" width="20%">
   <img src = "https://github.com/dhana2403/GTEx_sample/blob/main/images/TMM%2BCPM%2BSVA.png?raw=true" width="40%">
   </div>
     
**Citation**

Jothi, D., GTEx_Pro: A robust and accurate preprocessing pipeline for GTEx data: TMM+CPM normalization and SVA batch correction for enhanced multi-tissue analysis. bioRxiv, 2025: p. 2025.04.26.650748.





