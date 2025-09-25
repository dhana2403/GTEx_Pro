<table width="100%">
  <tr>
    <td style="vertical-align: middle;">
      <h1>GTEx Pro Pipeline</h1>
    </td>
    <td style="text-align: right; padding-right: 60px; vertical-align: middle;">
      <img src="https://raw.githubusercontent.com/dhana2403/GTEx_sample/main/2.png" alt="GTEx Pro Logo" width="200" />
    </td>
  </tr>
</table>

### Overview
This repository contains the academic version of the GTEx_pro pipeline for GTEx data analysis. It provides a comprehensive workflow for analyzing gene expression data from GTEx tissues, and improves the precision and reliability of the downstream gene expression analysis by stabilizing tissue-specific gene signals through the integration of robust normalization and batch correction methods

<p align="center">
  <img src="https://github.com/dhana2403/GTEx_sample/blob/main/images/GTEx_Pro_pipeline.png?raw=true" width="100%">
</p>

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

7. Build the Docker image:
   docker build -t myproject:latest .

8. Run the workflow:
   nextflow run workflow.nf -with-docker myproject:latest
   
8. Results in the output directory:
   The results will be saved in the directories defined by output_dir and processed_dir within the workflow.nf file.
   Specifically, you will find:
         Raw data and metadata in {output_dir}/data/raw/ and ${output_dir}/data/metadata/.
         Normalized expression data in {processed_dir}/expression/readcounts_tmm_all/.
         Adjusted (SVA-corrected) expression data in {processed_dir}/expression/adjusted_sva_all/.
   
**Results**

   <div style="display: flex; justify-content: space-between; gap: 90px;">
   <img src = "https://github.com/dhana2403/GTEx_sample/blob/main/images/TMM%2BCPM.png?raw=true" width="30%">
   <img src = "https://github.com/dhana2403/GTEx_sample/blob/main/images/TMM%2BCPM%2BSVA.png?raw=true" width="38%">
   </div>
     
**Citation**

Jothi, D. GTEx pro enables accurate multi-tissue gene expression analysis using robust normalization and batch correction. Sci Rep 15, 32684 (2025). https://doi.org/10.1038/s41598-025-20697-0.

**News**

Interactive web interface is under development for non-programmers. Acess it via https://gtexprov1.streamlit.app/
