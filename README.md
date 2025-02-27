[![License: All rights reserved]

# GTEx_Pro

Visualizing gene expression in GTEx

## Project Overview

This project contains data acquisition, preprocessing, normalization, batch correction, visualization, and organization code for GTEx dataset. It provides a comprehensive workflow for analyzing gene expression data from various tissues, facilitating downstream analysis and interpretation. This workflow implements a robust and accurate preprocessing pipeline for analyzing gene expression data of GTEx effectively correcting for batch effects while retaining biological variations associated with sex (view the 3D interactive PCA plot below) across 44 GTEx tissues. The code is currently being automated using Nextflow.

## Table of Contents

- [Project Overview](#project-overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [License](#license)
- [Contributing](#contributing)
- [Acknowledgments](#acknowledgments)

## Requirements

To run this project, you will need:

- R version >= 4.0
- R packages:
  - `dplyr`
  - `ggplot2`
  - `edgeR`
  - `limma`
  - `sva`
  - `officer`
  - `magrittr`
  - Other dependencies as necessary

## Installation

Install required R packages (if not already installed):

### Usage

1.  ![Flowchart](GTEx_Pro_pipeline.drawio.pdf)

2. Follow the order of the scripts detailed steps on:

[Data acquisition](https://github.com/dhana2403/GTEx_Pro/blob/main/1_Data_acquisition/scripts/01_data_acquisition.R)
[Data preprocessing](https://github.com/dhana2403/GTEx_Pro/blob/main/2_Data_preprocessing/scripts/02_data_preprocessing.R)
[Normalization](https://github.com/dhana2403/GTEx_Pro/blob/main/3_Data_normalization/scripts/03_data_normalization_TMM.R)
[Batch correction](https://github.com/dhana2403/GTEx_Pro/blob/main/4_Batch%20correction/scripts/04_01_SVA_batch_correction.R)Visualization
Organizing the data for further analysis

The final results and plots will be saved in the results.../ directory.

## PCA 

[View the 3D PCA plot](https://dhana2403.github.io/3D_plots/pca_plot_all_tissues_sva_3D.html) 
[View 3D UMAP plot](https://dhana2403.github.io/3D_plots/umap_plot_all_tissues_sva_3D.html)
[View 3D TNSE plot](https://dhana2403.github.io/3D_plots/tsne_plot_all_tissues_sva_3D.html)

#### License

This project is licensed under the **All Rights Reserved** license for the time being.
If you'd like to request permission for use or modification, please contact [dhana2403](https://github.com/dhana2403)

##### Contributing

If you'd like to offer feedback on the workflow/code, please contact [dhana2403](https://github.com/dhana2403)

###### Acknowledgments

GTEx Project for providing the data.
Partial preprocessing code was adapted from [mdonertas](https://github.com/mdonertas).

