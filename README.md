[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# GTEx_Pro

Visualizing gene expression in GTEx

## Project Overview

This project contains data acquisition, preprocessing, normalization, batch correction, visualization, and organization code for GTEx dataset. It provides a comprehensive pipeline for analyzing gene expression data from various tissues, facilitating downstream analysis and interpretation. This code implements a robust pipeline for analyzing gene expression data of GTEx effectively correcting for batch effects while retaining biological variations associated with sex. The code isn't fully automated and requires manual inputs at specific points. Future updates will be focused on automating the code with less manual inputs

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
Normalization
Batch correction (sva-recommended)
Visualization
Organizing the data for further analysis

The final results and plots will be saved in the results.../ directory.

## PCA Analysis

[View the 3D PCA plot](https://dhana2403.github.io/GTEx_Pro/pca_plot_all_tissues_sva_3D.html)

#### License

This project is licensed under the MIT License. See the LICENSE file for more details.

##### Contributing

Contributions are welcome! 
Please feel free to submit a Pull Request or open an Issue to discuss potential changes or improvements.
Fork the repository.
Open a Pull Request.
Offer your feedback

###### Acknowledgments

GTEx Project for providing the data.
Partial preprocessing code was adapted from [mdonertas](https://github.com/mdonertas).

