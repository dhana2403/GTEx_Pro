# GTEx_Pro

Visualizing gene expression in GTEx

## Project Overview

This project contains data acquisition, preprocessing, normalization, batch correction, visualization, and organization code for GTEx dataset. It provides a comprehensive pipeline for analyzing gene expression data from various tissues, facilitating downstream analysis and interpretation.

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

2. Follow the scripts in the src/ folder for detailed steps on:

Data acquisition
Data preprocessing
Normalization
Batch correction
Visualization
Organizing the data for further analysis

The final results and plots will be saved in the results.../ directory.

#### License

This project is licensed under the MIT License. See the LICENSE file for more details.

##### Contributing

Contributions are welcome! 
Please feel free to submit a Pull Request or open an Issue to discuss potential changes or improvements.
Fork the repository.
Open a Pull Request.

###### Acknowledgments

GTEx Project for providing the data.
Partial preprocessing code was adapted from [mdonertas](https://github.com/mdonertas).

