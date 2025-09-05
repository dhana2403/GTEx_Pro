# Install CRAN packages
cran_packages <- c("data.table", "tidyverse", "dplyr")
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos='https://cloud.r-project.org')
  }
}

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager", repos='https://cloud.r-project.org')
}

bioc_packages <- c("edgeR", "limma", "sva")
BiocManager::install(bioc_packages, ask = FALSE)
