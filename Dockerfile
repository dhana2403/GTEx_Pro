# Use stable R image with tidyverse pre-installed
FROM rocker/tidyverse:4.3.1

# Set working directory
WORKDIR /project

# Install minimal system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends unzip wget gzip && \
    rm -rf /var/lib/apt/lists/*

# Install Bioconductor packages
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install(c('edgeR','limma','sva'))"

# Copy Nextflow workflow and R scripts
COPY . .

# Default command (Nextflow will handle execution)
CMD ["R"]
