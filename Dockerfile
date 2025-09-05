# Use stable Rocker R image based on Debian stable
FROM rocker/r-ver:4.3.1

# Set working directory inside container
WORKDIR /project

# Install essential system dependencies for R packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        wget \
        gzip \
        unzip \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy your R requirements file into the container
COPY r_requirements.R .

# Install all R packages from your requirement file
RUN Rscript r_requirements.R

# Copy your Nextflow workflow and R scripts
COPY . .

# Default command (can be overridden when running container)
CMD ["nextflow", "run", "workflow.nf"]
