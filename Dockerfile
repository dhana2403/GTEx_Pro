# Base image with PyTorch (CUDA or CPU-only)
# Start from an official R image
FROM rocker/r-base:4.3.1

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y unzip r-base && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the Nextflow workflow and R scripts
COPY . .

