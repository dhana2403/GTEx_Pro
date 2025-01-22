

// Process for Data Acquisition
process data_acquisition {
    input:
    path output_dir  // The output directory where raw and metadata will be saved

    output:
    path 'data/raw/*', emit: raw_data
    path 'data/metadata/*', emit: metadata

    script:
    """
    # Ensure that the data/raw and data/metadata directories exist within the output directory
    mkdir -p ${output_dir}/data/raw
    mkdir -p ${output_dir}/data/metadata

    # Check if directories are created
    echo "Checking raw directory: ${output_dir}/data/raw"
    ls ${output_dir}/data/raw

    # Run the R script to acquire raw data
    Rscript ${output_dir}/1_Data_acquisition/scripts/01_data_acquisition.R --output ${output_dir}/data/raw

    # Check if raw data was downloaded
    ls ${output_dir}/data/raw

    # Run the R script to acquire metadata
    Rscript ${output_dir}/1_Data_acquisition/scripts/01_data_acquisition.R --output ${output_dir}/data/metadata

    # Check if metadata was downloaded
    ls ${output_dir}/data/metadata
    """
}

// Workflow definition
workflow {
<<<<<<< HEAD
    output_dir = file("/Users/dhanalakshmijothi/GTEx_Pro")  // Replace with your actual local path 
=======
    output_dir = file("/Users/dhanalakshmijothi/GTEx_Pro")  // Path to your local repository or directory
>>>>>>> origin/main

    // Run the data acquisition process
    data_acquisition(output_dir)
}
