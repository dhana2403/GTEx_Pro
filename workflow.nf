#!/usr/bin/env nextflow

// Process for Data Acquisition
process data_acquisition {
    input:
    path output_dir  // The output directory where raw and metadata will be saved

    output:
    path 'data/raw/*', emit: raw_data
    path 'data/metadata/*', emit: metadata

    script:
    """
    # Ensure that the data/raw and data/metadata directories exist
    mkdir -p ${output_dir}/data/raw
    mkdir -p ${output_dir}/data/metadata

    # Code to acquire your raw data (adjust as needed)
    python ${output_dir}/1_Data_acquisition/scripts/01_data_acquisition.R --output ${output_dir}/data/raw

    # Code to acquire your metadata (adjust as needed)
    python ${output_dir}/1_Data_acquisition/scripts/01_data_acquisition.R --output ${output_dir}/data/metadata
    """
}

// Workflow definition
workflow {
    output_dir = file("/Users/dhanalakshmijothi/GTEx_Pro")  // Path to your local repository or directory

    // Run the data acquisition process
    data_acquisition(output_dir)
}
