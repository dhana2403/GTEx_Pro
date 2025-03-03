process data_acquisition {
    input:
    path output_dir  // The output directory where raw and metadata will be saved

    output:
    path "${output_dir}/data/raw/*", emit: raw_data
    path "${output_dir}/data/metadata/*", emit: metadata

    script:
    """
    # Set the OUTPUT_DIR environment variable
    export OUTPUT_DIR=${output_dir}

    # Run the R script to acquire raw data and metadata
    Rscript ${output_dir}/01_data_acquisition/01_data_acquisition.R
    """
}

// Workflow definition
workflow {
    output_dir = file("/Users/dhanalakshmijothi/Desktop/GTEx_Pro_final")  // Replace with your actual local path

    // Run the data acquisition process
    data_acquisition(output_dir)
}
