#!/usr/bin/env nextflow

params.tissue = 'Brain - Cortex', 'Lung', 'Liver' // Tissue of interest
params.gene_id = 'ENSG00000158793' // Gene of interest

// Process for Data Acquisition
process data_acquisition {
    input:
    path output_dir  // The output directory where raw and metadata will be saved

    output:
    path 'data/raw/*', emit: raw_data
    path 'data/metadata/*', emit: metadata

    script:
    """
    # Code to acquire your raw and metadata
    """
}

// Workflow definition
workflow {
    output_dir = file("/Users/dhanalakshmijothi/GTEx_Pro")  // Replace with your actual local path 

    // Run the data acquisition process
    data_acquisition(output_dir)
}
