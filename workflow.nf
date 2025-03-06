// Process for data acquisition
process data_processing {
    input:
    path output_dir  // The output directory where raw and metadata will be saved
    val genes_of_interest  // The genes of interest passed from Nextflow workflow

    output:
    path "${output_dir}/data/raw/*", emit: raw_data
    path "${output_dir}/data/metadata/*", emit: metadata

    script:
    """
    # Set the OUTPUT_DIR environment variable
    export OUTPUT_DIR=${output_dir}

    # Run the R script to acquire raw data and metadata, passing the genes_of_interest
   Rscript ${output_dir}/01_data_processing/01_data_processing.R ${genes_of_interest.join(' ')}
    """
}

// Process for data normalization
process data_normalization {
    input:
    path processed_dir  // Directory containing processed expression data
    path output_dir

    output:
    path "${processed_dir}/expression/readcounts_tmm_all/*", emit: normalized_data
    path "${processed_dir}/sample_counts.csv", emit: sample_counts

    script:

    """
    # Set OUTPUT_DIR environment variable for the R script
    export OUTPUT_DIR=${output_dir}
    export PROCESSED_DIR=${processed_dir}


    # Run the R script for TMM normalization
    Rscript ${output_dir}/02_data_normalization/02_data_normalization.R
    """
}

process sva_batch_correction {

input:
path processed_dir
path output_dir

output:
path "${processed_dir}/expression/adjusted_sva_all/*", emit: adjusted_data

script:

"""
   # Set OUTPUT_DIR environment variable for the R script
    export PROCESSED_DIR=${processed_dir}

# Run the R script for SVA batch correction
  Rscript ${output_dir}/03_sva_batch_correction/03_sva_batch_correction.R
"""

}


// Workflow definition
workflow {
    // Define the output directory and genes of interest
    output_dir = file("/Users/dhanalakshmijothi/Desktop/GTEx_Pro_final")  // Replace with your actual local path
    processed_dir = file("${output_dir}/data/processed")


    // Define the genes_of_interest (can be passed as a string or list from the command line)
    genes_of_interest = [
        'ENSG00000198793', 'ENSG00000118689', 'ENSG00000096717', 'ENSG00000142082', 
        'ENSG00000133818', 'ENSG00000121691', 'ENSG00000017427', 'ENSG00000140443',  
        'ENSG00000141510', 'ENSG00000077463', 'ENSG00000130203', 'ENSG00000126458', 
        'ENSG00000142168', 'ENSG00000133116'
    ]

    // Run the data acquisition process, passing the output_dir and genes_of_interest
    data_processing(output_dir, genes_of_interest)
    data_normalization(processed_dir, output_dir)
    sva_batch_correction(processed_dir, output_dir)

}