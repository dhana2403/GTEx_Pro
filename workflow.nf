

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

process data_preprocessing {
    input:
    path output_dir
    val geneid_of_choice

    output:
    path 'data/processed/expression/readcounts/*', emit: processed_data
    path 'data/processed/attphe.rds', emit: attphe_data

    script:
    """
    # Run the R script for each gene in the list
    for gene in ${geneid_of_choice} ; do
        Rscript ${output_dir}/2_Data_preprocessing/scripts/02_data_preprocessing.R \
            --output ${output_dir}/data/processed \
            --gene "${gene}"
    done
    """
}

workflow {
    output_dir = file("/Users/dhanalakshmijothi/GTEx_Pro")  // Replace with your actual local path 
    geneid_of_choice = ['ENSG00000198793', 'ENSG00000118689', 'ENSG00000096717', 'ENSG00000142082', 'ENSG00000133818', 
                        'ENSG00000121691', 'ENSG00000017427', 'ENSG00000140443', 'ENSG00000141510', 'ENSG00000077463', 
                        'ENSG00000130203', 'ENSG00000126458', 'ENSG00000142168', 'ENSG00000133116']
    tissue_of_choice = 'Brain'  // Example of a defined tissue of choice

    // Run the data acquisition process
    data_acquisition(output_dir)

    // Run the data preprocessing process after acquisition
    data_preprocessing(output_dir, geneid_of_choice)
}

