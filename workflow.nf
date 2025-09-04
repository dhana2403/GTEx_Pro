// -------------------------------
// Channels
// -------------------------------

// Define static input channels
output_dir_ch = Channel.value(file("/Users/dhanalakshmijothi/Desktop/GTEx_sample"))

// Genes of interest as a value channel
genes_of_interest_ch = Channel.value([
    'ENSG00000198793', 'ENSG00000118689', 'ENSG00000096717', 'ENSG00000142082',
    'ENSG00000133818', 'ENSG00000121691', 'ENSG00000017427', 'ENSG00000140443',
    'ENSG00000141510', 'ENSG00000077463', 'ENSG00000130203', 'ENSG00000126458',
    'ENSG00000142168', 'ENSG00000133116'
])

// -------------------------------
// Processes
// -------------------------------

process data_processing {
    input:
    path output_dir
    val genes_of_interest

    output:
    path "${output_dir}/data/raw/*"      , emit: raw_data
    path "${output_dir}/data/metadata/*" , emit: metadata
    path "${output_dir}/data/processed"  , emit: processed_dir

    script:
    """
    export OUTPUT_DIR=${output_dir}
    Rscript ${output_dir}/01_data_processing/01_data_processing.R ${genes_of_interest.join(' ')}
    """
}

process data_normalization {
    input:
    path processed_dir
    path output_dir

    output:
    path "${processed_dir}/expression/readcounts_tmm_all/*" , emit: normalized_data
    path "${processed_dir}/sample_counts.csv"               , emit: sample_counts

    script:
    """
    export OUTPUT_DIR=${output_dir}
    export PROCESSED_DIR=${processed_dir}

    Rscript ${output_dir}/02_data_normalization/02_data_normalization.R
    """
}

process sva_batch_correction {
    input:
    path processed_dir
    path output_dir

    output:
    path "${processed_dir}/expression/adjusted_sva_all/*" , emit: adjusted_data

    script:
    """
    export PROCESSED_DIR=${processed_dir}
    Rscript ${output_dir}/03_sva_batch_correction/03_sva_batch_correction.R
    """
}

// -------------------------------
// Workflow
// -------------------------------
workflow {
    // Run data acquisition
    data_processing(output_dir_ch, genes_of_interest_ch)

    // Feed processed_dir into normalization
    data_normalization(data_processing.out.processed_dir, output_dir_ch)

    // Feed processed_dir into batch correction
    sva_batch_correction(data_processing.out.processed_dir, output_dir_ch)
}
