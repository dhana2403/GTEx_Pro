\name{01_data_processing}
\alias{01_data_processing}
\title{Download and Preprocess GTEx RNA-seq Data}
\description{
  This function downloads GTEx RNA-seq metadata and gene expression data,
  processes it, and saves the preprocessed data in the specified output directory.
  The function includes the following steps:
  \itemize{
    \item Creates directories for raw, processed, and metadata data.
    \item Downloads metadata and gene expression data using \code{wget}.
    \item Unzips downloaded files.
    \item Preprocesses metadata by merging different files and filtering missing data.
    \item Filters gene expression data based on a provided list of genes of interest.
    \item Saves the processed gene expression data and metadata to the specified directories.
  }
}
\usage{
01_data_processing(output_dir, genes_of_interest)
}
\arguments{
  \item{output_dir}{(character) The path to the output directory where the data will be saved. This is typically an environment variable (e.g., \code{OUTPUT_DIR}).}
  \item{genes_of_interest}{(character vector) A vector of gene identifiers of interest to filter the RNA-seq expression data.}
}
\value{
  A message indicating that the data has been downloaded, extracted, processed, and saved.
  No value is returned by the function.
}
\examples{
# Example usage of the function
output_dir <- "/path/to/output"
genes_of_interest <- c("ENSG00000198793", "ENSG00000118689")
01_data_processing(output_dir, genes_of_interest)
}
\seealso{
  \code{\link{read_tsv}}, \code{\link{fread}}, \code{\link{saveRDS}}
}
\author{
  Dhanalakshmi Jothi
}
\keyword{data processing}

