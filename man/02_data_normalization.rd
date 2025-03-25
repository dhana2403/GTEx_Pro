\name{02_data_normalization}
\alias{02_data_normalization}
\title{Normalize Gene Expression Data Using TMM from edgeR}
\description{
  This function performs normalization of gene expression data using the TMM (Trimmed Mean of M-values) method from the edgeR package.
  The function processes tissue-specific read count data, applies TMM normalization, and saves the normalized expression data.
  It also generates a CSV file summarizing the number of valid samples processed for each tissue.
}
\usage{
02_data_normalization(processed_dir)
}
\arguments{
  \item{processed_dir}{(character) The path to the directory containing processed data. This is typically an environment variable (e.g., \code{PROCESSED_DIR}). The directory should contain `expression/readcounts_all/` with the raw read count data and `attphe_all.rds` for metadata.}
}
\value{
  A CSV file containing the sample counts for each tissue after normalization is saved at the specified location in \code{processed_dir}. Additionally, the normalized data for each tissue is saved in the directory \code{expression/readcounts_tmm_all/}.
}
\examples{
# Example usage of the function
processed_dir <- "/path/to/processed"
02_data_normalization(processed_dir)
}
\seealso{
  \code{\link{DGEList}}, \code{\link{calcNormFactors}}, \code{\link{cpm}}
}
\author{
  Dhanalakshmi Jothi
}
\keyword{data normalization, RNA-seq, TMM, edgeR}
