\name{03_sva_batch_correction}
\alias{03_sva_batch_correction}
\title{Batch Effect Correction Using Surrogate Variable Analysis (SVA)}
\description{
  This function performs batch effect correction using Surrogate Variable Analysis (SVA) for RNA-seq data. 
  It removes unwanted batch effects using surrogate variables estimated by SVA, optionally skipping the process for sex-specific tissues.
  The corrected data is saved for each tissue in the `adjusted_sva_all` directory.
}
\usage{
03_sva_batch_correction(processed_dir)
}
\arguments{
  \item{processed_dir}{(character) The path to the directory containing processed data. This directory should have the following subdirectories:
    - `expression/readcounts_tmm_all/`: Contains normalized expression data for tissues in `.rds` format.
    - `attphe_all.rds`: Contains metadata related to sample attributes.}
}
\value{
  This function saves the batch-corrected RNA-seq expression data for each tissue in the directory `expression/adjusted_sva_all/`.
  It also skips tissues with fewer than 20 samples and provides a log of which tissues were processed.
}
\details{
  The function reads in normalized RNA-seq expression data for each tissue from `readcounts_tmm_all/`. It uses SVA to adjust for unwanted batch effects. If there are fewer than 20 samples in a tissue, it is skipped.
  For sex-specific tissues, batch effect correction can be skipped by setting `skip_sva = TRUE`. The adjusted expression data is then saved as `.rds` files.
}
\examples{
# Example usage of the function
processed_dir <- "/path/to/processed"
03_sva_batch_correction(processed_dir)
}
\seealso{
  \code{\link{sva}}, \code{\link{removeBatchEffect}}, \code{\link{num.sv}}
}
\author{
  Dhanalakshmi Jothi
}
\keyword{batch effect correction, RNA-seq, SVA, limma}
