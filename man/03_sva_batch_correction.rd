\name{03_sva_batch_correction}
\alias{03_sva_batch_correction}
\title{Batch Effect Correction Using Surrogate Variable Analysis (SVA)}
\description{
  This script performs batch effect correction on normalized RNA-seq data using Surrogate Variable Analysis (SVA). It identifies and removes unwanted sources of variation such as batch effects, while skipping SVA for sex-specific tissues. Adjusted expression data is saved per tissue in the `adjusted_sva_all` directory.
}
\usage{
03_sva_batch_correction(processed_dir)
}
\arguments{
  \item{processed_dir}{(character) The path to the directory containing processed data. This directory must contain:
    \itemize{
      \item \code{expression/readcounts_tmm_all/}: Directory with TMM-normalized expression matrices (.rds files).
      \item \code{attphe_all.rds}: Metadata file with sample information including \code{sex} and \code{batch1}.
    }
  }
}
\value{
  For each tissue, a batch-adjusted expression matrix is saved in the subdirectory \code{expression/adjusted_sva_all/}. Tissues with fewer than 20 samples are skipped. SVA is skipped for predefined sex-specific tissues, where only batch correction is applied.
}
\details{
  For each tissue, the function loads the normalized read counts and corresponding metadata. It then estimates the number of surrogate variables using the \code{be} (beating the eigenvalue) method and applies SVA to estimate hidden confounders. These are then used with \code{removeBatchEffect} from the \pkg{limma} package to adjust the expression data. 
  
  Sex-specific tissues (e.g., \code{Testis}, \code{Ovary}, \code{Prostate}, etc.) are excluded from SVA adjustment and instead corrected only for known batch effects.
}
\examples{
# Example usage
Sys.setenv(PROCESSED_DIR = "/path/to/processed")
03_sva_batch_correction(Sys.getenv("PROCESSED_DIR"))
}
\seealso{
  \code{\link[sva]{sva}}, \code{\link[sva]{num.sv}}, \code{\link[limma]{removeBatchEffect}}
}
\author{
  Dhanalakshmi Jothi
}
\keyword{RNA-seq, batch effect correction, surrogate variable analysis, SVA, limma, expression data}
