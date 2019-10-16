\name{RBN_Replicates}
\alias{RBN_Replicates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{RBN_Replicates}
\description{
Do RBN correction using exact replicates.
}
\usage{
RBN_Replicates(theInvariantMatrix, theVariantMatrix, theInvariantGroupId = "",
theVariantGroupId = "", theMatchedReplicatesFlag = TRUE,
theCombineOnlyFlag = FALSE, thePath = NULL, theWriteToFile = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theInvariantMatrix}{Matrix with sample names in colnames and features (like genes) in rownames. This matrix is invariant.}
  \item{theVariantMatrix}{Matrix with sample names in colnames and features (like genes) in rownames. This matrix is variant.}
  \item{theInvariantGroupId}{Group name used for labelling invariant features when combining matrixes.}
  \item{theVariantGroupId}{Group name used for labelling variant features when combining matrixes.}
  \item{theMatchedReplicatesFlag}{If TRUE, indicates that NAs should be added for missing replicates. Defaults to FALSE.}
  \item{theCombineOnlyFlag}{If TRUE, only combined the matrixes, do not correct. Defaults to FALSE.}
  \item{thePath}{Location for output.}
  \item{theWriteToFile}{TRUE means write corrected data to thePath.}
}
\details{Uses RBN (Replicates Based Normalization) to correct one data group to match another.
}
\value{Returned the corrected matrix unless thePath is given and theWriteToFile is TRUE, in which case the full path and filename to the corrected data file is returned.}
