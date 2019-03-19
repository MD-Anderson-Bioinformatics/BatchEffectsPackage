\name{RBN_Pseudoreplicates}
\alias{RBN_Pseudoreplicates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{RBN_Pseudoreplicates}
\description{
Do RBN correction using replicates with different identifiers (pseudoreplicates).
}
\usage{
RBN_Pseudoreplicates(theInvariantMatrix, theVariantMatrix,
theInvariantReplicates, theVariantReplicates, theInvariantGroupId = "",
theVariantGroupId = "", theMatchedReplicatesFlag = TRUE,
theCombineOnlyFlag = FALSE, thePath = NULL, theWriteToFile = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theInvariantMatrix}{Matrix with sample names in colnames and features (like genes) in rownames. This matrix is invariant.}
  \item{theVariantMatrix}{Matrix with sample names in colnames and features (like genes) in rownames. This matrix is variant.}
  \item{theInvariantReplicates}{Vector of feature ids indicating replicates for variant data.}
  \item{theVariantReplicates}{Vector of feature ids indicating replicates for invariant data.}
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
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
