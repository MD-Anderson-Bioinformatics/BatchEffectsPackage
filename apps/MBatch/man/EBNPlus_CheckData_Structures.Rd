\name{EBNPlus_CheckData_Structures}
\alias{EBNPlus_CheckData_Structures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{EBNPlus_CheckData_Structures}
\description{
Take two matrices and check that they will work as arguments to MBatch EBNPlus functions.
}
\usage{
EBNPlus_CheckData_Structures(theDataMatrix1, theDataMatrix2, theDataReplicates1 = NULL,
theDataReplicates2 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theDataMatrix1}{A matrix for data set 1 containing numeric values with columns being
  sample ids and rows being gene ids.}
  \item{theDataMatrix2}{A matrix for data set 2 containing numeric values with columns being
  sample ids and rows being gene ids.}
  \item{theDataReplicates1}{A vector of "replicates" in data set 1 used for corrections.
  Defaults to NULL. See "details" below.}
  \item{theDataReplicates2}{A vector of "replicates" in data set 2 used for corrections.
  Defaults to NULL. See "details" below.}
}
\details{
Performs a stop if conditions are not met.
This function checks the following:
Both matrix arguments pass is.matrix test.
Both matrices have column names.
Both matrices have row names and they interect at least once.
Both matrices have column names that intersect at least once or the replicate value vectors
are the same size and exist in the column names.
All data in the matrices is numeric.

Please note, support for performing corrections using replicate lists is not yet supported.
}
\value{
No meaningful value. Performs a stop if conditions are not met.
}
\references{}
\author{}
\note{}
\seealso{}
\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
