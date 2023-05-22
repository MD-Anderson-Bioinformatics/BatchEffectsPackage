\name{writeAsGenericMatrix}
\alias{writeAsGenericMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{writeAsGenericMatrix}
\description{
This function takes a matrix and writes it out as a cross-platform tab-delimited matrix file with column and row labels.
The matrix data should be doubles.
If row labels are included, there is an initial tab on the column label row.
}
\usage{
writeAsGenericMatrix(theFile, theMatrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theFile}{Full path to the output file.}
  \item{theMatrix}{The matrix to write out--requires column names and double (numeric) data.}
}
\value{
Returns TRUE if the save succeeded.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
#writeAsGenericMatrix <- function(theFile, theMatrix)
#{
#	write.table(theMatrix, file=theFile, quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)
#	return(TRUE)
#}
}
