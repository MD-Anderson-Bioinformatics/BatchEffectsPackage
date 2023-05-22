\name{writeAsGenericDataframe}
\alias{writeAsGenericDataframe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{writeAsGenericDataframe}
\description{
This function takes a data.frame and writes it out as a cross-platform
tab-delimited matrix file with column and row labels.
The data.frame data is handled as strings.
If row labels are included, there is an initial tab on the column label row.
}
\usage{
writeAsGenericDataframe(theFile, theDataframe)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theFile}{Full path to the output file.}
  \item{theDataframe}{The data.frame to write out--requires column names.}
}
\value{
Returns TRUE if the save succeeded.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
#writeAsGenericDataframe <- function(theFile, theDataframe)
#{
#	write.table(theDataframe, file=theFile, quote=FALSE,
#         sep="\t", col.names=TRUE, row.names=FALSE)
#	return(TRUE)
#}
}
