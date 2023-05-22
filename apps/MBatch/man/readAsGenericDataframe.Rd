\name{readAsGenericDataframe}
\alias{readAsGenericDataframe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{readAsGenericDataframe}
\description{
This function takes a cross-platform tab-delimited matrix file with column
and optional row labels and reads it in as a data.frame.
The data.frame data is handled as strings.
If row labels are included, there is an initial tab on the column label row.
}
\usage{
readAsGenericDataframe(theFile, theNaString=NULL, theUnknownString="Unknown")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theFile}{The cross-platform tab-delimited data.frame file with
  columns and optional row labels to be read. If row labels are included,
  there is an initial tab on the column label row.}
  \item{theNaString}{String to consider NA passed directly to read.csv in na.strings.}
  \item{theUnknownString}{String to replace empty strings with. Defaults to 'Unknown'.}
}
\value{
Returns the data.frame read from file.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
#readAsGenericDataframe <- function(theFile, theNaString=NULL, theUnknownString="Unknown")
#{
#	df <- read.csv(theFile, header=TRUE, sep="\t", as.is=TRUE, check.names=FALSE,
#	               stringsAsFactors=FALSE, colClasses="character", na.strings=theNaString,
#	               allowEscapes=TRUE, quote="" )
#	# check for colons in the batch type names (incompatible with MANOVA tests)
#	coln  <- colnames(df)
#	coln <- unlist(lapply(coln, function(theName)
#	{
#	  if (grepl(":", theName, fixed = TRUE))
#	  {
#	    logWarn("Replacing colon in batch type name ", theName, " with tilde")
#	    theName <- gsub(":", "~", theName, fixed = TRUE)
#	  }
#	  theName
#	}))
#	colnames(df) <- coln
#	df[df==""] <- theUnknownString
#  return(df)
#}
}
