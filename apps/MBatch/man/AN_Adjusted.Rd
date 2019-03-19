\name{AN_Adjusted}
\alias{AN_Adjusted}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{AN_Adjusted}
\description{
Perform ANOVA Adjusted corrections.
}
\usage{
AN_Adjusted(theBeaData, theBatchType, thePath = NULL, theWriteToFile = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theBeaData}{An instance of BEA_DATA \code{\link{BEA_DATA-class}}}
  \item{theBatchType}{A string identifying the batch type to correct.}
  \item{thePath}{Output path for any files.}
  \item{theWriteToFile}{TRUE to write the corrected data to file and return
  the file pathname instead of the corrected matrix.}
}
\details{
Perform ANOVA adjusted correction for the given batch type.
Returning the corrected matrix or a path to the filename.
}
\value{
Return the corrected matrix (columns are sample ids and rows are features)
or a path to the filename.
}
\references{}
\author{}
\note{}
\seealso{}
\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
