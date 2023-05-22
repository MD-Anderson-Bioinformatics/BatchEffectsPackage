\name{AN_Unadjusted}
\alias{AN_Unadjusted}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{AN_Unadjusted}
\description{
Perform ANOVA Unadjusted corrections.
}
\usage{
AN_Unadjusted(theBeaData, theBatchType, thePath = NULL,
theDataVersion=NULL, theTestVersion=NULL,
theWriteToFile = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theBeaData}{An instance of BEA_DATA \code{\link{BEA_DATA-class}}}
  \item{theBatchType}{A string identifying the batch type to correct.}
  \item{thePath}{Output path for any files.}
	\item{theDataVersion}{Object of class \code{"character"} Version of the data prefaced
	with DATA_, such as DATA_2022-09-09-1600. Use empty string to not include. }
	\item{theTestVersion}{Object of class \code{"character"} Version of the test prefaced
	with TEST_, such as TEST_2022-10-10-1300. Use empty string to not include. }
  \item{theWriteToFile}{TRUE to write the corrected data to file and return
  the cleanFilePathname instead of the corrected matrix.}
}
\details{
Perform ANOVA unadjusted correction for the given batch type.
Returning the corrected matrix or a path to the filename.
}
\value{
Return the corrected matrix (columns are sample ids and rows are features)
or a path to the filename.
}
