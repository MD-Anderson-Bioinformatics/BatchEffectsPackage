\name{MP_ByBatch}
\alias{MP_ByBatch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MP_ByBatch}
\description{
Perform Median Polish by Batch Corrections
}
\usage{
MP_ByBatch(theBeaData, theBatchType, thePath = NULL, theWriteToFile = FALSE)
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
Perform Median Polish correction by Batch.
Returning the corrected matrix or a path to the filename.
}
\value{
Return the corrected matrix (columns are sample ids and rows are features)
or a path to the filename.
}
