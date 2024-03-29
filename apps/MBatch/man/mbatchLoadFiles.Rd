\name{mbatchLoadFiles}
\alias{mbatchLoadFiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mbatchLoadFiles}
\description{
Load the given files into a BEA_DATA object, sorting and checking the data for consistency.
}
\usage{
mbatchLoadFiles(theGeneDataFile, theBatchFile, theCovariateFile=NULL, theNaStrings=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theGeneDataFile}{Object of class \code{"character"} This is a string containing the full path to the gene data file which is tab delimited (without quotes) with sample ids across the top and gene-equivalents on the left. The sample ids start with a leading tab. }
  \item{theBatchFile}{Object of class \code{"character"} This is a string containing the full path to the batch data file which is tab delimited (without quotes). The first row gives the batch types. The first batch "type" is "Sample". }
  \item{theCovariateFile}{Object of class \code{"character"} This is a string containing the full path to the covariate data file which is tab delimited (without quotes). The first row gives the covariate types. The first covariate "type" is "Sample". }
  \item{theNaStrings}{This is the na.strings argument to read.csv. This is used because some data contains "NA" values that are code strings and not NA values. }
}
\value{
An instance of BEA_DATA \code{\link{BEA_DATA-class}}
}
\author{
Tod Casasent \email{tdcasasent@mdanderson.org}
}
