\name{mbatchIncludeExcludeData}
\alias{mbatchIncludeExcludeData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mbatchIncludeExcludeData}
\description{Makes sure that the appropriate list(s) of included/excluded
samples/genes is enforced and returned the BEA Data object.}
\usage{
mbatchIncludeExcludeData(theBeaData, theIncludeSamples = NULL,
theIncludeGenes = NULL, theExcludeSamples = NULL, theExcludeGenes = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theBeaData}{An instance of BEA_DATA \code{\link{BEA_DATA-class}}}
  \item{theIncludeSamples}{A list of sample ids.
  Any sample ids not in this list are removed from the data.}
  \item{theIncludeGenes}{A list of gene ids.
  Any gene ids not in this list are removed from the data.}
  \item{theExcludeSamples}{A list of sample ids.
  Any sample ids in this list are removed from the data.}
  \item{theExcludeGenes}{A list of gene ids.
  Any gene ids in this list are removed from the data.}
}
\value{
Returns an instance of BEA_DATA \code{\link{BEA_DATA-class}}
with the updated data and batch information.
}
