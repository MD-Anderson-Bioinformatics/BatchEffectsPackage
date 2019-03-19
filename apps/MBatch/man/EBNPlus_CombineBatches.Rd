\name{EBNPlus_CombineBatches}
\alias{EBNPlus_CombineBatches}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{EBNPlus_CombineBatches}
\description{
Take two different data set's batch data.frames, combine them, and add an EBNPlus batch
giving original data sets.
}
\usage{
EBNPlus_CombineBatches(theBeaBatches1, theBeaBatches2, theEBNP_Data1BatchId,
theEBNP_Data2BatchId, theBarcodeTrimFunction = NULL, theSep=".")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theBeaBatches1}{The data.frame containing batch information for data set 1.
  The "Sample" column should contain barcodes and is required.}
  \item{theBeaBatches2}{The data.frame containing batch information for data set 2.
  The "Sample" column should contain barcodes and is required.}
  \item{theEBNP_Data1BatchId}{The Batch Id for data set 1, as passed to one of the
  other EBNPlus functions (for example, RNASeqV2)}
  \item{theEBNP_Data2BatchId}{The Batch Id for data set 2, as passed to one of the
  other EBNPlus functions (for example, Agilent4502)}
  \item{theBarcodeTrimFunction}{A function applied to trim barcodes if needed.
  This defaults to NULL (indicating no trimming) and should not be needed for TCGA data.}
  \item{theSep}{Separator used when adding ids to existing barcodes.
  This defaults to ".""}
}
\details{
Combine the two data.frames.
Add the new "batch names" at batches for the two data sets and the EBNPlus batch type.
}
\value{
Return a data.frame.
Names for the data.frame are original names, with a new column "EBNPlus".
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
