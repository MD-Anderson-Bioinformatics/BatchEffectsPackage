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
Use EBNPlus_Correction_Structures for Corrections

Use the EBNPlus_Correction_Structures function to performed 
corrections.

For most users, the function EBNPlus_Correction_Structures 
is what you want to use for processing. It is designed 
specifically to do corrections. The EBNPlus_TrainAndValidate 
functions are for researchers interested in the internal 
workings of the EBNPlus algorithm.

See the tests/EBNPlus_Correction_Structures.R or 
tests/EBNPlus_Correction_Files.R for details.

Combine the two data.frames.
Add the new "batch names" at batches for the two data sets and the EBNPlus batch type.
}
\value{
Return a data.frame.
Names for the data.frame are original names, with a new column "EBNPlus".
}
