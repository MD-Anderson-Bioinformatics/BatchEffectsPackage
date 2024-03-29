\name{EBNPlus_CheckData_Structures}
\alias{EBNPlus_CheckData_Structures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{EBNPlus_CheckData_Structures}
\description{
Take two matrices and check that they will work as arguments to MBatch EBNPlus functions.
}
\usage{
EBNPlus_CheckData_Structures(theDataMatrix1, theDataMatrix2, theDataReplicates1 = NULL,
theDataReplicates2 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theDataMatrix1}{A matrix for data set 1 containing numeric values with columns being
  sample ids and rows being gene ids.}
  \item{theDataMatrix2}{A matrix for data set 2 containing numeric values with columns being
  sample ids and rows being gene ids.}
  \item{theDataReplicates1}{A vector of "replicates" in data set 1 used for corrections.
  Defaults to NULL. See "details" below.}
  \item{theDataReplicates2}{A vector of "replicates" in data set 2 used for corrections.
  Defaults to NULL. See "details" below.}
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

Performs a stop if conditions are not met.
This function checks the following:
Both matrix arguments pass is.matrix test.
Both matrices have column names.
Both matrices have row names and they interect at least once.
Both matrices have column names that intersect at least once or the replicate value vectors
are the same size and exist in the column names.
All data in the matrices is numeric.

Please note, support for performing corrections using replicate lists is not yet supported.
}
\value{
No meaningful value. Performs a stop if conditions are not met.
}
