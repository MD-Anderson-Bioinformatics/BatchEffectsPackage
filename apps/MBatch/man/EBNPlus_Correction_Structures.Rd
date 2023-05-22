\name{EBNPlus_Correction_Structures}
\alias{EBNPlus_Correction_Structures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{EBNPlus_Correction_Structures}
\description{
Take two matrices and correct them into a single data set using
the EBNPlus algorithm.
}
\usage{
EBNPlus_Correction_Structures(theDataMatrix1, theDataMatrix2,
theBatchId1, theBatchId2,
theEBNP_BatchWithZero, theEBNP_FixDataSet, theEBNP_CorrectForZero,
theEBNP_ParametricPriorsFlag,
theSeed = NULL, theOutputDir=NULL,
theDataVersion=NULL, theTestVersion=NULL,
thePriorFile=NULL, theEBNP_MinSampleNum = 3,
theEBNP_AddData1Rows = FALSE, theEBNP_AddData2Rows = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theDataMatrix1}{A matrix for data set 1 containing numeric values with
  columns being sample ids and rows being gene ids.}
  \item{theDataMatrix2}{A matrix for data set 2 containing numeric values with
  columns being sample ids and rows being gene ids.}
  \item{theBatchId1}{The group name for data set 1 (for example, RNASeqV2)}
  \item{theBatchId2}{The group name for data set 2 (for example, Agilent4502)}
  \item{theEBNP_BatchWithZero}{Batch With Zero indicates which data set
  contains zero values. This is a string "1", "2", or "both".
  These values will be changed in a future release.}
  \item{theEBNP_FixDataSet}{Fix Data Set indicates which data set should
  be set as invariate during the corrections. Value are numeric, 1 or 2 for
  the appropriate data set, or as.numeric(NA) for neither.}
  \item{theEBNP_CorrectForZero}{Correct For Zero indicates whether
  or not data sets
  marked as "Batch With Zero" should be modified from zero (TRUE or FALSE)}
  \item{theEBNP_ParametricPriorsFlag}{Use parametric adjustments
  for corrections (TRUE or FALSE)}
  \item{theSeed}{NULL or a seed to use. Defaults to NULL.
  This is passed to the standard R set.seed function.}
  \item{theOutputDir}{Defaults to NULL and does not create PNG.
  Non-NULL path for where to create the prior plots PNG.}
  \item{theDataVersion}{Object of class \code{"character"} Version of the data prefaced
	with DATA_, such as DATA_2022-09-09-1600. Use empty string to not include. }
	\item{theTestVersion}{Object of class \code{"character"} Version of the test prefaced
	with TEST_, such as TEST_2022-10-10-1300. Use empty string to not include. }
	\item{thePriorFile}{Defaults to NULL and does not create PNG.
  Filename for prior plots PNG.}
  \item{theEBNP_MinSampleNum}{Defaults to 3. Any row (gene) with less than
  this number of samples is dropped.}
  \item{theEBNP_AddData1Rows}{Defaults to FALSE. TRUE indicates rows in data
  set 2 not in 1 should be added as all NA.}
  \item{theEBNP_AddData2Rows}{Defaults to FALSE. TRUE indicates rows in data
  set 1 not in 2 should be added as all NA.}
}
\details{
Take two matrices and correct them into a single data set using the EBNPlus
algorithm and the provided arguments, using all replicates for training.

The function takes the variant data set and corrects it to become part of
the invariant data set using replicates shared between the sets. Replicates
are determined by finding matching sample ids. It uses a validation ratio
value of 0 to use all replicates for training. Any duplicates in genes or
barcodes are removed. Genes that start with "?" are also removed.
See EBNPlus_CheckData_Structures and EBNPlus_CombineBatches for
details on data prep.
}
\value{
BeaEBNplusFiles returns a matrix. The columns of the matrix are the
barcodes with "." and the group name attached to the barcodes. This extra
string is attached, since the correction algorithm requires replicates,
some barcode will be identical, so some method preserve distinctness is needed.
The row names are gene symbols.
}
