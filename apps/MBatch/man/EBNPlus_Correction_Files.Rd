\name{EBNPlus_Correction_Files}
\alias{EBNPlus_Correction_Files}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{EBNPlus_Correction_Files}
\description{
Take two Standardized TCGA Data matrix files and two batch files,
apply the normal TCGA Batch Effects Correction values to those files to
correct the data into one dataset.
}
\usage{
EBNPlus_Correction_Files(theDataFile1, theDataFile2,
theOutputDir, theDataVersion, theTestVersion,
theBatchId1, theBatchId2,
theSeed = NULL, theEBNP_PriorPlotsFlag=FALSE, theEBNP_MinSampleNum = 3,
theEBNP_AddData1Rows = FALSE, theEBNP_AddData2Rows = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theDataFile1}{A string containing the full path to a
  Standardized data matrix file for data set 1.}
  \item{theDataFile2}{A string containing the full path to a
  Standardized data matrix file for data set 2.}
  \item{theOutputDir}{Defaults to NULL and does not create PNG.
  Non-NULL path for where to create the prior plots PNG.}
  \item{theDataVersion}{Object of class \code{"character"} Version of the data prefaced
	with DATA_, such as DATA_2022-09-09-1600. Use empty string to not include. }
	\item{theTestVersion}{Object of class \code{"character"} Version of the test prefaced
	with TEST_, such as TEST_2022-10-10-1300. Use empty string to not include. }
  \item{theBatchId1}{The group name for data set 1 (for example, RNASeqV2)}
  \item{theBatchId2}{The group name for data set 2 (for example, Agilent4502)}
  \item{theSeed}{NULL or a seed to use. Defaults to NULL.
  This is passed to the standard R set.seed function.}
  \item{theEBNP_PriorPlotsFlag}{Defaults to FALSE.
  TRUE indicates prior plots PNG should be created.}
  \item{theEBNP_MinSampleNum}{
  Defaults to 3. Any row (gene) with less than this
  number of samples is dropped.}
  \item{theEBNP_AddData1Rows}{Defaults to FALSE.
  TRUE indicates rows in data set 2 not in 1 should be added as all NA.}
  \item{theEBNP_AddData2Rows}{Defaults to FALSE.
  TRUE indicates rows in data set 1 not in 2 should be added as all NA.}
}
\details{
EBNPlus_Correction_Files does corrections and creates files as output
for internal (MD Anderson) MBatch runs. While the batch files are not
used here, the more general functions wrapped by this
function require batch files.

The function EBNPlus_Correction_Files takes data set 2 and corrects it to
become part of data set 1 using replicates shared between the sets.
EBNPlus_Correction_Files uses batch with zero as 1, fix data set as NA,
and correct for zero as TRUE. It uses a validation ratio value of 0
to use all replicates for training. Barcodes are kept uncut. Genes
have pipes and any strings after the pipes removed. Any duplicates
in genes or barcodes are removed.
Genes that start with "?" are also removed.

Details on what the arguments mean can be found in
EBNPlus_Correction_Structures and
EBNPlus_TrainAndValidateReplicates_Structures.
}
\value{
EBNPlus_Correction_Files returns a list containing full path to the
matrix file with the corrected data set. The columns of the matrix are
the barcodes with "." and the group name attached to the barcodes.
This extra string is attached, since the correction algorithm requires
replicates, some barcode will be identical, so some method preserve
distinctness is needed. The row names are gene symbols.
}
