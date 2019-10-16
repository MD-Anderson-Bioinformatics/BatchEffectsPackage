\name{EBNPlus_TrainAndValidateReplicates_Structures}
\alias{EBNPlus_TrainAndValidateReplicates_Structures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{EBNPlus_TrainAndValidateReplicates_Structures}
\description{Create a training and validation set and return
the corrected matrices for those data sets.}
\usage{
EBNPlus_TrainAndValidateReplicates_Structures(theDataMatrix1, theDataMatrix2,
theBatchId1, theBatchId2, theEBNP_BatchWithZero, theEBNP_FixDataSet,
theEBNP_CorrectForZero, theEBNP_ParametricPriorsFlag,
theEBNP_ValidationRatio, theEBNP_TestRatio=0, theSeed = NULL,
theEBNP_PriorPlotsFile=NULL, theEBNP_MinSampleNum = 3,
theEBNP_AddData1Rows = FALSE, theEBNP_AddData2Rows = FALSE,
theTestSeed=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theDataMatrix1}{A matrix for data set 1 containing numeric values
  with columns being sample ids and rows being gene ids.}
  \item{theDataMatrix2}{A matrix for data set 2 containing numeric values
  with columns being sample ids and rows being gene ids.}
  \item{theBatchId1}{The group name for data set 1 (for example, RNASeqV2)}
  \item{theBatchId2}{The group name for data set 2 (for example, Agilent4502)}
  \item{theEBNP_BatchWithZero}{Batch With Zero indicates which data set
  contains zero values. This is a string "1", "2", or "both".
  These values will be changed in a future release.}
  \item{theEBNP_FixDataSet}{Fix Data Set indicates which data set should
  be set as invariate during the corrections. Value are numeric, 1 or 2
  for the appropriate data set, or as.numeric(NA) for neither.}
  \item{theEBNP_CorrectForZero}{Correct For Zero indicates whether or
  not data sets marked
  as "Batch With Zero" should be modified from zero (TRUE or FALSE)}
  \item{theEBNP_ParametricPriorsFlag}{
  Use parametric adjustments for corrections (TRUE or FALSE)}
  \item{theEBNP_ValidationRatio}{A numeric. If 0, do not use a training set.
  For this function, the value should be less than 1. 0.33 means 33 percent
  of the replicates available are used
  for validation with 67 percent for training. }
  \item{theEBNP_TestRatio}{A numeric. If 0, do not use a test set.
  For this function, the value should be less than 1. 0.33 means 33 percent
  of the replicates available are used
  for the test set with 67 percent for training. }
  \item{theSeed}{NULL or a seed to use. Defaults to NULL.
  This is passed to the standard R set.seed function.}
  \item{theEBNP_PriorPlotsFile}{Defaults to NULL and does not create PNG.
  Non-NULL should be
  the full path and filename for where to create the prior plots PNG.}
  \item{theEBNP_MinSampleNum}{Defaults to 3. Any row (gene) with less than
  this number of samples is dropped.}
  \item{theEBNP_AddData1Rows}{Defaults to FALSE. TRUE indicates rows in
  data set 2 not in 1 should be added as all NA.}
  \item{theEBNP_AddData2Rows}{Defaults to FALSE. TRUE indicates rows in
  data set 1 not in 2 should be added as all NA.}
  \item{theTestSeed}{NULL or a seed to use. Defaults to NULL.
  This is passed to the standard R set.seed function.
  This is used only when creating Test Set.}
}
\details{
Take two matrices and determine training and validation data sets.
Using the EBNPlus algorithm and the provided arguments,
performing the training and validation corrections.
}
\value{
A list of matrices.
  \item{TestSet1}{a matrix of the test set from data set 1}
  \item{TestSet2}{a matrix of the test set from data set 2}
  \item{TrainingSet1}{a matrix of the training set from data set 1}
  \item{TrainingSet2}{a matrix of the training set from data set 2}
  \item{TrainingResults}{a matrix of the corrected training data}
  \item{ValidationSet1}{a matrix of the validation set from data set 1}
  \item{ValidationSet2}{a matrix of the validation set from data set 2}
  \item{ValidationResults}{a matrix of the corrected validation data}
  \item{CorrectedResults}{If Validation Ratio is 0, the validation set values
  will be 0 x 0 matrices, but this list member will contain a matrix of
  corrected values, using all replicates.}
}
