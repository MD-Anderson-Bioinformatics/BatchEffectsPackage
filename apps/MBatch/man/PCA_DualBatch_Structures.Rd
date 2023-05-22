\name{PCA_DualBatch_Structures}
\alias{PCA_DualBatch_Structures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{PCA_DualBatch_Structures}
\description{
A function to perform PCA dual batch analysis from data in files.
}
\usage{
PCA_DualBatch_Structures(theData, theTitle,
theBatchTypeAndValuePairsToRemove, theBatchTypeAndValuePairsToKeep,
theListForDoCentroidDualBatchType,
theOutputDir, theDataVersion, theTestVersion,
theIsPcaTrendFunction=function(...) {FALSE},
theDoDSCFlag=TRUE, theDoDscPermsFileFlag=FALSE, theDoSampleLocatorFlag=TRUE,
theListOfComponentsToPlot=c(1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4),
theDSCPermutations=2000, theDSCThreads=1, theMinBatchSize=2,
theSeed=NULL, theMaxGeneCount=20000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theData}{An instance of BEA_DATA \code{\link{BEA_DATA-class}}}
  \item{theTitle}{Object of class \code{"character"}
  Title to use in PNG files. }
  \item{theBatchTypeAndValuePairsToRemove}{Object of class \code{"list"}
  A list of vectors containing the batch type (or * for all types)
  and the value to remove.
  list() indicates none while NULL will cause an error. }
  \item{theBatchTypeAndValuePairsToKeep}{Object of class \code{"list"}
  A list of vectors containing the batch type (or * for all types)
  and a vector
  of the the value(s) to keep. list() indicates none while NULL
  will cause an error. }
  \item{theListForDoCentroidDualBatchType}{Object of class \code{"vector"}
  A vector (list) of strings. Each pair of strings should be batch types
  and are used to generate dual-batch plots. }
  \item{theOutputDir}{Object of class \code{"character"}
  Directory in which to place output PNG files. }
	\item{theDataVersion}{Object of class \code{"character"} Version of the data prefaced
	with DATA_, such as DATA_2022-09-09-1600. Use empty string to not include. }
	\item{theTestVersion}{Object of class \code{"character"} Version of the test prefaced
	with TEST_, such as TEST_2022-10-10-1300. Use empty string to not include. }
  \item{theIsPcaTrendFunction}{Object of class \code{"function"}
  This function is used by the PCA Centroids algorithm to determine if
  trending lines should be added. If the function returns true, a line
  is drawn between the centroids, using the string sorted order. This
  means, for date or time fields, the batch ids should be formatted such
  that when sorted as strings, the batch ids are in chronological order
  from oldest to newest. The function receives arguments
  theBatchTypeName (a string identifying the batch type) and
  theListOfBatchIds (a list of strings, which are the batch ids).
  See below for an example.
  Default to returning false.}
  \item{theDoDSCFlag}{Object of class \code{"logical"}
  A flag indicating whether or not to
  perform the DSC calculation for PCA-Plus plots. }
  \item{theDoDscPermsFileFlag}{Object of class \code{"logical"}
  A flag indicating whether
  or not to perform the DSC calculation for PCA-Plus plots. }
  \item{theDoSampleLocatorFlag}{Object of class \code{"logical"}
  A flag indicating whether or not to create the ANY_DSCPerms.txt and
  ANY_CompX_CompY_DSCPerms.txt files, which contains contain DSC, Db,
  Dw, a list of DSC, list of Dw, and list of Db for each permutation.
  File format consists of six lines per permutation, each line of the
  format "<perm number><tab><label><tab><value1><tab>...<tab><valueN>".
  The labels provided are DSC, Db, Dw, DSCList, DwList and DbList. }
  \item{theListOfComponentsToPlot}{Object of class \code{"vector"}
  A vector (list) of integers. Each pair of integers is used to
  create a PCA plot of those two components. Defaults to
  c(1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4)}
  \item{theDSCPermutations}{Object of class \code{"numeric"}
  The number of permutations to perform for DSC calculations.
  Defaults to 2000}
  \item{theDSCThreads}{Object of class \code{"numeric"}
  The number of threads (or CPUs or cores)
  to use to perform DSC calculations.
  Defaults to 1}
  \item{theMinBatchSize}{Object of class \code{"numeric"}
  Minimum batch to allow in PCA calculations. Defaults to 2.}
  \item{theSeed}{Object of class \code{"numeric"} If non-zero,
  a seed to use when
  calling sample to reduce genes used in PCA plots. Default 0.}
  \item{theMaxGeneCount}{Object of class \code{"numeric"} If non-zero,
  maximum number of genes allowed for PCA plots. If count is greater
  than this, sample is used to select a subset of genes. Default 20,000.}
}
