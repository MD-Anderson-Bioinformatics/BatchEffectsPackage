\name{PCA_Regular_Structures}
\alias{PCA_Regular_Structures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{PCA_Regular_Structures}
\description{
A function to perform PCA analysis from data in files.
}
\usage{
PCA_Regular_Structures(theData, theTitle, theOutputPath,
theBatchTypeAndValuePairsToRemove, theBatchTypeAndValuePairsToKeep,
theIsPcaTrendFunction=function(...) {FALSE},
theDoCentroidsMtoMFlag=TRUE, theDoPlainMtoMFlag=FALSE,
theDoDSCFlag=TRUE, theDoDscPermsFileFlag=FALSE, theDoSampleLocatorFlag=TRUE,
theListOfComponentsToPlot=c(1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4),
theDSCPermutations=2000, theDSCThreads=1, theMinBatchSize=2,
theJavaParameters="-Xms1200m", theSeed=NULL, theMaxGeneCount=20000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theData}{An instance of BEA_DATA \code{\link{BEA_DATA-class}}}
  \item{theTitle}{Object of class \code{"character"}
  Title to use in PNG files. }
  \item{theOutputPath}{Object of class \code{"character"}
  Directory in which to place output PNG files. }
  \item{theBatchTypeAndValuePairsToRemove}{Object of class \code{"list"}
  A list of vectors containing the batch type (or * for all types)
  and the value to remove.
  list() indicates none while NULL will cause an error. }
  \item{theBatchTypeAndValuePairsToKeep}{Object of class \code{"list"}
  A list of vectors containing the batch type (or * for all types) and a
  vector of the
  value(s) to keep. list() indicates none while NULL will cause an error. }
  \item{theIsPcaTrendFunction}{Object of class \code{"function"}
  This function is used by the PCA Centroids algorithm to determine if
  trending lines should be added. If the function returns true, a line is
  drawn between the centroids, using the string sorted order. This means,
  for date or time fields, the batch ids should be formatted such that when
  sorted as strings, the batch ids are in chronological order from oldest
  to newest. The function receives arguments theBatchTypeName (a string
  identifying the batch type) and theListOfBatchIds (a list of strings,
  which are the batch ids). See below for an example. }
  \item{theDoCentroidsMtoMFlag}{Object of class \code{"logical"}
  A flag indicating whether or
  not to perform the many to many comparison with centroids. }
  \item{theDoPlainMtoMFlag}{Object of class \code{"logical"}
  A flag indicating whether or not
  to perform the many to many comparison without centroids. }
  \item{theDoDSCFlag}{Object of class \code{"logical"}
  A flag indicating whether or not to
  perform the DSC calculation for PCA-Plus plots. }
  \item{theDoDscPermsFileFlag}{Object of class \code{"logical"}
  A flag indicating whether or
  not to perform the DSC calculation for PCA-Plus plots. }
  \item{theDoSampleLocatorFlag}{Object of class \code{"logical"}
  A flag indicating whether or
  not to create the ANY_Scores.txt file, which contains the
  list of the PCA scores for samples versus components. }
  \item{theListOfComponentsToPlot}{Object of class \code{"vector"}
  A vector (list) of integers.
  Each pair of integers is used to create a PCA plot of those two components. }
  \item{theDSCPermutations}{Object of class \code{"numeric"}
  The number of permutations to
  perform for DSC calculations. }
  \item{theDSCThreads}{Object of class \code{"numeric"}
  The number of threads (or CPUs or cores)
  to use to perform DSC calculations. }
  \item{theMinBatchSize}{Object of class \code{"numeric"}
  Minimum batch to allow in PCA calculations. }
  \item{theJavaParameters}{Object of class \code{"character"}
  Command line arguments to
  Java 7 for DSC calculations--for example if memory needs to be
  more than "-Xms1200m". }
  \item{theSeed}{Object of class \code{"numeric"} If non-zero, a
  seed to use when calling
  sample to reduce genes used in PCA plots. }
  \item{theMaxGeneCount}{Object of class \code{"numeric"} If non-zero,
  maximum number of
  genes allowed for PCA plots. If count is greater than this,
  sample is used to select a subset of genes.}
}
\examples{
#	#############################################################################
#	mydir <- file.path(getwd(), "working", "PCA_Regular_Structures")
#	dir.create(mydir, recursive=TRUE)
#	setLogging(new("Logging", theFile=file.path(mydir, "mbatch.log")))
#	## path for output
#	myOutputPath <- file.path(mydir, "output")
#	dir.create(myOutputPath, recursive=TRUE)
#	#############################################################################
#	## this is the title used for the images
#	myTitle <- "Disease/Data Type/Platform/Data Level"
#	#############################################################################
#	# mbatchLoadFiles
#	#############################################################################
#	## This is the single gene data file to be read and processed
#	myGeneDataFile <- system.file("extdata", "matrix_data.tsv", package="MBatch")
#	## This is the single batch file to be read and processed
#	myBatchFile <- system.file("extdata", "batches.tsv", package="MBatch")
#	## this loads the files
#	myData <- mbatchLoadFiles(myGeneDataFile, myBatchFile)
#	###
#	isTrendBatch<-function(theBatchTypeName, theListOfBatchIds)
#	{
#		return(is.element(theBatchTypeName, c("ShipDate")))
#	}
#	###
#
#	PCA_Regular_Structures(myData, myTitle, myOutputPath, list(), list(),
# theIsPcaTrendFunction=isTrendBatch, theDoCentroidsMtoMFlag=TRUE,
# theDoPlainMtoMFlag=TRUE,
# theDoCentroidsOtoMFlag=TRUE, theDoPlainOtoMFlag=TRUE,
# theDoDSCFlag=TRUE, theDoDscPermsFileFlag=TRUE,
# theDoSampleLocatorFlag=TRUE,
#	theListOfComponentsToPlot=c(1, 2, 1, 3), theDSCPermutations=100,
#	theDSCThreads=5, theMinBatchSize=0,
#	theJavaParameters="-Xms1000m", theSeed=0, theMaxGeneCount=0)
}
