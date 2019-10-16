\name{SupervisedClustering_Batches_Structures}
\alias{SupervisedClustering_Batches_Structures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{SupervisedClustering_Batches_Structures}
\description{
A function to perform Supervised Clustering for batch types from data in files.
}
\usage{
SupervisedClustering_Batches_Structures(theData, theTitle, theOutputPath,
theDoHeatmapFlag,
theBatchTypeAndValuePairsToRemove = list(),
theBatchTypeAndValuePairsToKeep = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theData}{An instance of BEA_DATA \code{\link{BEA_DATA-class}}}
  \item{theTitle}{Object of class \code{"character"}
  Title to use in PNG files. }
  \item{theOutputPath}{Object of class \code{"character"}
  Directory in which to place output PNG files. }
  \item{theDoHeatmapFlag}{Object of class \code{"logical"}
  A flag indicating whether or not
  to create HC heatmap, where TRUE meants to create heatmap. }
  \item{theBatchTypeAndValuePairsToRemove}{Object of class \code{"list"}
  A list of vectors containing the batch type (or * for all types)
  and the value to remove. list() indicates
  none while NULL will cause an error. }
  \item{theBatchTypeAndValuePairsToKeep}{Object of class \code{"list"}
  A list of vectors containing the batch type (or * for all types) and a
  vector of the the value(s) to keep.
  list() indicates none while NULL will cause an error. }
}
\examples{
#	#############################################################################
#	mydir<-file.path(getwd()"working","SupervisedClustering_Batches_Structures")
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
#	SupervisedClustering_Batches_Structures(myData, myTitle, myOutputPath,
#list(), list(),
# theDoHeatmapFlag=TRUE)
}
