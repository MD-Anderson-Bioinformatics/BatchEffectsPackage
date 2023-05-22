\name{HierarchicalClustering_Structures}
\alias{HierarchicalClustering_Structures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{HierarchicalClustering_Structures}
\description{
A function to perform Hierarchical Clustering analysis from data in files.
}
\usage{
HierarchicalClustering_Structures(theData, theTitle, theOutputDir,
theDataVersion, theTestVersion,
theBatchTypeAndValuePairsToRemove, theBatchTypeAndValuePairsToKeep)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theData}{An instance of BEA_DATA \code{\link{BEA_DATA-class}}}
  \item{theTitle}{Object of class \code{"character"} Title to use in PNG files. }
  \item{theOutputDir}{Object of class \code{"character"}
  Directory in which to place output PNG files. }
  \item{theDataVersion}{Object of class \code{"character"} Version of the data prefaced
  with DATA_, such as DATA_2022-09-09-1600. Use empty string to not include. }
  \item{theTestVersion}{Object of class \code{"character"} Version of the test prefaced
  with TEST_, such as TEST_2022-10-10-1300. Use empty string to not include. }
  \item{theBatchTypeAndValuePairsToRemove}{Object of class \code{"list"} A list of vectors
  containing the batch type (or * for all types) and the value to remove. list() indicates
  none while NULL will cause an error. }
  \item{theBatchTypeAndValuePairsToKeep}{Object of class \code{"list"} A list of vectors
  containing the batch type (or * for all types) and a vector of the the value(s) to keep.
  list() indicates none while NULL will cause an error. }
}
\examples{
#	#############################################################################
#	mydir <- cleanFilePath(getwd(), "working", "HierarchicalClustering_Structures")
#	dir.create(mydir, recursive=TRUE)
#	setLogging(new("Logging", theFile=cleanFilePath(mydir, "mbatch.log")))
#	## path for output
#	myOutputPath <- cleanFilePath(mydir, "output")
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
#	## This is the single batch file to be read and processed
#	myCovariateFile <- system.file("extdata", "covariates.tsv", package="MBatch")
#	## this loads the files
#	myData <- mbatchLoadFiles(myGeneDataFile, myBatchFile, myCovariateFile)
#	###
#
#	HierarchicalClustering_Structures(myData, myTitle, myOutputPath,
# theDataVersion, theTestVersion, list(), list())
}
