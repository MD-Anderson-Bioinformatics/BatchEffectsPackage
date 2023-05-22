\name{Boxplot_AllSamplesData_Structures}
\alias{Boxplot_AllSamplesData_Structures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Boxplot_AllSamplesData_Structures}
\description{
A function to perform Box Plot all samples and all data analysis from data in Structures.
}
\usage{
Boxplot_AllSamplesData_Structures (theData, theTitle, theOutputDir,
theBatchTypeAndValuePairsToRemove, theBatchTypeAndValuePairsToKeep,
theDataVersion, theTestVersion,
theMaxGeneCount=20000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theData}{An instance of BEA_DATA \code{\link{BEA_DATA-class}}}
  \item{theTitle}{Object of class \code{"character"} Title to use in PNG files. }
  \item{theOutputDir}{Object of class \code{"character"} Directory in which to place output PNG
  files. }
  \item{theBatchTypeAndValuePairsToRemove}{Object of class \code{"list"} A list of vectors
  containing the batch type (or * for all types) and the value to remove. list() indicates
  none while NULL will cause an error. }
  \item{theBatchTypeAndValuePairsToKeep}{Object of class \code{"list"} A list of vectors
  containing the batch type (or * for all types) and a vector of the the value(s) to keep.
  list() indicates none while NULL will cause an error. }
	\item{theDataVersion}{Object of class \code{"character"} Version of the data prefaced
	with DATA_, such as DATA_2022-09-09-1600. Use empty string to not include. }
	\item{theTestVersion}{Object of class \code{"character"} Version of the test prefaced
	with TEST_, such as TEST_2022-10-10-1300. Use empty string to not include. }
  \item{theMaxGeneCount}{Integer giving maximum number of features (genes) to keep. Default is 20000. 0 means keep all.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The All Samples Boxplots plot the value for each feature (genes or probes) for a sample, with the samples grouped and colored by batch. So the vertical axis is based on the values of the original data and the points plotted are features. The actual meaning of the data used, such as expression, read counts, and the like, will vary based on the data being processed.
}
\examples{
#	#############################################################################
#	mydir <- cleanFilePath(getwd(), "working", "Boxplot_AllSamplesData_Structures")
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
#	Boxplot_AllSamplesData_Structures(myData, myTitle, myOutputPath, list(), list())
}
