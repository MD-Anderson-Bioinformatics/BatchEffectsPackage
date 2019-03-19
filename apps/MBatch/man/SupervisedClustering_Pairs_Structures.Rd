\name{SupervisedClustering_Pairs_Structures}
\alias{SupervisedClustering_Pairs_Structures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{SupervisedClustering_Pairs_Structures}
\description{
A function to perform Supervised Clustering for pairs of batch types
from data in files.
}
\usage{
SupervisedClustering_Pairs_Structures(theData, theTitle, theOutputPath,
theDoHeatmapFlag, theListOfBatchPairs,
theBatchTypeAndValuePairsToRemove=list(),
theBatchTypeAndValuePairsToKeep=list()
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theData}{An instance of BEA_DATA \code{\link{BEA_DATA-class}}}
  \item{theTitle}{Object of class \code{"character"}
  Title to use in PNG files. }
  \item{theOutputPath}{Object of class \code{"character"}
  Directory in which to place output PNG files. }
  \item{theDoHeatmapFlag}{Object of class \code{"logical"}
  A flag indicating whether or
  not to create HC heatmap, where TRUE meants to create heatmap. }
  \item{theListOfBatchPairs}{Object of class \code{"vector"}
  A list of strings, where pairs of strings give batch types
  to use for pairs assessment. Use c("") for empty set. }
  \item{theBatchTypeAndValuePairsToRemove}{Object of class \code{"list"}
  A list of vectors containing the batch type (or * for all types)
  and the value to remove.
  list() indicates none while NULL will cause an error.
  Default to list().}
  \item{theBatchTypeAndValuePairsToKeep}{Object of class \code{"list"}
  A list of vectors containing the batch type (or * for all types) and a
  vector of the the value(s) to keep.
  list() indicates none while NULL will cause an error.
  Default to list().}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#	#############################################################################
#	mydir <- file.path(getwd(),"working","SupervisedClustering_Pairs_Structures")
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
#	SupervisedClustering_Pairs_Structures(myData, myTitle, myOutputPath, list(),
#     list(),
#			theDoHeatmapFlag=TRUE,
#     theListOfBatchPairs=c("PlateId", "TSS", "BatchId", "TSS"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
