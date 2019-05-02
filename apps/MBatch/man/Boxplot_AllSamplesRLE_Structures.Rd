\name{Boxplot_AllSamplesRLE_Structures}
\alias{Boxplot_AllSamplesRLE_Structures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Boxplot_AllSamplesRLE_Structures}
\description{
A function to perform Box Plot RLE all data analysis from data in files.
}
\usage{
Boxplot_AllSamplesRLE_Structures(theData, theTitle, theOutputPath,
theBatchTypeAndValuePairsToRemove, theBatchTypeAndValuePairsToKeep,
theMaxGeneCount=20000, theJavaParameters = c("-Xms8000m", "-Djava.awt.headless=true"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theData}{An instance of BEA_DATA \code{\link{BEA_DATA-class}}}
  \item{theTitle}{Object of class \code{"character"} Title to use in PNG files. }
  \item{theOutputPath}{Object of class \code{"character"}
  Directory in which to place output PNG files. }
  \item{theBatchTypeAndValuePairsToRemove}{Object of class \code{"list"} A list of vectors
  containing the batch type (or * for all types) and the value to remove. list() indicates
  none while NULL will cause an error. }
  \item{theBatchTypeAndValuePairsToKeep}{Object of class \code{"list"} A list of vectors
  containing the batch type (or * for all types) and a vector of the the value(s) to keep.
  list() indicates none while NULL will cause an error. }
  \item{theMaxGeneCount}{Integer giving maximum number of features (genes) to keep. Default is 20000. 0 means keep all.}
  \item{theJavaParameters}{Object of class \code{"character"} String for initializing JVM.
  Defaults to Xms8000m and Djava.awt.headless=true. }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The All Samples RLE Boxplots plot the value for each feature (genes or probes) for a sample, with the samples grouped and colored by batch. In this case, RLE is used to move the mean of each sample to zero on the vertical axis. So the vertical axis is based on the values of the original data and the points plotted are features. The actual meaning of the data used, such as expression, read counts, and the like, will vary based on the data being processed.
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
#	mydir <- file.path(getwd(), "working", "Boxplot_AllSamplesRLE_Structures")
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
#	## This is the single batch file to be read and processed
#	myCovariateFile <- system.file("extdata", "covariates.tsv", package="MBatch")
#	## this loads the files
#	myData <- mbatchLoadFiles(myGeneDataFile, myBatchFile, myCovariateFile)
#	###
#
#	Boxplot_AllSamplesRLE_Structures(myData, myTitle, myOutputPath, list(), list())
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
