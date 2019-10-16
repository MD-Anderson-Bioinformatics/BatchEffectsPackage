\name{CDP_Files}
\alias{CDP_Files}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CDP_Files}
\description{Make a Correlation Density Plot for given data.}
\usage{
CDP_Files(theFilePath, theDataFile1, theDataFile2, theSubTitle,
theUnmatchedCount = 1000, theMethod = "pearson",
theUse = "pairwise.complete.obs", theSeed = NULL,
theUseReplicatesUnpaired=FALSE,
theLinePlot=TRUE, theHistPlot=TRUE, theBinWidth=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theFilePath}{Full path and filename for PNG output file}
  \item{theDataFile1}{Matrix file with columns as samples
  (with initial tab) and rows as features.}
  \item{theDataFile2}{Matrix file with columns as samples
  (with initial tab) and rows as features.}
  \item{theSubTitle}{Subtitle for image, giving data type being displayed.}
  \item{theUnmatchedCount}{Number of iterations for unpaired samples.}
  \item{theMethod}{Defaults to pearson.
  Valid values are: concordance, pearson, kendall, spearman.}
  \item{theUse}{Defaults to pairwise.complete.obs.
  Valid values are accepted by the method parameter to cor.}
  \item{theSeed}{Default to NULL.}
  \item{theUseReplicatesUnpaired}{Defaults to FALSE.
  If TRUE, use both the replicates and non-replicates for the unpaired plot.}
  \item{theLinePlot}{Default to TRUE. TRUE means plot the lines for Correlation Density Plots.}
  \item{theHistPlot}{Default to TRUE. TRUE means plot the histogram for Correlation Density Plots.}
  \item{theBinWidth}{Default to NULL. Non-null means to use the given wide for bins. Otherwise, use default for hist.}
}
\details{Make a Correlation Density Plot for given data.}
\value{No useful values}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (theFilePath, theDataFile1, theDataFile2, theSubTitle,
    theUnmatchedCount = 1000, theMethod = "pearson",
    theUse = "pairwise.complete.obs",
    theSeed = NULL)
{
    logInfo(paste("CDP_Files -- theDataFile1=", theDataFile1))
    myData1 <- readAsGenericMatrix(theDataFile1)
    logInfo(paste("CDP_Files -- theDataFile2=", theDataFile2))
    myData2 <- readAsGenericMatrix(theDataFile2)
    CDP_Structures(theFilePath, myData1, myData2, theSubTitle = theSubTitle,
        theUnmatchedCount = theUnmatchedCount, theMethod = theMethod,
        theUse = theUse, theSeed = theSeed)
  }
}
