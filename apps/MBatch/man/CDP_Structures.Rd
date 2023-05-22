\name{CDP_Structures}
\alias{CDP_Structures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CDP_Structures}
\description{Make a Correlation Density Plot for given data.}
\usage{
CDP_Structures(thePath, theDataVersion, theTestVersion,
theFileName, theData1, theData2, theSubTitle,
theUnmatchedCount = 1000, theMethod = "pearson",
theUse = "pairwise.complete.obs", theSeed = NULL,
theUseReplicatesUnpaired=FALSE,
theLinePlot=TRUE, theHistPlot=TRUE, theBinWidth=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{thePath}{Full path for output}
	\item{theDataVersion}{Object of class \code{"character"} Version of the data prefaced
	with DATA_, such as DATA_2022-09-09-1600. Use empty string to not include. }
	\item{theTestVersion}{Object of class \code{"character"} Version of the test prefaced
	with TEST_, such as TEST_2022-10-10-1300. Use empty string to not include. }
  \item{theFileName}{PNG file name}
  \item{theData1}{Matrix with columns as samples and rows as features.}
  \item{theData2}{Matrix with columns as samples and rows as features.}
  \item{theSubTitle}{Subtitle for image, giving data type being displayed.}
  \item{theUnmatchedCount}{Number of iterations for unpaired samples.}
  \item{theMethod}{Defaults to pearson.
  Valid values are: pearson, kendall, spearman.}
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
