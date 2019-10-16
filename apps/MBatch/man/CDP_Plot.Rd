\name{CDP_Plot}
\alias{CDP_Plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CDP_Plot}
\description{Make a Correlation Density Plot for given data.}
\usage{
CDP_Plot(theFilePath, theData1, theData2,
theData1PairedReplicates, theData2PairedReplicates,
theData1UnmatchedReplicates, theData2UnmatchedReplicates,
theSubTitle, theMethod = "pearson",
theUse = "pairwise.complete.obs", theSeed = NULL,
theLinePlot=TRUE, theHistPlot=TRUE, theBinWidth=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theFilePath}{Full path and filename for PNG output file}
  \item{theData1}{Matrix with columns as samples and rows as features.}
  \item{theData2}{Matrix with columns as samples and rows as features.}
  \item{theData1PairedReplicates}{List of replicates to use
  for matched pairs from theData1.}
  \item{theData2PairedReplicates}{List of replicates to use
  for matched pairs from theData2.}
  \item{theData1UnmatchedReplicates}{List of replicates to use
  for unmatched pairs from theData1.}
  \item{theData2UnmatchedReplicates}{List of replicates to use
  for unmatched pairs from theData2.}
  \item{theSubTitle}{Subtitle for image, giving data type being displayed.}
  \item{theMethod}{Defaults to pearson.
  Valid values are: concordance, pearson, kendall, spearman.}
  \item{theUse}{Defaults to pairwise.complete.obs.
  Valid values are accepted by the method parameter to cor.}
  \item{theSeed}{Default to NULL.}
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
function (theFilePath, theData1, theData2, theData1PairedReplicates,
    theData2PairedReplicates,
    theData1UnmatchedReplicates, theData2UnmatchedReplicates,
    theSubTitle, theMethod = "pearson",
    theUse = "pairwise.complete.obs", theSeed = NULL)
{
    theData1 <- theData1[rownames(theData1)[rownames(theData1) \%in\%
        rownames(theData2)], ]
    theData2 <- theData2[rownames(theData1)[rownames(theData1) \%in\%
        rownames(theData2)], ]
    pairedDensity <- calculateDensityCorrFromList(theData1, theData2,
        theData1PairedReplicates, theData2PairedReplicates, theMethod,
        theUse, theSeed)
    unmatchedDensity <- calculateDensityCorrFromList(theData1,
        theData2, theData1UnmatchedReplicates, theData2UnmatchedReplicates,
        theMethod, theUse, theSeed)
    CairoPNG(filename = theFilePath, width = 480, height = 480)
    on.exit(dev.off())
    xRange <- c(min(c(pairedDensity$x, unmatchedDensity$x)),
        max(c(pairedDensity$x, unmatchedDensity$x)))
    plot(pairedDensity, type = "l", ann = FALSE, col = "deepskyblue",
        xlim = xRange,
        main = paste("Paired and Unmatched Correlation Density Plots\n",
            theSubTitle), xlab = convertCDPsubtitle(theMethod),
        ylab = "Density")
    par(new = TRUE)
    plot(unmatchedDensity, type = "l", ann = TRUE, col = "firebrick2",
        axes = FALSE, xlim = xRange,
        main = paste("Paired and Unmatched Correlation Density Plots\n",
            theSubTitle), xlab = convertCDPsubtitle(theMethod),
        ylab = "Density")
    legend("topleft", c(paste("paired (n =",
    length(unique(sort(c(theData1PairedReplicates,
        theData2PairedReplicates)))), ")"), paste("unpaired (n =",
        length(unique(sort(c(theData1UnmatchedReplicates,
        theData2UnmatchedReplicates)))),
        ")")), lty = 1, col = c("deepskyblue", "firebrick2"),
        bty = "n")
  }
}
