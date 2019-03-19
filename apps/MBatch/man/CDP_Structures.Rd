\name{CDP_Structures}
\alias{CDP_Structures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CDP_Structures}
\description{Make a Correlation Density Plot for given data.}
\usage{
CDP_Structures(theFilePath, theData1, theData2, theSubTitle,
theUnmatchedCount = 1000, theMethod = "pearson",
theUse = "pairwise.complete.obs", theSeed = NULL,
theUseReplicatesUnpaired=FALSE,
theLinePlot=TRUE, theHistPlot=TRUE, theBinWidth=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theFilePath}{Full path and filename for PNG output file}
  \item{theData1}{Matrix with columns as samples and rows as features.}
  \item{theData2}{Matrix with columns as samples and rows as features.}
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (theFilePath, theData1, theData2, theSubTitle,
theUnmatchedCount = 1000,
    theMethod = "pearson", theUse = "pairwise.complete.obs",
    theSeed = NULL)
{
    theData1 <- theData1[rownames(theData1)[rownames(theData1) \%in\%
        rownames(theData2)], ]
    theData2 <- theData2[rownames(theData1)[rownames(theData1) \%in\%
        rownames(theData2)], ]
    pairedSamples <- colnames(theData1)[colnames(theData1) \%in\%
        colnames(theData2)]
    pairedDensity <- calculateDensityCorrFromList(theData1, theData2,
        pairedSamples, pairedSamples, theMethod, theUse, theSeed)
    data1samples <- colnames(theData1)[!colnames(theData1) \%in\%
        pairedSamples]
    data2samples <- colnames(theData2)[!colnames(theData2) \%in\%
        pairedSamples]
    unpairedSamples1 <- sample(data1samples, theUnmatchedCount,
        replace = TRUE)
    unpairedSamples2 <- sample(data2samples, theUnmatchedCount,
        replace = TRUE)
    unmatchedDensity <- calculateDensityCorrFromList(theData1,
        theData2, unpairedSamples1, unpairedSamples2, theMethod,
        theUse, theSeed)
    CairoPNG(filename = theFilePath, width = 480, height = 480)
    on.exit(dev.off())
    xRange <- c(min(c(pairedDensity$x, unmatchedDensity$x)),
        max(c(pairedDensity$x, unmatchedDensity$x)))
    plot(pairedDensity, type = "l", ann = FALSE, col = "deepskyblue",
        xlim = xRange, main =
        paste("Paired and Unmatched Correlation Density Plots\n",
            theSubTitle), xlab = convertCDPsubtitle(theMethod),
        ylab = "Density")
    par(new = TRUE)
    plot(unmatchedDensity, type = "l", ann = TRUE, col = "firebrick2",
        axes = FALSE, xlim = xRange, main =
        paste("Paired and Unmatched Correlation Density Plots\n",
            theSubTitle), xlab = convertCDPsubtitle(theMethod),
        ylab = "Density")
    legend("topleft", c(paste("paired N =",
    length(unique(sort(pairedSamples)))),
        paste("unmatched N =", length(unique(sort(c(data1samples,
            data2samples)))))), lty = 1, col = c("deepskyblue",
        "firebrick2"), bty = "n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
