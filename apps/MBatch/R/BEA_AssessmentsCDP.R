# MBatch Copyright (c) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021 University of Texas MD Anderson Cancer Center
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
# MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>

calculateCDP_correllation <- function(theData1, theData2, theSamples1, theSamples2, theMethod, theUse)
{
  result <- NULL
  if (length(theSamples1)>0)
  {
    result <- sapply(1:length(theSamples1), function(theIndex)
    {
      myCorr <- NA
      if ("concordance"==theMethod)
      {
        myCorr <- epi.ccc(theData1[,theSamples1[theIndex]], theData2[,theSamples2[theIndex]])
        myCorr <- myCorr$rho.c$est
      }
      else
      {
        myCorr <- cor(theData1[,theSamples1[theIndex]], theData2[,theSamples2[theIndex]], method=theMethod, use=theUse)
      }
      myCorr
    })
  }
  result
}

calculateCDP_density <- function(theCorrelations)
{
  result <- NULL
  if (!is.null(theCorrelations))
  {
    result <- density(theCorrelations, na.rm=TRUE)
  }
  result
}

convertCDPsubtitle <- function(theName)
{
  if ("concordance"==theName)
  {
    theName <- "Concordance Correlation"
  }
  else if ("pearson"==theName)
  {
    theName <- "Pearson's Correlation"
  }
  else if ("kendall"==theName)
  {
    theName <- "Kendall's Correlation"
  }
  else if ("spearman"==theName)
  {
    theName <- "Spearman's Correlation"
  }
  theName
}

CDP_Plot <- function(theFilePath, theData1, theData2,
                     theData1PairedReplicates, theData2PairedReplicates,
                     theData1UnmatchedReplicates, theData2UnmatchedReplicates,
                     theSubTitle,
                     theMethod="pearson", theUse="pairwise.complete.obs", theSeed=NULL,
                     theLinePlot=TRUE, theHistPlot=TRUE, theBinWidth=NULL)
{
  stopifnotWithLogging("paired replicates must match in length", (length(theData1PairedReplicates)==length(theData2PairedReplicates)))
  stopifnotWithLogging("unpaired replicates must match in length", (length(theData1UnmatchedReplicates)==length(theData2UnmatchedReplicates)))
  logInfo("CDP_Plot theFilePath=", theFilePath)
  logInfo("CDP_Plot theData1PairedReplicates=", length(theData1PairedReplicates))
  #logInfo("CDP_Plot theData1PairedReplicates=", paste(theData1PairedReplicates, collapse=", ", sep=","))
  logInfo("CDP_Plot theData2PairedReplicates=", length(theData2PairedReplicates))
  #logInfo("CDP_Plot theData2PairedReplicates=", paste(theData2PairedReplicates, collapse=", ", sep=","))
  logInfo("CDP_Plot theData1UnmatchedReplicates=", length(theData1UnmatchedReplicates))
  logInfo("CDP_Plot theData2UnmatchedReplicates=", length(theData2UnmatchedReplicates))
  tryCatch({
    theSubTitle <- breakIntoTitle(theSubTitle)
    # make rows the same
    theData1 <- theData1[rownames(theData1)[rownames(theData1) %in% rownames(theData2)],]
    theData2 <- theData2[rownames(theData1)[rownames(theData1) %in% rownames(theData2)],]
    set.seed(theSeed)
    pairedCorr <- calculateCDP_correllation(theData1, theData2, theData1PairedReplicates, theData2PairedReplicates, theMethod, theUse)
    unmatchedCorr <- calculateCDP_correllation(theData1, theData2, theData1UnmatchedReplicates, theData2UnmatchedReplicates, theMethod, theUse)
    logInfo("CDP_Plot pairedCorr=", length(pairedCorr))
    #print(pairedCorr)
    logInfo("CDP_Plot unmatchedCorr=", length(unmatchedCorr))
    #print(unmatchedCorr)
    pairedHist <- NULL
    if(!is.null(pairedCorr))
    {
      if(!is.null(theBinWidth))
      {
        pairedHist <- hist(pairedCorr, plot=FALSE, breaks=seq(from=-1, to=1, by=theBinWidth))
      }
      else
      {
        pairedHist <- hist(pairedCorr, plot=FALSE)
      }
    }
    unmatchedHist <- NULL
    if(!is.null(unmatchedCorr))
    {
      if(!is.null(theBinWidth))
      {
        unmatchedHist <- hist(unmatchedCorr, plot=FALSE, breaks=seq(from=-1, to=1, by=theBinWidth))
      }
      else
      {
        unmatchedHist <- hist(unmatchedCorr, plot=FALSE)
      }
    }
    pairedDensity <- calculateCDP_density(pairedCorr)
    unmatchedDensity <- calculateCDP_density(unmatchedCorr)
    logInfo("CDP_Plot pairedDensity$x=", length(pairedDensity$x))
    logInfo("CDP_Plot pairedDensity$y=", length(pairedDensity$y))
    logInfo("CDP_Plot pairedDensity$bw=", pairedDensity$bw)
    logInfo("CDP_Plot unmatchedDensity$x=", length(unmatchedDensity$x))
    logInfo("CDP_Plot unmatchedDensity$y=", length(unmatchedDensity$y))
    logInfo("CDP_Plot unmatchedDensity$bw=", unmatchedDensity$bw)
    xRange <- c(max(-1, min(c(pairedDensity$x,unmatchedDensity$x))), min(1, max(c(pairedDensity$x,unmatchedDensity$x))))
    logInfo("CDP_Plot xRange=", xRange)
    #xRange <- c(-1, 1)
    # plot!
    pdf(NULL)
    CairoPNG(filename=theFilePath, width = 480, height = 480 )
    on.exit(dev.off())
    on.exit(par("oma"), add=TRUE)
    par(oma=c(3,0,0,0))
    myAnn <- TRUE
    myAxes <- TRUE
    if ((TRUE==theHistPlot)&&(!is.null(pairedHist)))
    {
      plot(pairedHist, ann=myAnn, col=rgb(0,0,1,1/4), border=rgb(0,0,1,1/4),  axes=myAxes, xlim=xRange,
           main=paste("Correlation Density Plots\n",theSubTitle),
           xlab=convertCDPsubtitle(theMethod), ylab="Density")
      par(new=TRUE)
      myAnn <- FALSE
      myAxes <- FALSE
    }
    if ((TRUE==theHistPlot)&&(!is.null(unmatchedHist)))
    {
      plot(unmatchedHist, ann=myAnn, col=rgb(1,0,0,1/4), border=rgb(1,0,0,1/4), axes=myAxes, xlim=xRange,
           main=paste("Correlation Density Plots\n",theSubTitle),
           xlab=convertCDPsubtitle(theMethod), ylab="Density")
      par(new=TRUE)
      myAnn <- FALSE
      myAxes <- FALSE
    }
    if ((TRUE==theLinePlot)&&(!is.null(pairedDensity)))
    {
      plot(pairedDensity, type='l', ann=myAnn, col=rgb(0,0,1), axes=myAxes, xlim=xRange,
           main=paste("Correlation Density Plots\n",theSubTitle),
           xlab=convertCDPsubtitle(theMethod), ylab="Density")
      par(new=TRUE)
      myAnn <- FALSE
      myAxes <- FALSE
    }
    if ((TRUE==theLinePlot)&&(!is.null(unmatchedDensity)))
    {
      plot(unmatchedDensity, type='l', ann=myAnn, col=rgb(1,0,0), axes=myAxes, xlim=xRange,
           main=paste("Correlation Density Plots\n",theSubTitle),
           xlab=convertCDPsubtitle(theMethod), ylab="Density")
      par(new=TRUE)
      myAnn <- FALSE
      myAxes <- FALSE
    }
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    nPaired <- 0
    if (is.null(theData1PairedReplicates))
    {
      nPaired <- 0
    }
    else if (is.null(theData2PairedReplicates))
    {
      nPaired <- 0
    }
    else
    {
      nPaired <- length(unique(sort(c(theData1PairedReplicates, theData2PairedReplicates))))
    }
    nUnmatched <- 0
    if (is.null(theData1UnmatchedReplicates))
    {
      nUnmatched <- 0
    }
    else if (is.null(theData2UnmatchedReplicates))
    {
      nUnmatched <- 0
    }
    else
    {
      nUnmatched <- length(unique(sort(c(theData1UnmatchedReplicates, theData2UnmatchedReplicates))))
    }
    legend("bottom",
           c(paste("Paired (n = ", nPaired, ")", sep="", collapse=""),
             paste("Unpaired (n = ", nUnmatched, ")", sep="", collapse="")),
           lty=1, col=c(rgb(0,0,1),rgb(1,0,0)), bty='n',
           xpd = TRUE)
  },
  warning=function(e)
  {
    logWarn("1 Unable to calculate CDP--too many NAs, Infinities or NaNs in data")
    unlink(theFilePath)
    openAndWriteIssuesLogFileCDP(dirname(theFilePath))
  },
  error=function(e)
  {
    logWarn("2 Unable to calculate CDP--too many NAs, Infinities or NaNs in data")
    unlink(theFilePath)
    openAndWriteIssuesLogFileCDP(dirname(theFilePath))
  })
}

openAndWriteIssuesLogFileCDP<-function(theOutputDir)
{
  myFile <- file(cleanFilePath(theOutputDir, "error.log"), "w+")
  on.exit(close(myFile))
  cat("Correlation not possible\n", file=myFile, append=TRUE)
}

CDP_Structures <- function(theFilePath, theData1, theData2, theSubTitle, theUnmatchedCount=1000,
                           theMethod="pearson", theUse="pairwise.complete.obs", theSeed=NULL, theUseReplicatesUnpaired=FALSE,
                           theLinePlot=TRUE, theHistPlot=TRUE, theBinWidth=NULL)
{
  # get list of natural replicates
  pairedSamples <- colnames(theData1)[colnames(theData1) %in% colnames(theData2)]
  # get list of unmatched replicates
  unpairedSamples1 <- NULL
  unpairedSamples2 <- NULL
  if (TRUE==theUseReplicatesUnpaired)
  {
    unpairedSamples1 <- sample(colnames(theData1), theUnmatchedCount, replace=TRUE)
    unpairedSamples2 <- sample(colnames(theData2), theUnmatchedCount, replace=TRUE)
  }
  else
  {
    data1samples <- colnames(theData1)[!colnames(theData1) %in% pairedSamples]
    data2samples <- colnames(theData2)[!colnames(theData2) %in% pairedSamples]
    if (is.null(pairedSamples))
    {
      data1samples <- colnames(theData1)
      data2samples <- colnames(theData2)
    }
    if ((0==length(data1samples))||(0==length(data2samples)))
    {
      unpairedSamples1 <- c()
      unpairedSamples2 <- c()
    }
    else
    {
      unpairedSamples1 <- sample(data1samples, theUnmatchedCount, replace=TRUE)
      unpairedSamples2 <- sample(data2samples, theUnmatchedCount, replace=TRUE)
    }
  }
  CDP_Plot(theFilePath, theData1, theData2, pairedSamples, pairedSamples, unpairedSamples1, unpairedSamples2,
           theSubTitle, theMethod, theUse, theSeed, theLinePlot=theLinePlot, theHistPlot=theHistPlot, theBinWidth=theBinWidth)
}

CDP_Files <- function(theFilePath, theDataFile1, theDataFile2, theSubTitle, theUnmatchedCount=1000,
                      theMethod="pearson", theUse="pairwise.complete.obs", theSeed=NULL, theUseReplicatesUnpaired=FALSE,
                      theLinePlot=TRUE, theHistPlot=TRUE, theBinWidth=NULL)
{
  logInfo(paste("CDP_Files -- theDataFile1=", theDataFile1))
  myData1 <- readAsGenericMatrix(theDataFile1)
  logInfo(paste("CDP_Files -- theDataFile2=", theDataFile2))
  myData2 <- readAsGenericMatrix(theDataFile2)
  CDP_Structures(theFilePath, myData1, myData2, theSubTitle=theSubTitle, theUnmatchedCount=theUnmatchedCount,
                 theMethod=theMethod, theUse=theUse, theSeed=theSeed, theUseReplicatesUnpaired=theUseReplicatesUnpaired,
                 theLinePlot=theLinePlot, theHistPlot=theHistPlot, theBinWidth=theBinWidth)
}
