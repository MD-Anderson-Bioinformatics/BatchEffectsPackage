# MBatch Copyright (c) 2011-2022 University of Texas MD Anderson Cancer Center
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
# MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>


boxplot_openAndWriteIssuesLogFile<-function(theOutputDir)
{
  myFile <- file(cleanFilePath(theOutputDir, "error.log"), "w+")
  on.exit(close(myFile))
  cat("Unable to calculate boxplot\n", file=myFile, append=TRUE)
}

###########################################################################################
###########################################################################################

createBatchEffectsOutput_BoxPlot_AllSampleRLE<-function(theMatrixGeneData, theDataframeBatchData,
																												theTitle, theOutputDir,
																												theDataVersion, theTestVersion, thePngFlag)
{
  logDebug("createBatchEffectsOutput_BoxPlot_AllSampleRLE - theOutputDir", theOutputDir)
  checkDirForCreation(theOutputDir)
  logDebug("dim(theMatrixGeneData)", dim(theMatrixGeneData))
  logDebug("length(colnames(theMatrixGeneData))", length(colnames(theMatrixGeneData)))
  logDebug("length(rownames(theMatrixGeneData))", length(rownames(theMatrixGeneData)))
  logDebug("dim(theDataframeBatchData)", dim(theDataframeBatchData))
  logDebug("length(names(theDataframeBatchData))", length(names(theDataframeBatchData)))
  theOutputDir <- cleanFilePath(theOutputDir, "AllSample-RLE")
  theOutputDir <- addVersionsIfNeeded(theOutputDir, theDataVersion, theTestVersion)
  for(batchTypeIndex in c(2:length(theDataframeBatchData)))
  {
    ### compile data and information for display
    batchTypeName <- names(theDataframeBatchData)[batchTypeIndex]
    logDebug("batchTypeName = ", batchTypeName)
    #	  samplesForBatches <- as.character(as.vector(unlist(theDataframeBatchData$Sample)))
    #	  batchIdsForSamples <- as.character(as.vector(unlist(theDataframeBatchData[batchTypeIndex])))
    # add " / batch-id" to sample ids for theMatrixGeneData
    # index of sorted values from theDataframeBatchData
    #	  labelsSorted <- theDataframe$Sample[order(batchIdsForSamples)]
    #	  batchIdsSorted <- batchIdsForSamples[order(batchIdsForSamples)]
    #	  labelWithBatch <- paste(labelsSorted, batchIdsSorted, sep=" / ")
    #
    # build names for output files
    theBoxDataFile <- cleanFilePath(theOutputDir, paste("BoxPlot_", "AllSample-RLE", "_BoxData-", batchTypeName, ".tsv", sep=""))
    theHistogramFile <- cleanFilePath(theOutputDir, paste("BoxPlot_", "AllSample-RLE", "_Histogram-", batchTypeName, ".tsv", sep=""))
    theAnnotationsFile <- cleanFilePath(theOutputDir, paste("BoxPlot_", "AllSample-RLE", "_Annotations-", batchTypeName, ".tsv", sep=""))
    thePngFile <- cleanFilePath(theOutputDir, paste("BoxPlot_", "AllSample-RLE", "_Diagram-", batchTypeName, ".png", sep=""))
    title <- paste(theTitle, "/ Boxplot AllSample-RLE", batchTypeName, sep=" ")
    boxplotLabels <- buildBoxplotLabels(theDataframeBatchData, batchTypeName)
    calcAndWriteBoxplot(theMatrixGeneData, title, theBoxDataFile, thePngFile, theHistogramFile, theAnnotationsFile, boxplotLabels)
  }
  # let this return the output dir, used by MBatchUtils
  theOutputDir
}

buildBoxplotLabels <- function(theDataframeBatchData, theBatchType)
{
  logDebug("theBatchType=", theBatchType)
  samplesForBatches <- as.character(as.vector(unlist(theDataframeBatchData$Sample)))
  batchIdsForSamples <- as.character(as.vector(unlist(theDataframeBatchData[theBatchType])))
  # add " / batch-id" to sample ids for theMatrixGeneData
  # index of sorted values from theDataframeBatchData
  labelWithBatch <- paste(batchIdsForSamples, samplesForBatches, sep=" / ")
  labelWithBatch
}

createBatchEffectsOutput_BoxPlot_AllSampleData<-function(theMatrixGeneData, theDataframeBatchData,
																							 theTitle, theOutputDir, theDataVersion, theTestVersion, thePngFlag)
{
  checkDirForCreation(theOutputDir)
	logDebug("dim(theMatrixGeneData)", dim(theMatrixGeneData))
	logDebug("length(colnames(theMatrixGeneData))", length(colnames(theMatrixGeneData)))
	logDebug("length(rownames(theMatrixGeneData))", length(rownames(theMatrixGeneData)))
	logDebug("dim(theDataframeBatchData)", dim(theDataframeBatchData))
	logDebug("length(names(theDataframeBatchData))", length(names(theDataframeBatchData)))
	theOutputDir <- cleanFilePath(theOutputDir, "AllSample-Data")
	theOutputDir <- addVersionsIfNeeded(theOutputDir, theDataVersion, theTestVersion)
	for(batchTypeIndex in c(2:length(theDataframeBatchData)))
	{
	  ### compile data and information for display
	  batchTypeName <- names(theDataframeBatchData)[batchTypeIndex]
	  logDebug("batchTypeName = ", batchTypeName)
	  # build names for output files
	  theBoxDataFile <- cleanFilePath(theOutputDir, paste("BoxPlot_", "AllSample-Data", "_BoxData-", batchTypeName, ".tsv", sep=""))
	  theHistogramFile <- cleanFilePath(theOutputDir, paste("BoxPlot_", "AllSample-Data", "_Histogram-", batchTypeName, ".tsv", sep=""))
	  theAnnotationsFile <- cleanFilePath(theOutputDir, paste("BoxPlot_", "AllSample-Data", "_Annotations-", batchTypeName, ".tsv", sep=""))
	  thePngFile <- cleanFilePath(theOutputDir, paste("BoxPlot_", "AllSample-Data", "_Diagram-", batchTypeName, ".png", sep=""))
	  title <- paste(theTitle, "/ Boxplot AllSample-Data", batchTypeName, sep=" ")
	  boxplotLabels <- buildBoxplotLabels(theDataframeBatchData, batchTypeName)
	  #
	  calcAndWriteBoxplot(theMatrixGeneData, title, theBoxDataFile, thePngFile, theHistogramFile, theAnnotationsFile, boxplotLabels)
	}
	# let this return the output dir, used by MBatchUtils
	theOutputDir
}

createBatchEffectsOutput_BoxPlot_Group<-function(theMatrixGeneData, theDataframeBatchData,
																								 theTitle, theOutputDir,
																								 theListOfGroupBoxFunction, theListOfGroupBoxLabels,
																								 theDataVersion, theTestVersion,
																								 thePngFlag)
{
	stopifnotWithLogging("Group box list and group box labels should be the same length", length(theListOfGroupBoxFunction)==length(theListOfGroupBoxLabels))
  checkDirForCreation(theOutputDir)
	logDebug("groupFunction")
	logDebug("dim(theMatrixGeneData)", dim(theMatrixGeneData))
	logDebug("length(colnames(theMatrixGeneData))", length(colnames(theMatrixGeneData)))
	logDebug("length(rownames(theMatrixGeneData))", length(rownames(theMatrixGeneData)))
	logDebug("dim(theDataframeBatchData)", dim(theDataframeBatchData))
	logDebug("length(names(theDataframeBatchData))", length(names(theDataframeBatchData)))
	vectorOfOutputDirs <- c()
	for(index in 1:length(theListOfGroupBoxFunction))
	{
	  groupLabel <- theListOfGroupBoxLabels[[index]]
	  logDebug("groupLabel = ", groupLabel)
	  groupAndLabel <- paste("Group-", groupLabel, sep="")
	  myOutputPath <- cleanFilePath(theOutputDir, groupAndLabel)
	  myOutputPath <- addVersionsIfNeeded(myOutputPath, theDataVersion, theTestVersion)
	  vectorOfOutputDirs <- c(vectorOfOutputDirs, myOutputPath)
	  for(batchTypeIndex in c(2:length(theDataframeBatchData)))
	  {
	    ### compile data and information for display
	    batchTypeName <- names(theDataframeBatchData)[batchTypeIndex]
	    logDebug("batchTypeName = ", batchTypeName)
	    ###############################################################
	    # build new matrix, with
	    # columns being batches
	    # rows being sample ids
	    # sample ids not in group are set to NA
	    ###############################################################
	    newMatrix <- buildGroupMatrix(theMatrixGeneData, theDataframeBatchData, batchTypeName, theListOfGroupBoxFunction[[index]])
	    ####
	    # build names for output files
	    theBoxDataFile <- cleanFilePath(myOutputPath, paste("BoxPlot_", groupAndLabel, "_BoxData-", batchTypeName, ".tsv", sep=""))
	    theHistogramFile <- cleanFilePath(myOutputPath, paste("BoxPlot_", groupAndLabel, "_Histogram-", batchTypeName, ".tsv", sep=""))
	    theAnnotationsFile <- cleanFilePath(myOutputPath, paste("BoxPlot_", groupAndLabel, "_Annotations-", batchTypeName, ".tsv", sep=""))
	    thePngFile <- cleanFilePath(myOutputPath, paste("BoxPlot_", groupAndLabel, "_Diagram-", batchTypeName, ".png", sep=""))
	    title <- paste(theTitle, "/ Boxplot", groupAndLabel, batchTypeName, sep=" ")
	    boxplotLabels <- colnames(newMatrix)
	    #
	    calcAndWriteBoxplot(newMatrix, title, theBoxDataFile, thePngFile, theHistogramFile, theAnnotationsFile, boxplotLabels)
	  }
	}
	# let this return vector of output dirs, used by MBatchUtils
	vectorOfOutputDirs
}

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

calcAndWriteBoxplot <- function(theData, theTitle, theBoxDataFile, thePngFile, theHistogramFile,
                                theAnnotationsFile, theBoxplotLabels)
{
  #logDebug("theBoxplotLabels=", theBoxplotLabels)
  logDebug("calcAndWriteBoxplot - theBoxDataFile=", theBoxDataFile)
  logDebug("calcAndWriteBoxplot - dim(theData)[1]=", dim(theData)[1])
  logDebug("calcAndWriteBoxplot - dim(theData)[2]=", dim(theData)[2])
  checkDirForCreation(dirname(theBoxDataFile))
  tryCatch({
    checkIfTestError()
    writeTitleFile(theTitle, thePngFile)
    logDebug("calcAndWriteBoxplot - before calcAndWriteBoxDataFile")
    calcAndWriteBoxDataFile(theData, theBoxDataFile, thePngFile, theBoxplotLabels, theTitle)
    logDebug("calcAndWriteBoxplot - after calcAndWriteBoxDataFile")
    logDebug("calcAndWriteBoxplot - before calcAndWriteHistogramFile")
    calcAndWriteHistogramFile(theData, theHistogramFile)
    logDebug("calcAndWriteBoxplot - after calcAndWriteHistogramFile")
    logDebug("calcAndWriteBoxplot - before calcAndWriteAnnotationsFile")
    calcAndWriteAnnotationsFile(theData, theAnnotationsFile)
    logDebug("calcAndWriteBoxplot - after calcAndWriteAnnotationsFile")
  },error=function(myError)
  {
    logError("calcAndWriteBoxplot(): error in calcAndWriteAnnotationsFile, ERROR= ", myError)
    boxplot_openAndWriteIssuesLogFile(dirname(theBoxDataFile))
  })
}

###########################################################################################
###########################################################################################

getBatchFactorsForBoxplot <- function(theBoxplotLabels)
{
  splitted <- strsplit(x=theBoxplotLabels, split=" / ", fixed=TRUE)
  batches <- sapply(splitted, function(x) x[1])
  factor(batches)
}

getColorsForBoxplot <- function(theBoxplotLabels)
{
  batchLevels <- getBatchFactorsForBoxplot(theBoxplotLabels)
  batchCount <- length(unique(sort(batchLevels)))
  colors <- beaRainbow(batchCount)
  colors[batchLevels]
}

calcAndWriteBoxDataFile <- function(theData, theFile, thePngFile, theBoxplotLabels, theTitle)
{
  logDebug("calcAndWriteBoxDataFile theFile=", theFile)
  logDebug("calcAndWriteBoxDataFile thePngFile=", thePngFile)
  listOfFivenumResults <- list()
  ######################################
  # theBoxDataFile Id	LowerOutMax	LowerOutMin	LowerNotch	LowerWhisker	LowerHinge	Median	UpperHinge	UpperWhisker	UpperNotch	UpperOutMin	UpperOutMax
  cat("Id	LowerOutMax	LowerOutMin	LowerNotch	LowerWhisker	LowerHinge	Median	UpperHinge	UpperWhisker	UpperNotch	UpperOutMin	UpperOutMax", "\n", file=theFile, sep="", append=FALSE)
  for(i in 1:ncol(theData))
  {
    # do not subtract median of entire matrix
    # boxNumbers <- fivenum(as.vector(unlist(theData[,i]-theMedian)), na.rm=TRUE)
    boxNumbers <- fivenum(as.vector(unlist(theData[,i])), na.rm=TRUE)
    listOfFivenumResults[[length(listOfFivenumResults)+1]] <- boxNumbers
    cat(paste(
      colnames(theData)[i], # Id
      minEpsilon(theData[,i], boxNumbers[1]),     # LowerOutMax
      NA,                   #	LowerOutMin
      NA,                   # LowerNotch
      paste(boxNumbers, collapse="\t"), # LowerWhisker	LowerHinge	Median	UpperHinge	UpperWhisker
      NA,                   #	UpperNotch
      NA,                   # UpperOutMin
      maxEpsilon(theData[,i], boxNumbers[5]),     # UpperOutMax
      sep="\t"), "\n", sep="", file=theFile, append=TRUE)
  }
  propWidth = 20*length(theBoxplotLabels)
  if (propWidth<500)
  {
    propWidth <-500
  }
  if (propWidth>2000)
  {
    propWidth <-2000
  }
  logDebug("calcAndWriteBoxDataFile CairoPNG=", thePngFile)
  CairoPNG(filename=thePngFile, pointsize=12, width = propWidth, height = 600)
  on.exit(dev.off())
  on.exit(par("mar"), add=TRUE)
  par(mar=c(20,4,4,2))
  # , theBoxplotLabels
  listOfFivenumResults <- listOfFivenumResults[order(theBoxplotLabels)]
  logDebug("calcAndWriteBoxDataFile call boxplot")
  boxplot(listOfFivenumResults, main=theTitle,
          col=getColorsForBoxplot(sort(theBoxplotLabels)),
          xaxt="n", xlab="")
  logDebug("calcAndWriteBoxDataFile call text")
  text(1:length(listOfFivenumResults), par("usr")[3] - 0.05, labels=sort(theBoxplotLabels), xpd=TRUE, srt=45, adj=1)
  logDebug("calcAndWriteBoxDataFile done")
}

###########################################################################################
###########################################################################################

calcAndWriteHistogramFile <- function(theData, theFile)
{
  logDebug("calcAndWriteHistogramFile ", theFile)
  # count number of breaks
  maxBins <- 0
  #logDebug("calcAndWriteHistogramFile 1 ncol(theData) =", ncol(theData))
  #logDebug("calcAndWriteHistogramFile 1 maxBins =", maxBins)
  for(myCheckIndex in 1:ncol(theData))
  {
    #logDebug("calcAndWriteHistogramFile 1 hist myCheckIndex =", myCheckIndex)
    # do not subtract median of entire matrix
    # myH <- hist(theData[,myCheckIndex]-theMedian, plot=FALSE)
    myH <- hist(theData[,myCheckIndex], plot=FALSE)
    tmpBins <- length(myH$mids)
    # use gt, not lt, to get max bin size
    if (tmpBins>maxBins)
    {
      maxBins <- tmpBins
    }
    #logDebug("calcAndWriteHistogramFile 1 tmpBins =", tmpBins)
    #logDebug("calcAndWriteHistogramFile 1b maxBins =", maxBins)
  }
  #logDebug("calcAndWriteHistogramFile ncol(theData)=", ncol(theData))
  #logDebug("calcAndWriteHistogramFile maxBins=", maxBins)
  ######################################
  # theHistogramFile entry	size	x0	y0	x1	y1  ...	xN	yN
  cat("entry	size", file=theFile, sep="", append=FALSE)
  for (myXYindex in 0:(maxBins-1))
  {
    cat(paste("\tx", myXYindex, "\ty", myXYindex, sep=""), file=theFile, sep="", append=TRUE)
  }
  cat("\n", file=theFile, sep="", append=TRUE)
  # and then the values
  for(myHistIndex in 1:ncol(theData))
  {

    #logDebug("calcAndWriteHistogramFile 2 hist myHistIndex =", myHistIndex)
    # do not subtract median of entire matrix
    #myH <- hist(theData[,i]-theMedian, plot=FALSE)
    myH <- hist(theData[,myHistIndex], plot=FALSE)
    cat(paste(
      colnames(theData)[myHistIndex], # entry
      length(myH$counts),      # size
      paste(myH$density, myH$mids, sep="\t", collapse="\t"),
      sep="\t"), "\n", sep="", file=theFile, append=TRUE)
  }
}

###########################################################################################
###########################################################################################

calcAndWriteAnnotationsFile <- function(theData, theFile)
{
  logDebug("calcAndWriteAnnotationsFile theFile=", theFile)
  ######################################
  # theAnnotationsFile
  # key	value
  # Total-Data-Points	5000
  # Non-NA-Points-TCGA-06-0650-01A-02D-1697-05	4980
  cat("key	value", "\n", file=theFile, sep="", append=FALSE)
  cat(paste("Total-Data-Points", ncol(theData),sep="\t"), "\n", file=theFile, sep="", append=TRUE)
  for(i in 1:ncol(theData))
  {
    cat(paste(paste("Non-NA-Points-", colnames(theData)[i], sep=""), sum(is.na(theData[,i])), sep="\t"), "\n", file=theFile, sep="", append=TRUE)
  }
}

###########################################################################################
###########################################################################################

buildGroupMatrix <- function(theDataMatrix, theBatchDF, theBatchType, theFunction)
{
  logDebug("buildGroupMatrix - dim(theDataMatrix)[1]=", dim(theDataMatrix)[1])
  logDebug("buildGroupMatrix - dim(theDataMatrix)[2]=", dim(theDataMatrix)[2])
  logDebug("buildGroupMatrix - dim(theBatchDF)[1]=", dim(theBatchDF)[1])
  logDebug("buildGroupMatrix - dim(theBatchDF)[2]=", dim(theBatchDF)[2])
  logDebug("buildGroupMatrix theBatchType=", theBatchType)
  ###############################################################
  # build new matrix, with
  # columns being batches
  # rows being sample ids
  # sample ids not in group are set to NA
  ###############################################################
  fullBatchList <- as.vector(unlist(theBatchDF[theBatchType]))
  logDebug("buildGroupMatrix fullBatchList=", length(fullBatchList))
  batchNameList <- unique(sort(fullBatchList))
  logDebug("buildGroupMatrix batchNameList=", length(batchNameList))
  sampleIdList <- theBatchDF$Sample
  logDebug("buildGroupMatrix sampleIdList=", length(sampleIdList))
  # build matrix of NAs
  newMatrix <- matrix(nrow=length(sampleIdList), ncol=length(batchNameList), dimnames=list(sampleIdList, batchNameList))
  for (batchName in batchNameList)
  {
    for (sampleId in sampleIdList[which(fullBatchList==batchName)])
    {
      val <- theDataMatrix[,sampleId]
      res <- theFunction(val)
      newMatrix[sampleId, batchName] <- res
    }
  }
  newMatrix
}

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
