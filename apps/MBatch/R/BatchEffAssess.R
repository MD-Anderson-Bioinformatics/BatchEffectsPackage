#MBatch Copyright ? 2011, 2012, 2013, 2014, 2015, 2016, 2017 University of Texas MD Anderson Cancer Center
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(gtools, warn.conflicts=FALSE, verbose=FALSE)
library(rJava, warn.conflicts=FALSE, verbose=FALSE)

getMBatchVersion<-function()
{
	return("MBatch Version: 2019-09-04-1100")
}

mbatchWriteSuccessfulLog <- function()
{
  logInfo("mbatchAssess Finishing")
}

sortMatrix <- function(theMatrix)
{
  if (!is.null(theMatrix))
  {
    if (dim(theMatrix)[1]>0)
    {
      if (dim(theMatrix)[2]>0)
      {
        # order rows, columns
        theMatrix <- theMatrix[order(rownames(theMatrix)), order(colnames(theMatrix))]
      }
    }
  }
  theMatrix
}

sortDataframe <- function(theDataframe)
{
  if (!is.null(theDataframe))
  {
    if (dim(theDataframe)[1]>0)
    {
      if (dim(theDataframe)[2]>0)
      {
        theDataframe <- theDataframe[order(theDataframe$Sample),]
      }
    }
  }
  theDataframe
}

mbatchLoadStructures<-function(theGeneMatrix, theBatchDataframe, theCovariatedataframe=NULL)
{
	########################################################
	collateOrigValue<-Sys.getlocale("LC_COLLATE")
	on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
	Sys.setlocale("LC_COLLATE","C")
	logDebug("Changing LC_COLLATE to C for duration of run")
	on.exit(warnings(), add=TRUE)
	### set options to display warnings after the script finishes.
	options(warn=0) # warnings shown after script finishes
	### if there is a warning, show the calls leading up to it
	options(showWarnCalls=TRUE)
	### if there is an error, show the calls leading up to it
	options(showErrorCalls=TRUE)
	########################################################
	logInfo("\\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ ")
	logInfo("Starting mbatchLoadStructures")
	logInfo(getMBatchVersion())
	#######################################################################
	logInfo("sort batches by gene file samples")
	batSamples <- theBatchDataframe$Sample
	genSamples <- colnames(theGeneMatrix)
	theBatchDataframe <- theBatchDataframe[order(batSamples, genSamples),]
	logInfo("filter samples in batches using gene samples")
	theBatchDataframe <- filterSamplesToBatchesDataframe(theBatchDataframe, theGeneMatrix)
	logInfo("sort covariates by gene file samples")
	if(!is.null(theCovariatedataframe))
	{
		logInfo("filter samples in covariates using gene samples")
		theCovariatedataframe <- filterSamplesToBatchesDataframe(theCovariatedataframe, theGeneMatrix)
		covSamples <- theCovariatedataframe$Sample
		genSamples <- colnames(theGeneMatrix)
		theCovariatedataframe <- theCovariatedataframe[order(covSamples, genSamples),]
	}
	else
	{
		theCovariatedataframe <- data.frame()
	}
	myData<-new("BEA_DATA", sortMatrix(theGeneMatrix), sortDataframe(theBatchDataframe), sortDataframe(theCovariatedataframe))
	logInfo("Finishing mbatchLoadStructures")
	logInfo("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
	return(myData)
}

mbatchLoadFiles<-function(theGeneDataFile, theBatchFile, theCovariateFile=NULL, theNaStrings=NULL)
{
	########################################################
	collateOrigValue<-Sys.getlocale("LC_COLLATE")
	on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
	Sys.setlocale("LC_COLLATE","C")
	logDebug("Changing LC_COLLATE to C for duration of run")
	on.exit(warnings(), add=TRUE)
	### set options to display warnings after the script finishes.
	options(warn=0) # warnings shown after script finishes
	### if there is a warning, show the calls leading up to it
	options(showWarnCalls=TRUE)
	### if there is an error, show the calls leading up to it
	options(showErrorCalls=TRUE)
	########################################################
	logInfo("\\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ ")
	logInfo("Starting mbatchLoadFiles")
	logInfo(getMBatchVersion())
	#######################################################################
	### load the data
	#######################################################################
	### load the gene data file into matrix
	### colnames(matrixGeneData) gives sample ids
	### rownames(matrixGeneData) gives gene symbols
	matrixGeneData <- NULL
	dataframeSamplesToBatches <- NULL
	dataframeSamplesToCovariates <- data.frame()
	logInfo("read batch file=", theBatchFile)
	#dataframeSamplesToBatches <- readBatchFile(theBatchFile)
	dataframeSamplesToBatches <- readAsGenericDataframe(theBatchFile, theNaStrings)
	logInfo("read gene file=", theGeneDataFile)
	#matrixGeneData <- readGeneFile(theGeneDataFile)
	matrixGeneData <- readAsGenericMatrix(theGeneDataFile)
	#logDebug("mbatchLoadFiles samples from batch = ", paste(dataframeSamplesToBatches$Sample, collapse=", "))
	#logDebug("mbatchLoadFiles samples from genes = ", paste(colnames(matrixGeneData), collapse=", "))
	logInfo("filter samples in batches using gene samples")
	dataframeSamplesToBatches <- filterSamplesToBatchesDataframe(dataframeSamplesToBatches, matrixGeneData)
	#logDebug("mbatchLoadFiles samples from batch = ", paste(dataframeSamplesToBatches$Sample, collapse=", "))
	#logDebug("mbatchLoadFiles samples from genes = ", paste(colnames(matrixGeneData), collapse=", "))
	logInfo("sort batches by gene file samples")
	batSamples <- dataframeSamplesToBatches$Sample
	genSamples <- colnames(matrixGeneData)
	dataframeSamplesToBatches <- dataframeSamplesToBatches[order(batSamples, genSamples),]
	#logDebug("mbatchLoadFiles samples from batch = ", paste(dataframeSamplesToBatches$Sample, collapse=", "))
	#logDebug("mbatchLoadFiles samples from genes = ", paste(colnames(matrixGeneData), collapse=", "))
	if(!is.null(theCovariateFile))
	{
		logInfo("read covariate file")
		#dataframeSamplesToCovariates <- readBatchFile(theCovariateFile)
		dataframeSamplesToCovariates <- readAsGenericDataframe(theCovariateFile)
		logInfo("filter samples in covariates using gene samples")
		dataframeSamplesToCovariates <- filterSamplesToBatchesDataframe(dataframeSamplesToCovariates, matrixGeneData)
		logInfo("sort covariates by gene file samples")
		genSamples <- colnames(matrixGeneData)
		covSamples <- dataframeSamplesToCovariates$Sample
		dataframeSamplesToCovariates <- dataframeSamplesToCovariates[order(covSamples, genSamples),]
	}
	#logDebug("mbatchLoadFiles samples from batch = ", paste(dataframeSamplesToBatches$Sample, collapse=", "))
	#logDebug("mbatchLoadFiles samples from genes = ", paste(colnames(matrixGeneData), collapse=", "))
	myData<-new("BEA_DATA", sortMatrix(matrixGeneData), sortDataframe(dataframeSamplesToBatches), sortDataframe(dataframeSamplesToCovariates))
	logInfo("Finishing mbatchLoadFiles")
	logInfo("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
	return(myData)
}

#c("Type"),
mbatchFilterData<-function(theBeaData,
													 theBatchTypeAndValuePairsToRemove=list(c("*", "unknown"), c("*", "Unknown")),
													 theBatchTypeAndValuePairsToKeep=list(list("Type", c("01", "03", "05"))),
													 theBatchTypesToRemove=NULL,
													 theMinIqr=0, theMinSd=0, theMinMad=0)
{
	########################################################
	collateOrigValue<-Sys.getlocale("LC_COLLATE")
	on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
	Sys.setlocale("LC_COLLATE","C")
	logDebug("Changing LC_COLLATE to C for duration of run")
	on.exit(warnings(), add=TRUE)
	### set options to display warnings after the script finishes.
	options(warn=0) # warnings shown after script finishes
	### if there is a warning, show the calls leading up to it
	options(showWarnCalls=TRUE)
	### if there is an error, show the calls leading up to it
	options(showErrorCalls=TRUE)
	########################################################
	logInfo("\\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ ")
	logInfo("mbatchFilterData Starting")
	logInfo(getMBatchVersion())
	matrixFilteredGeneData <- filterGeneDataMatrix(theBeaData@mData, theMinIqr, theMinSd, theMinMad)
	logDebug("mbatchFilterData Prefilter, gene data had ", nrow(theBeaData@mData), " while post filter ", nrow(matrixFilteredGeneData))
	stopifnotWithLogging("mbatchFilterData Filtered gene data contains no data to process.",nrow(matrixFilteredGeneData)>0)
	### load the batch data into matrix with (unfiltered) samples
	### dataframeSamplesToBatches$Sample gives list of samples
	### as.vector(unlist(names(dataframeSamplesToBatches)[2:length(dataframeSamplesToBatches)])) gives list of batch types
	### batchIdsForSamples <- as.vector(unlist(dataframeSamplesToBatches[batchTypeIndex])) gives batch ids for all samples
	### sortedListOfBatchIds <- sort(unique(batchIdsForSamples)) gives unique set of batch ids for current samples
	#logDebug("mbatchFilterData dataframe theBeaData@mBatches$Sample = ", paste(theBeaData@mBatches$Sample, collapse=", "))
	dataframeFilteredSamplesToBatches <- filterSamplesToBatchesDataframe(theBeaData@mBatches, matrixFilteredGeneData)
	#logDebug("mbatchFilterData dataframe $Sample = ", paste(dataframeFilteredSamplesToBatches$Sample, collapse=", "))
	logDebug("mbatchFilterData Prefilter, batch data had ", nrow(theBeaData@mBatches), " while post filter ", nrow(dataframeFilteredSamplesToBatches))
	stopifnotWithLogging("mbatchFilterData Filtered sample batches contain no data to process.", nrow(dataframeFilteredSamplesToBatches)>0)
	#####################################################################################################################
	### remove any samples with given batch types and values
	if (length(theBatchTypeAndValuePairsToRemove)>0)
	{
		logDebug("mbatchFilterData Before removing requested batches, batch data has ", length(dataframeFilteredSamplesToBatches$Sample), " samples")
		#logDebug("mbatchFilterData dataframe $Sample = ", paste(dataframeFilteredSamplesToBatches$Sample, collapse=", "))
		logDebug("mbatchFilterData Before removing requested batches, gene data has ", length(colnames(matrixFilteredGeneData)), " samples")
		logDebug("mbatchFilterData Remove requested batches (", paste(theBatchTypeAndValuePairsToRemove, sep="", collapse=", "), ")")
		dataframeFilteredSamplesToBatches <- removeSamplesWithListedBatches(dataframeFilteredSamplesToBatches, theBatchTypeAndValuePairsToRemove)
		logDebug("mbatchFilterData After removing requested batches, batch data has ", length(dataframeFilteredSamplesToBatches$Sample), " samples")
		#logDebug("mbatchFilterData dataframe $Sample = ", paste(dataframeFilteredSamplesToBatches$Sample, collapse=", "))
		matrixFilteredGeneData <- setGenesToSameAsSamples(matrixFilteredGeneData, dataframeFilteredSamplesToBatches)
		logDebug("mbatchFilterData After removing requested batches, gene data has ", length(colnames(matrixFilteredGeneData)), " samples")
	}
	#####################################################################################################################
	### only keep samples with given batch types and values
	if (length(theBatchTypeAndValuePairsToKeep)>0)
	{
		logDebug("mbatchFilterData Before keeping requested batches, batch data has ", length(dataframeFilteredSamplesToBatches$Sample), " samples")
		#logDebug("mbatchFilterData dataframe $Sample = ", paste(dataframeFilteredSamplesToBatches$Sample, collapse=", "))
		logDebug("mbatchFilterData Before keeping requested batches, gene data has ", length(colnames(matrixFilteredGeneData)), " samples")
		logDebug("mbatchFilterData Remove keeping batches (", paste(theBatchTypeAndValuePairsToKeep, sep="", collapse=", "), ")")
		dataframeFilteredSamplesToBatches <- keepSamplesWithListedBatches(dataframeFilteredSamplesToBatches, theBatchTypeAndValuePairsToKeep)
		logDebug("mbatchFilterData After keeping requested batches, batch data has ", length(dataframeFilteredSamplesToBatches$Sample), " samples")
		#logDebug("mbatchFilterData dataframe $Sample = ", paste(dataframeFilteredSamplesToBatches$Sample, collapse=", "))
		matrixFilteredGeneData <- setGenesToSameAsSamples(matrixFilteredGeneData, dataframeFilteredSamplesToBatches)
		logDebug("mbatchFilterData After keeping requested batches, gene data has ", length(colnames(matrixFilteredGeneData)), " samples")
	}
	#####################################################################################################################
	### remove the given batch types
	if (FALSE==is.null(theBatchTypesToRemove))
	{
		logDebug("mbatchFilterData Before removing batch types, batch data has ", length(dataframeFilteredSamplesToBatches$Sample), " samples")
		#logDebug("mbatchFilterData dataframe $Sample = ", paste(dataframeFilteredSamplesToBatches$Sample, collapse=", "))
		logDebug("mbatchFilterData Before removing batch types, gene data has ", length(colnames(matrixFilteredGeneData)), " samples")
		logDebug("mbatchFilterData removing batch types (", paste(theBatchTypesToRemove, sep="", collapse=", "), ")")
		dataframeFilteredSamplesToBatches <- removeListedBatchTypes(dataframeFilteredSamplesToBatches, theBatchTypesToRemove)
		logDebug("mbatchFilterData After removing batch types, batch data has ", length(dataframeFilteredSamplesToBatches$Sample), " samples")
		#logDebug("mbatchFilterData dataframe $Sample = ", paste(dataframeFilteredSamplesToBatches$Sample, collapse=", "))
		matrixFilteredGeneData <- setGenesToSameAsSamples(matrixFilteredGeneData, dataframeFilteredSamplesToBatches)
		logDebug("mbatchFilterData After removing batch types, gene data has ", length(colnames(matrixFilteredGeneData)), " samples")
	}
	#####################################################################################################################
	theBeaData@mData <- matrixFilteredGeneData
	theBeaData@mBatches <- dataframeFilteredSamplesToBatches
	logInfo("mbatchFilterData Finishing")
	logInfo("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
	return(theBeaData)
}

mbatchTrimData<-function(theMatrix, theMaxSize=15000000)
{
	########################################################
	collateOrigValue<-Sys.getlocale("LC_COLLATE")
	on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
	Sys.setlocale("LC_COLLATE","C")
	logDebug("Changing LC_COLLATE to C for duration of run")
	on.exit(warnings(), add=TRUE)
	### set options to display warnings after the script finishes.
	options(warn=0) # warnings shown after script finishes
	### if there is a warning, show the calls leading up to it
	options(showWarnCalls=TRUE)
	### if there is an error, show the calls leading up to it
	options(showErrorCalls=TRUE)
	########################################################
	logInfo("\\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ ")
	logInfo("mbatchTrimData Starting")
	logInfo(getMBatchVersion())
	if ((theMaxSize>0) && ((nrow(theMatrix)*ncol(theMatrix))>theMaxSize))
	{
		numberOfGenesToUse = floor(theMaxSize/ncol(theMatrix))
		iqr<-apply(theMatrix, 1, IQR, na.rm=TRUE)
		iqr <- sort(iqr, decreasing=TRUE, na.last=TRUE)
		genesToKeep <- sort(names(iqr)[1:numberOfGenesToUse])
		theMatrix <- theMatrix[genesToKeep, ]
	}
	logInfo("mbatchTrimData theMaxSize=", theMaxSize)
	logInfo("mbatchTrimData ncol(theMatrix)=", ncol(theMatrix))
	logInfo("mbatchTrimData nrow(theMatrix)=", nrow(theMatrix))
	logInfo("mbatchTrimData Finishing")
	logInfo("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
	return(theMatrix)
}

mbatchIncludeExcludeData<-function(theBeaData,
																	 theIncludeSamples=NULL, theIncludeGenes=NULL,
																	 theExcludeSamples=NULL, theExcludeGenes=NULL)
{
	########################################################
	collateOrigValue<-Sys.getlocale("LC_COLLATE")
	on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
	Sys.setlocale("LC_COLLATE","C")
	logDebug("Changing LC_COLLATE to C for duration of run")
	on.exit(warnings(), add=TRUE)
	### set options to display warnings after the script finishes.
	options(warn=0) # warnings shown after script finishes
	### if there is a warning, show the calls leading up to it
	options(showWarnCalls=TRUE)
	### if there is an error, show the calls leading up to it
	options(showErrorCalls=TRUE)
	########################################################
	logInfo("\\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ \\/ ")
	logInfo("mbatchIncludeExcludeData Starting")
	logInfo(getMBatchVersion())
	#logInfo("mbatchIncludeExcludeData theIncludeSamples ", theIncludeSamples)
	#logInfo("mbatchIncludeExcludeData theIncludeGenes ", theIncludeGenes)
	#logInfo("mbatchIncludeExcludeData theExcludeSamples ", theExcludeSamples)
	#logInfo("mbatchIncludeExcludeData theExcludeGenes ", theExcludeGenes)
	### colnames(matrixGeneData) gives sample ids
	### rownames(matrixGeneData) gives gene symbols
	if (!is.null(theIncludeSamples))
	{
		logInfo("mbatchIncludeExcludeData include only these samples")
		# include only these samples
		matrixFilteredGeneData <- theBeaData@mData[,colnames(theBeaData@mData)%in%theIncludeSamples]
		dataframeFilteredSamplesToBatches <- filterSamplesToBatchesDataframe(theBeaData@mBatches, matrixFilteredGeneData)
		theBeaData@mData <- matrixFilteredGeneData
		theBeaData@mBatches <- dataframeFilteredSamplesToBatches
	}
	if (!is.null(theIncludeGenes))
	{
		logInfo("mbatchIncludeExcludeData include only these genes")
		# include only these genes
		matrixFilteredGeneData <- theBeaData@mData[rownames(theBeaData@mData)%in%theIncludeGenes,]
		dataframeFilteredSamplesToBatches <- filterSamplesToBatchesDataframe(theBeaData@mBatches, matrixFilteredGeneData)
		theBeaData@mData <- matrixFilteredGeneData
		theBeaData@mBatches <- dataframeFilteredSamplesToBatches
	}
	if (!is.null(theExcludeSamples))
	{
		logInfo("mbatchIncludeExcludeData exclude these samples")
		# exclude these samples
		matrixFilteredGeneData <- theBeaData@mData[,!(colnames(theBeaData@mData)%in%theExcludeSamples)]
		dataframeFilteredSamplesToBatches <- filterSamplesToBatchesDataframe(theBeaData@mBatches, matrixFilteredGeneData)
		theBeaData@mData <- matrixFilteredGeneData
		theBeaData@mBatches <- dataframeFilteredSamplesToBatches
	}
	if (!is.null(theExcludeGenes))
	{
		logInfo("mbatchIncludeExcludeData exclude these genes")
		# exclude these genes
		matrixFilteredGeneData <- theBeaData@mData[!(rownames(theBeaData@mData)%in%theExcludeGenes),]
		dataframeFilteredSamplesToBatches <- filterSamplesToBatchesDataframe(theBeaData@mBatches, matrixFilteredGeneData)
		theBeaData@mData <- matrixFilteredGeneData
		theBeaData@mBatches <- dataframeFilteredSamplesToBatches
	}
	logInfo("mbatchIncludeExcludeData Finishing")
	logInfo("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
	return(theBeaData)
}

################################################################################

compressIntoFilename<-function(theString)
{
	### listing whole list of characters out looks wrong, but is locale independent
	theString <- gsub("[^ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789/\\]", "", theString)
	theString <- gsub("\\", "_", theString, fixed=TRUE)
	theString <- gsub("/", "_", theString, fixed=TRUE)
	return(theString)
}

###breakIntoTitle<-function(theString, theOldChar="\\", theNewChar=" ", theWidth=40)
breakIntoTitle<-function(theString, theOldChar=.Platform$file.sep, theNewChar=" ", theWidth=40)
{
  if (nchar(theString)>theWidth)
  {
    theString <- gsub(theOldChar, theNewChar, theString, fixed=TRUE)
  	theString <- strwrap(theString, width=theWidth)
  	theString <- paste(theString, collapse="\n")
  	theString <- gsub(theNewChar, theOldChar, theString, fixed=TRUE)
  }
	return(theString)
}

beaRainbow<-function(theNumberOfColors, v=1.0, s=1.0, shuffle=FALSE)
{
	### v=1.0 - for HierClust
	### s=0.7 - for BatchCorr
	### v=0.7 - for PCA
	colors <- rainbow(theNumberOfColors, v=v, s=s)
	if ((length(colors)>4)&&(shuffle==TRUE))
	{
		colorsAB <- colors[even(1:theNumberOfColors)]
		colorsA <- colorsAB[even(1:length(colorsAB))]
		colorsB <- colorsAB[odd(1:length(colorsAB))]
		colorsCD <- colors[odd(1:theNumberOfColors)]
		colorsC <- colorsCD[even(1:length(colorsCD))]
		colorsD <- colorsCD[odd(1:length(colorsCD))]
		colors <- c(colorsA, colorsC, colorsB, colorsD)
	}
	return(colors)
}

################################################################################

matrixWithIssues<-function(...)
{
	warnLevel<-getOption("warn")
	on.exit(options(warn=warnLevel))
	options(warn=3) # warnings are errors
	return(matrix(...))
}

asmatrixWithIssues<-function(...)
{
	warnLevel<-getOption("warn")
	on.exit(options(warn=warnLevel))
	options(warn=3) # warnings are errors
	return(as.matrix(...))
}

################################################################################

filterGeneDataMatrix<-function(theMatrixGeneData, theMinIqr, theMinSd, theMinMad)
{
	logDebug("rows pre filter ", nrow(theMatrixGeneData))
  ### clear NaN and INF
  theMatrixGeneData[is.nan(theMatrixGeneData)] <- NA
  theMatrixGeneData[is.infinite(theMatrixGeneData)] <- NA
  ### clear NAs
  #na.gene<-as.vector(unlist(apply(is.na(theMatrixGeneData),1,all)))
  #geneData<-theMatrixGeneData[!na.gene,]
  naLengene<-as.vector(unlist(apply(!is.na(theMatrixGeneData),1,function(row)
  {
    return (sum(row)>1)
  })))
  geneData<-theMatrixGeneData[naLengene,]
  ### do IQR, SD and MAD
  iqr<-as.vector(unlist(apply(geneData,1, IQR, na.rm=TRUE)))
  Sd<-as.vector(unlist(apply(geneData,1, sd, na.rm=TRUE)))
  Mad<-as.vector(unlist(apply(geneData, 1, mad, na.rm=TRUE)))
  temp<-geneData[(iqr>=theMinIqr) & (Sd>=theMinSd) & (Mad>=theMinMad),]
  logDebug("rows post filter ", nrow(temp))
	return(temp)
}

removeListedBatchTypes<-function(theDataframeFilteredSamplesToBatches, theBatchTypesToRemove)
{
	return(theDataframeFilteredSamplesToBatches[,!(names(theDataframeFilteredSamplesToBatches) %in% theBatchTypesToRemove), drop=FALSE])
}

removeSamplesWithListedBatches<-function(theDataframeFilteredSamplesToBatches, theBatchTypeAndValuePairsToRemove)
{
	for(myPair in theBatchTypeAndValuePairsToRemove)
	{
		theDataframeFilteredSamplesToBatches <- removeBatchAndValue(theDataframeFilteredSamplesToBatches, myPair[[1]], myPair[[2]])
	}
	return(theDataframeFilteredSamplesToBatches)
}

keepSamplesWithListedBatches<-function(theDataframeFilteredSamplesToBatches, theBatchTypeAndValuePairsToKeep)
{
	for(myPair in theBatchTypeAndValuePairsToKeep)
	{
		theDataframeFilteredSamplesToBatches <- keepBatchAndValues(theDataframeFilteredSamplesToBatches, myPair[[1]], myPair[[2]])
	}
	return(theDataframeFilteredSamplesToBatches)
}

removeBatchAndValue<-function(theDataframeFilteredSamplesToBatches, theBatch, theValue)
{
	if ("*"==theBatch)
	{
		for(batch in colnames(theDataframeFilteredSamplesToBatches)[2:length(theDataframeFilteredSamplesToBatches)])
		{
			theDataframeFilteredSamplesToBatches <- removeBatchAndValue(theDataframeFilteredSamplesToBatches, batch, theValue)
		}
	}
	else
	{
	  if (sum(theDataframeFilteredSamplesToBatches[,theBatch]!=theValue)>0)
	  {
	    theDataframeFilteredSamplesToBatches <- theDataframeFilteredSamplesToBatches[theDataframeFilteredSamplesToBatches[,theBatch]!=theValue,]
	  }
	}
	return(theDataframeFilteredSamplesToBatches)
}

keepBatchAndValues<-function(theDataframeFilteredSamplesToBatches, theBatch, theValues)
{
	if ("*"==theBatch)
	{
		for(batch in colnames(theDataframeFilteredSamplesToBatches)[2:length(theDataframeFilteredSamplesToBatches)])
		{
			theDataframeFilteredSamplesToBatches <- keepBatchAndValues(theDataframeFilteredSamplesToBatches, batch, theValues)
		}
	}
	else
	{
		theDataframeFilteredSamplesToBatches <- theDataframeFilteredSamplesToBatches[theDataframeFilteredSamplesToBatches[,theBatch]%in%theValues,]
	}
	return(theDataframeFilteredSamplesToBatches)
}

# old code to filter all values given
#tfList<-apply(theDataframeFilteredSamplesToBatches, 1, function(xlist)
#{
#	xlist<-as.vector(xlist)
#	xlist<-xlist[2:length(xlist)]
#	overlap<-is.element(xlist, theListOfBatchesToRemove)
#	return(sum(overlap)==0)
#})
#return(theDataframeFilteredSamplesToBatches[tfList,])


### load the batch data into matrix with (unfiltered) samples
### dataframeSamplesToBatches$Sample gives list of samples
### as.vector(unlist(names(dataframeSamplesToBatches)[2:length(dataframeSamplesToBatches)])) gives list of batch types
### batchIdsForSamples <- as.vector(unlist(dataframeSamplesToBatches[batchTypeIndex])) gives batch ids for all samples
### sortedListOfBatchIds <- sort(unique(batchIdsForSamples)) gives unique set of batch ids for current samples

### load the gene data file into matrix
### colnames(matrixGeneData) gives sample ids
### rownames(matrixGeneData) gives gene symbols

setGenesToSameAsSamples<-function(theMatrixFilteredGeneData, theDataframeFilteredSamplesToBatches)
{
	sampleidList <- theDataframeFilteredSamplesToBatches$Sample
	#logDebug("setGenesToSameAsSamples ($Sample from data frame) sampleidList=", paste(sampleidList, collapse=", "))
	#logDebug("setGenesToSameAsSamples colnames(theMatrixFilteredGeneData)=", paste(colnames(theMatrixFilteredGeneData), collapse=", "))
	tfList<-lapply(colnames(theMatrixFilteredGeneData), function(xsample)
			{
				return(is.element(xsample, sampleidList)[[1]])
			})
	tfList <- as.vector(unlist(tfList))
	return(theMatrixFilteredGeneData[,tfList, drop=FALSE])
}

filterSamplesToBatchesDataframe<-function(theDataframeSamplesToBatches, theMatrixFilteredGeneData)
{
	### colnames(theMatrixFilteredGeneData) - sample ids
	### rownames(theMatrixFilteredGeneData) - gene ids
	### names(theDataframeSamplesToBatches) - [1] "sample"     "PlateId"    "ShipDate"   "Batch"      "TSS"        "BCR"        "CGCCandGSC"
	#logDebug("filterSamplesToBatchesDataframe ($Sample from data frame) sampleidList=", paste(theDataframeSamplesToBatches$Sample, collapse=", "))
	#logDebug("filterSamplesToBatchesDataframe colnames(theMatrixFilteredGeneData)=", paste(colnames(theMatrixFilteredGeneData), collapse=", "))
	## DONE INEFFICIENTLY ON PURPOSE, AS THERE IS NO GUARENTEE TO ORDER
	tfList<-lapply(theDataframeSamplesToBatches$Sample, function(xsample)
	{
		return(is.element(xsample, colnames(theMatrixFilteredGeneData))[[1]])
	})
	tfList <- as.vector(unlist(tfList))
	return(theDataframeSamplesToBatches[tfList, , drop=FALSE])
}

################################################################################

closeDevices<-function(theDevices)
{
	for(myDev in theDevices)
	{
		dev.off(which=myDev)
	}
}

writeStringToImage<-function(theFile, theString)
{
	logDebug("writeStringToImage ", theFile, paste(theString, collapse=" "))
	CairoPNG(filename=theFile, pointsize=24, width = 1000, height = 1000)
	on.exit(dev.off())
	plot.new()
	mtext(theString)
}

writeStringToDevice<-function(theDevice, theString)
{
	logDebug("writeStringToDevice ", paste(theDevice, collapse=", "), paste(theString, collapse=" "))
	if (is.list(theDevice))
	{
		for(myDevice in theDevice)
		{
			writeStringToDevice(myDevice, theString)
		}
	}
	else
	{
		dev.set(which = theDevice)
		plot.new()
		mtext(theString, cex=0.5)
	}
}

################################################################################

sortBatchesBasedOnSize<-function(theUniqueBatchIdNames, theListOfAllBatchIds)
{
	countsPerBatch <- sapply(theUniqueBatchIdNames, function (batchId)
			{
				return(sum(as.vector(theListOfAllBatchIds)==batchId))
			}, USE.NAMES=FALSE)
	return(as.character(theUniqueBatchIdNames[order(countsPerBatch, decreasing=TRUE)]))
}

getListOfBatchNamesWithCounts<-function(theListOfAllBatchIds, theUniqueBatchIdNames)
{
	### if this changes the order of the list, the legends produced by PCA will be incorrect
	batchIdentifiersWithCount <- sapply(theUniqueBatchIdNames, function (batchId)
			{
				return(breakIntoTitle(paste(batchId," (",sum(as.vector(theListOfAllBatchIds)==batchId), ")", sep=""), theOldChar=" ", theNewChar=" ", theWidth=30))
			}, USE.NAMES=FALSE)
	return(batchIdentifiersWithCount)
}

getListOfBatchNamesWithCountsAndTotal<-function(theListOfAllBatchIds, theUniqueBatchIdNames)
{
	batchIdentifiersWithCount <- sapply(theUniqueBatchIdNames, function (batchId)
			{
				batchTotal <- sum(as.vector(theListOfAllBatchIds)==batchId)
				remTotal <- (length(theListOfAllBatchIds)-batchTotal)
				return(paste(batchId," (", batchTotal, " vs ", remTotal, ")", sep=""))
			}, USE.NAMES=FALSE)
	return(batchIdentifiersWithCount)
}

widthOfLines<-function(theListOfStrings)
{
	counter <- 1
	for(myString in theListOfStrings)
	{
		tempResult <- nchar(myString)
		if (tempResult>counter)
		{
			counter <- tempResult
		}
	}
	return(counter)
}

################################################################################

sortXBasedOnOrderOfY<-function(x, y)
{
	return(x[factor(order(y))])
}

################################################################################

combineStringPairs<-function(theLabelList, theValueList, theSep, thePadFlag, theRoundFlag, theLessThanReplacement)
{
	if (TRUE==thePadFlag)
	{
		maxLen <- max(nchar(theLabelList))
		maxLen <- maxLen + 1
		theLabelList <- unlist(as.vector(lapply(theLabelList, function(x)
								{
									while(nchar(x)<maxLen)
									{
										x <- paste(x, " ")
									}
									x
								})))
	}
	if (TRUE==theRoundFlag)
	{
		theValueList <- unlist(as.vector(lapply(theValueList, function(x)
								{
									if (is.numeric(x))
									{
										cast <- as.numeric(x)
										if (!is.na(cast))
										{
											return(round(cast, digits=3))
										}
										else
										{
											return(x)
										}
									}
									else
									{
										return(NA);
									}
								})))
	}
	###if (TRUE==theLessThanReplacement)
	###{
	###	theValueList <- unlist(as.vector(lapply(theValueList, function(x)
	###			{
	###				cast <- as.numeric(x)
	###				if (!is.na(cast))
	###				{
	###					if (cast<0.005)
	###					{
	###						return("<0.005")
	###					}
	###					else
	###					{
	###						return(x)
	###					}
	###				}
	###				else
	###				{
	###					return(x)
	###				}
	###			})))
	###}
	paste(theLabelList, theValueList, sep=theSep)
}

################################################################################

writeGeneralAnnotations <- function(theOutputDir, theFilename)
{
	checkCreateDir(file.path(theOutputDir))
	outputFile <- file.path(theOutputDir, theFilename)
	#####
	version <- getMBatchVersion()
	if (!is.null(installed.packages()[ ,"Version"]["MBatch"]))
	{
		version <- as.character(unlist(installed.packages()[ ,"Version"]["MBatch"]))
	}
	#####
	foo <- matrix(c("MBatch Version", version), ncol=2)
	colnames(foo) <- c("key", "value")
	writeAsGenericMatrixNoRows(outputFile, foo)
}

#trimTo5000AndWrite <- function(theMatrix, theOutputDir, theFilename)
#{
#	trimmedData <- mbatchTrimData(theMatrix, ncol(theMatrix)*5000)
#	checkCreateDir(file.path(theOutputDir))
#	outputFile <- file.path(theOutputDir, theFilename)
#	writeDataToFile(trimmedData, outputFile)
#}

writeDataToFile<-function(theData, theFile)
{
	logDebug("Write to file ", theFile)
	writeAsGenericMatrix(theFile, theData)
	logDebug("Finished write to file ", theFile)
}

################################################################################

###batchList: A vector with the list of batches, should be in the same
### order as the samples are in your matrix the cutoff number
batchTreshold <- function(batchList, threshold)
{
	rmBatches <- names(which(table(batchList) < threshold))
	bTreshold <- sapply(batchList, function(thres)
			{
				rmBatchIndex <- which(rmBatches == thres)
				length(rmBatchIndex)==0
			}, USE.NAMES=FALSE)
	return(bTreshold)
}

filterTheGeneDataOnBatchSize<-function(theMatrixGeneData, theMinBatchSize, theBatchIds)
{
	### filter by batch size
	###logDebug("filterTheGeneDataOnBatchSize start")
	if (theMinBatchSize>=1)
	{
		logDebug("filterTheGeneDataOnBatchSize theMinBatchSize=", theMinBatchSize)
		threshHoldList<-batchTreshold(theBatchIds, theMinBatchSize)
		###logDebug("filterTheGeneDataOnBatchSize length(theBatchIds)= ", length(theBatchIds))
		###logDebug("filterTheGeneDataOnBatchSize nrow(theMatrixGeneData)= ", nrow(theMatrixGeneData))
		###logDebug("filterTheGeneDataOnBatchSize ncol(theMatrixGeneData)= ", ncol(theMatrixGeneData))
		###logDebug("filterTheGeneDataOnBatchSize length(threshHoldList)= ", length(threshHoldList))
		###logDebug("filterTheGeneDataOnBatchSize threshHoldList = ", paste(threshHoldList, collapse=","))
		theMatrixGeneData <- theMatrixGeneData[,which(threshHoldList)]
	}
	return(theMatrixGeneData)
}

filterBatchIdsOnBatchSize<-function(theBatchIds, theMinBatchSize)
{
	###logDebug("filterBatchIdsOnBatchSize start")
	if (theMinBatchSize>=1)
	{
		logDebug("filterBatchIdsOnBatchSize theMinBatchSize=", theMinBatchSize)
		threshHoldList<-batchTreshold(theBatchIds, theMinBatchSize)
		###logDebug("filterBatchIdsOnBatchSize length(theBatchIds)= ", length(theBatchIds))
		###logDebug("filterBatchIdsOnBatchSize length(threshHoldList)= ", length(threshHoldList))
		###logDebug("filterBatchIdsOnBatchSize threshHoldList = ", paste(threshHoldList, collapse=","))
		return(theBatchIds[which(threshHoldList)])
	}
	else
	{
		return(theBatchIds)
	}
}

################################################################################

convertDataFrameToSi<-function(theDataFrame)
{
	logDebug("convertDataFrameToSi start")
	### data frame col 1=sample id, col 2+ are batch types
	### values are batch ids
	logDebug("convertDataFrameToSi asmatrixWithIssues")
	foo<-asmatrixWithIssues(theDataFrame[,2:ncol(theDataFrame)], nrow=nrow(theDataFrame), ncol=ncol(theDataFrame)-1)
	### rownames of SI are sample ids
	logDebug("convertDataFrameToSi rownames")
	rownames(foo)<-theDataFrame[,1]
	### colnames of SI are batch types
	logDebug("convertDataFrameToSi colnames")
	colnames(foo)<-colnames(theDataFrame)[2:ncol(theDataFrame)]
	logDebug("convertDataFrameToSi done")
	foo
}

################################################################################

buildList<-function(theValue, theLength)
{
	### looks weird so that it works with NA and NAN
	return(as.vector(unlist(lapply(c(1:theLength), function(x) {return(theValue)}))))
}

################################################################################

### Writes node name and session info and passed liste of objects to file
### settings_[nodename].txt in designated output directory unless directory
### and filename designated in method call
###
writeSettings<-function( theObjectList=c(),
		settingsOutputDir=getwd(),
		settingsFileName=sprintf("settings_%s.txt", Sys.info()['nodename'])
)
{
	nodename=Sys.info()['nodename']

	sink(type = c("output", "message"), file=file.path(settingsOutputDir,settingsFileName))
	cat("nodename:", nodename, "\n\n", sep="")
	cat("working directory:", getwd(), "\n\n", sep="")
	cat("MBatch Version:", getMBatchVersion(), "\n\n", sep="")
	cat("--------------------------------------\nsessionInfo:\n")

	print(sessionInfo())

	if (length(theObjectList)>0)
	{
		cat("\n======================================\n")
		cat("Startup Settings:\n")
		cat("--------------------------------------\n")

		lapply(theObjectList, displayObject)
	}

	sink()
}

### Used by method writeSettings to display an object followed by a dividing line
displayObject<-function(object)
{
	show(object)
	cat("--------------------------------------\n")
}

################################################################################

checkPackageSettings<-function()
{
	stopifnotWithLogging("Package requires locale LC_COLLATE to be set to C for consistency in sorting.", ("C"==Sys.getlocale("LC_COLLATE")))
}

###	collateOrigValue<-Sys.getlocale("LC_COLLATE")
###	on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue))
###	Sys.setlocale("LC_COLLATE","C")

################################################################################

pvalueZeroCheck<-function(theValue, thePermutationCount)
{
	result <- as.character(theValue)
	if(result<(1/thePermutationCount))
	{
		result <- paste("<", format(1/thePermutationCount, scientific=FALSE), sep="")
	}
	return(result)
}

epsilonZeroCheck<-function(theValue)
{
	### if theValue is less than or equal 1x10^-7, return 0, otherwise return theValue
	if (is.na(theValue))
	{
		return(theValue)
	}
	else if (theValue<=0.0000001)
	{
		return(0)
	}
	else
	{
		return(theValue)
	}
}

minEpsilon<-function(theValueVector, theValue)
{
  tmpMin <- min(theValueVector, theValue)
  epsVal <- epsilonZeroCheck(abs(abs(tmpMin)-abs(theValue)))
  if (!is.na(epsVal))
  {
    if (0==epsVal)
    {
      tmpMin <- NA
    }
  }
  tmpMin
}

maxEpsilon<-function(theValueVector, theValue)
{
  tmpMax <- max(theValueVector, theValue)
  epsVal <- epsilonZeroCheck(abs(abs(tmpMax)-abs(theValue)))
  if (!is.na(epsVal))
  {
    if (0==epsVal)
    {
      tmpMax <- NA
    }
  }
  tmpMax
}

################################################################################
####################################################################
###
####################################################################

writeBatchDataFiles<-function(theOutputDir, theDataframeSamplesToBatches, theFilename)
{
	logDebug("write writeBatchDataFiles")
	checkCreateDir(file.path(theOutputDir))
	outputFile <- file.path(theOutputDir, theFilename)
	writeAsGenericDataframe(outputFile, theDataframeSamplesToBatches)
	logDebug("finished writeBatchDataFiles")
}

####################################################################
###
####################################################################

createGraphic_PCAPlusDscPoints<-function(theOutputFile, theCompX, theCompY, theDataDir, theBatchDataDir)
{
	rasterizeJs <- system.file("js", "rasterizePCA.js", package="MBatch")
	pcaHtml <- system.file("js", "pca-plus-dsc-points.html", package="MBatch")
	args = paste("?PCX=PC", theCompX, "&PCY=PC", theCompY,"&BATCHDATA_DIR=", URLencode(theBatchDataDir),"&DATA_DIR=", URLencode(theDataDir),sep="")
	createJsGraphic(rasterizeJs, pcaHtml, theOutputFile, args)
}

createGraphic_PCAPlusDscPointsRays<-function(theOutputFile, theCompX, theCompY, theDataDir, theBatchDataDir)
{
	rasterizeJs <- system.file("js", "rasterizePCA.js", package="MBatch")
	pcaHtml <- system.file("js", "pca-plus-dsc-points-rays.html", package="MBatch")
	args = paste("?PCX=PC", theCompX, "&PCY=PC", theCompY,"&BATCHDATA_DIR=", URLencode(theBatchDataDir),"&DATA_DIR=", URLencode(theDataDir),sep="")
	createJsGraphic(rasterizeJs, pcaHtml, theOutputFile, args)
}

createJsGraphic<-function(theRasterizeFile, theHtmlFile, theOutputFile, theUrlArgs)
{
	system(paste("phantomjs \"", theRasterizeFile, "\" \"file:///", theHtmlFile, theUrlArgs, "\" \"",theOutputFile, sep=""))
}

####################################################################
###
####################################################################

filterBasedOnGeneLimit<-function(theMatrixGeneData, theSeed, theGeneLimit)
{
	if (theGeneLimit>=1)
	{
		logDebug("filterBasedOnGeneLimit theGeneLimit=", theGeneLimit)
		if(length(rownames(theMatrixGeneData))>theGeneLimit)
		{
	    set.seed(theSeed)
			logDebug("filterBasedOnGeneLimit trimming")
			theMatrixGeneData <- mbatchTrimData(theMatrixGeneData, length(colnames(theMatrixGeneData))*theGeneLimit)
		}
	}
	logDebug("filterBasedOnGeneLimit theMatrixGeneData nrow=", nrow(theMatrixGeneData))
	logDebug("filterBasedOnGeneLimit theMatrixGeneData ncol=", ncol(theMatrixGeneData))
	return(theMatrixGeneData)
}

####################################################################
###
####################################################################

mbatchStandardLegend <- function(theTitle, theVersion, theLegendNames, theLegendColors, theLegendSymbols, theFilenamePath)
{
	myClass1 <- system.file("LegendJava", "jcommon-1.0.17.jar", package="MBatch")
	myClass2 <- system.file("LegendJava", "jfreechart-1.0.14.jar", package="MBatch")
	myClass3 <- system.file("LegendJava", "LegendJava.jar", package="MBatch")
	myJavaJars <- file.path(myClass3, myClass1, myClass2, fsep=.Platform$path.sep)
	logDebug("mbatchStandardLegend - Calling .jinit ", myJavaJars)
	.jinit(classpath=myJavaJars, force.init = TRUE, parameters="-Xms4800m")
	logDebug("mbatchStandardLegend - .jinit complete")
	logDebug("mbatchStandardLegend - theTitle ", theTitle)
	logDebug("mbatchStandardLegend - theVersion ", theVersion)
	logDebug("mbatchStandardLegend - theFilenamePath ", theFilenamePath)
	logDebug("mbatchStandardLegend - theLegendNames ", paste(theLegendNames, collapse=", "))
	logDebug("mbatchStandardLegend - theLegendNames ", length(theLegendNames))
	logDebug("mbatchStandardLegend - theLegendColors ", length(theLegendColors))
	logDebug("mbatchStandardLegend - theLegendSymbols ", length(theLegendSymbols))
	myColors <- NULL
	if (!is.null(theLegendColors))
	{
		foo <- col2rgb(theLegendColors)
		for(myIndex in 1:(length(foo)/3))
		{
			myColors <- c(myColors, sprintf("#%02x%02x%02x", foo[,myIndex]["red"], foo[,myIndex]["green"], foo[,myIndex]["blue"]))
		}
	}
	logDebug("mbatchStandardLegend - myColors ", paste(myColors, collapse=","))
	logDebug("mbatchStandardLegend before java")
	result <- .jcall("org/mda/legendjava/LegendJava", returnSig = "Z",
													method='writeLegend',
													.jnew("java/lang/String",theTitle),
													.jnew("java/lang/String",theVersion),
													.jarray(as.vector(as.character(theLegendNames))),
									 				.jarray(as.vector(as.character(myColors))),
									 				.jarray(as.vector(as.integer(theLegendSymbols)), contents.class="[I"),
													.jnew("java/lang/String",theFilenamePath))
	if(FALSE==result)
	{
		stop("Failed to write legend")
	}
	logDebug("mbatchStandardLegend after java")
}

mbatchStandardCombineLegends<-function(theTitle, theFilenamePath, theListOfFiles)
{
	myClass1 <- system.file("LegendJava", "jcommon-1.0.17.jar", package="MBatch")
	myClass2 <- system.file("LegendJava", "jfreechart-1.0.14.jar", package="MBatch")
	myClass3 <- system.file("LegendJava", "LegendJava.jar", package="MBatch")
	myJavaJars <- file.path(myClass3, myClass1, myClass2, fsep=.Platform$path.sep)
	logDebug("mbatchStandardCombineLegends - Calling .jinit ", myJavaJars)
	.jinit(classpath=myJavaJars, force.init = TRUE, parameters="-Xms4800m")
	logDebug("mbatchStandardCombineLegends - .jinit complete")
	logDebug("mbatchStandardCombineLegends - theTitle ", theTitle)
	logDebug("mbatchStandardCombineLegends - theFilenamePath ", theFilenamePath)
	logDebug("mbatchStandardCombineLegends - theListOfFiles ", paste(theListOfFiles, collapse=", "))
	logDebug("mbatchStandardLegend before java")
	result <- .jcall("org/mda/legendjava/LegendJava", returnSig = "Z",
									 method='combineLegends',
									 .jnew("java/lang/String",theTitle),
									 .jarray(as.vector(as.character(theListOfFiles))),
									 .jnew("java/lang/String",theFilenamePath))
	if(FALSE==result)
	{
		stop("Failed to write combined legend")
	}
	logDebug("mbatchStandardLegend after java")
}

####################################################################
###
####################################################################

getTestInputDir <- function()
{
  value <- Sys.getenv("MBATCH_TEST_INPUT")
  if (!isTRUE(file.exists(value)))
  {
    value <- "/BatchEffectsPackage_data/testing_static/MATRIX_DATA"
  }
  value
}

getTestOutputDir <- function()
{
  value <- Sys.getenv("MBATCH_TEST_OUTPUT")
  if (!isTRUE(file.exists(value)))
  {
    value <- "/BatchEffectsPackage_data/testing_dynamic/MBatch"
  }
  value
}

getTestCompareDir <- function()
{
  value <- Sys.getenv("MBATCH_TEST_COMPARE")
  if (!isTRUE(file.exists(value)))
  {
    value <- "/BatchEffectsPackage_data/testing_static/COMPARE"
  }
  value
}

####################################################################
###
####################################################################

compareTwoMatrices <- function(theCorrected, theCompare)
{
  correctedRows <- dim(theCorrected)[1]
  correctedCols <- dim(theCorrected)[2]
  compareRows <- dim(theCompare)[1]
  compareCols <- dim(theCompare)[2]

  stopifnot(correctedRows==compareRows)
  stopifnot(correctedCols==compareCols)

  for(myCol in 1:correctedCols)
  {
    if (!(colnames(theCorrected)[myCol]==colnames(theCompare)[myCol]))
    {
      message("myCol=", myCol)
      message("colnames(theCorrected)[myCol]=", colnames(theCorrected)[myCol])
      message("colnames(theCompare)[myCol]=", colnames(theCompare)[myCol])
      return(FALSE)
    }
    for(myRow in 1:correctedRows)
    {
      if (!(rownames(theCorrected)[myRow]==rownames(theCompare)[myRow]))
      {
        message("myRow=", myRow)
        message("rownames(theCorrected)[myRow]=", rownames(theCorrected)[myRow])
        message("rownames(theCompare)[myRow]=", rownames(theCompare)[myRow])
        return(FALSE)
      }
      #message("checking myRow=", myRow, " and myCol=", myCol, " theCorrected[myRow, myCol]=", theCorrected[myRow, myCol], " theCompare[myRow, myCol]=",theCompare[myRow, myCol])
      if (!is.finite(theCorrected[myRow, myCol])&&!is.finite(theCompare[myRow, myCol]))
      {
        # ignore
      }
      else if (!(all.equal(theCorrected[myRow, myCol], theCompare[myRow, myCol])))
      {
        message("myRow=", myRow)
        message("myCol=", myCol)
        message("theCorrected[myRow, myCol]=", theCorrected[myRow, myCol])
        message("theCompare[myRow, myCol]=", theCompare[myRow, myCol])
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

####################################################################
###
####################################################################
