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

VALUE_TYPE_LIST <- c('metric','p-value')
CZ_NAME_LIST <- c('CZR', 'CZN', 'CZNL', 'CZNS','CZS')

doBatchcorrPreTest<-function(theBatchcorrData)
{
	if (!is.null(theBatchcorrData))
	{
		if((TRUE==theBatchcorrData@mDoOtoMFlag)||(TRUE==theBatchcorrData@mDoMtoMFlag))
		{
			### only test if batchcorr run requested
			stopifnotWithLogging("The min batch requested for batchcorr must be at least 4", (theBatchcorrData@mMinBatchSize>=4))
		}
	}
}

createBatchEffectsOutput_batchcorr<-function(theMatrixGeneData, theDataframeBatchData, theTitle,
		theMinNumberOfGenes, theNumberOfPermutatedGenes, theNumberOfPermutations,
		theMinBatchSize=20, theAdjustedFlag=TRUE, theNumberOfThreads=1,
		theContinueFunction=NULL, theContinueTestCount=0, theSeed=NULL,
		theOutputDir=getwd(), theBatchCorrBase="BatchCorr")
{
	logDebug("createBatchEffectsOutput_batchcorr 1")
	collateOrigValue<-Sys.getlocale("LC_COLLATE")
	on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
	Sys.setlocale("LC_COLLATE","C")
	logDebug("Changing LC_COLLATE to C for duration of run")
	checkPackageSettings()
	if(isNamespaceLoaded("batchcorr"))
	{
		resultList <- batchcorrForBatchType_calculate(theMatrixGeneData, theDataframeBatchData,
			theMinNumberOfGenes, theNumberOfPermutatedGenes, theNumberOfPermutations,
			theMinBatchSize, theAdjustedFlag, theNumberOfThreads, theSeed,
			theContinueFunction, theContinueTestCount)
		### RData output write resultList and theDataframeBatchData for batchcorr
		### TODO: D3 batchcorr data
		###checkCreateDir(theOutputDir, theBatchCorrBase)
		###outputFile <- cleanFilePath(cleanFilePath(theOutputDir, theBatchCorrBase), "BatchCorr_Many2Many_Calc.RData"))
		###logDebug("RData-output write resultList and theDataframeBatchData for batchcorr ", outputFile)
		###save(resultList, theDataframeBatchData, file=outputFile)
		logDebug("createBatchEffectsOutput_batchcorr 2")
		batchcorrForBatchType_writeImage(resultList, theTitle, theOutputDir, theBatchCorrBase)
	}
	logDebug("createBatchEffectsOutput_batchcorr 3")
}

createBatchEffectsOutput_batchcorr_one2many<-function(theMatrixGeneData, theDataframeBatchData, theTitle,
		theMinNumberOfGenes, theNumberOfPermutatedGenes, theNumberOfPermutations,
		theMinBatchSize=20, theAdjustedFlag=TRUE, theNumberOfThreads=1,
		theContinueFunction=NULL, theContinueTestCount=0, theSeed=NULL,
		theOutputDir=getwd(), theBatchCorrBase="BatchCorr")
{
	logDebug("createBatchEffectsOutput_batchcorr_one2many 1")
	collateOrigValue<-Sys.getlocale("LC_COLLATE")
	on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
	Sys.setlocale("LC_COLLATE","C")
	logDebug("Changing LC_COLLATE to C for duration of run")
	checkPackageSettings()
	if(require(batchcorr, warn.conflicts=FALSE))
	{
		resultList <- one2ManyBatchcorrForBatchType_calculate(theMatrixGeneData, theDataframeBatchData,
						theMinNumberOfGenes, theNumberOfPermutatedGenes, theNumberOfPermutations,
						theMinBatchSize, theAdjustedFlag, theNumberOfThreads, theSeed,
						theContinueFunction, theContinueTestCount)
		### RData output write resultList for batchcorr one2many
		### TODO: D3 batchcorr data
		### checkCreateDir(theOutputDir, theBatchCorrBase)
		### outputFile <- cleanFilePath(cleanFilePath(theOutputDir, theBatchCorrBase), "BatchCorr_One2Many_Calc.RData")
		### logDebug("RData output write resultList and theDataframeBatchData for batchcorr one2many ", outputFile)
		### save(resultList, theDataframeBatchData, file=outputFile)
		logDebug("createBatchEffectsOutput_batchcorr_one2many 2")
		one2ManyBatchcorrForBatchType_images(resultList, theTitle, theOutputDir, theBatchCorrBase)
	}
	logDebug("createBatchEffectsOutput_batchcorr_one2many 3")
}

####################################################################
####################################################################

one2ManyBatchcorrForBatchType_images<-function(theResultList, theTitle, theOutputDir, theBatchCorrBase)
{
	collateOrigValue<-Sys.getlocale("LC_COLLATE")
	on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
	Sys.setlocale("LC_COLLATE","C")
	logDebug("Changing LC_COLLATE to C for duration of run")
	checkPackageSettings()
	if(require(batchcorr, warn.conflicts=FALSE))
	{
		### loop for each batch type
		for(batchType in names(theResultList))
		{
			batchTypeList <- theResultList[[batchType]]
			### get metric results
			### loop for each unique batch id in this batch type
			metricList<-lapply(names(batchTypeList), function(uniqueBatchId)
					{
						results<-batchTypeList[[uniqueBatchId]]
						metricListResult<-buildList(NaN, 5)
						if (TRUE==is.array(results))
						{
							metricResults<-results["metric",]
							metricListResult<-as.vector(unlist(metricResults))
						}
						return(metricListResult)
					}
			)
			metricList<-as.vector(unlist(metricList))
			###
			pvalueList<-lapply(names(batchTypeList), function(uniqueBatchId)
					{
						results<-batchTypeList[[uniqueBatchId]]
						pvalueListResult<-buildList(NaN, 5)
						if (TRUE==is.array(results))
						{
							pvalueResults<-results["p-value",]
							pvalueListResult<-as.vector(unlist(pvalueResults))
						}
						return(pvalueListResult)
					}
			)
			pvalueList<-as.vector(unlist(pvalueList))
			###
			batchIds<-names(batchTypeList)
			###VALUE_TYPE_LIST <- c('metric','p-value')
			###CZ_NAME_LIST <- c('CZR', 'CZN', 'CZNL', 'CZNS','CZS')
			### dimnames = c(rowNames, colNames)
			metricMatrix <- matrixWithIssues(metricList, nrow=length(batchIds), ncol=length(CZ_NAME_LIST), dimnames=list(batchIds, CZ_NAME_LIST))
			pvalueMatrix <- matrixWithIssues(pvalueList, nrow=length(batchIds), ncol=length(CZ_NAME_LIST), dimnames=list(batchIds, CZ_NAME_LIST))
			outputDir <- checkCreateDir(checkCreateDir(checkCreateDir(theOutputDir, theBatchCorrBase), batchType), "OneToMany")
			###logDebug("batchcorrForBatchType_writeImage outputDir ", outputDir)
			metricFile <- createDirPlusFilename(outputDir, theBatchCorrBase, "_Metric", "_", "ALL", "_Diagram.png")
			pvalueFile <- createDirPlusFilename(outputDir, theBatchCorrBase, "_Pvalue", "_", "ALL", "_Diagram.png")
			###logDebug("batchcorrForBatchType_writeImage metricFile ", metricFile)
			###logDebug("batchcorrForBatchType_writeImage pvalueFile ", pvalueFile)
			writePvalueImage(metricFile, pvalueFile, "ALL", "metric", metricMatrix, paste(theTitle, "Metric", batchType, sep=" "))
			writePvalueImage(metricFile, pvalueFile, "ALL", "p-value", pvalueMatrix, paste(theTitle, "P-Value", batchType, sep=" "))
		}
		###batchcorrForBatchType_writeImage(resultList, theTitle, theBatchCorrMetricFile, theBatchCorrPValueFile)
	}
}

####################################################################
####################################################################

batchcorrForBatchType_calculate<-function(theMatrixGeneData, theDataframeBatchData,
		theMinNumberOfGenes, theNumberOfPermutatedGenes, theNumberOfPermutations,
		theMinBatchSize=20, theAdjustedFlag=TRUE, theNumberOfThreads=1, theSeed=NULL,
		theContinueFunction=NULL, theContinueTestCount=0)
{
	logDebug("batchcorrForBatchType_calculate")
	collateOrigValue<-Sys.getlocale("LC_COLLATE")
	on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
	Sys.setlocale("LC_COLLATE","C")
	logDebug("Changing LC_COLLATE to C for duration of run")
	checkPackageSettings()
	stopifnotWithLogging("The number of columns of gene data must equal the number of rows of batch data.", ncol(theMatrixGeneData)==nrow(theDataframeBatchData))
	if(require(batchcorr, warn.conflicts=FALSE))
	{
		### for each batch type
		resultList<-lapply(c(2:length(theDataframeBatchData)), function(batchTypeIndex)
				{
					return(batchcorrForBatchType_calcForSingleIndex(batchTypeIndex, theMatrixGeneData, theDataframeBatchData,
									theMinNumberOfGenes, theNumberOfPermutatedGenes, theNumberOfPermutations,
									theMinBatchSize, theAdjustedFlag, theNumberOfThreads, theSeed,
									theContinueFunction, theContinueTestCount))
				}
		)
		names(resultList) <- names(theDataframeBatchData)[2:length(theDataframeBatchData)]
		return(resultList)
	}
	else
	{
		return(NULL)
	}
}

one2ManyBatchcorrForBatchType_calculate<-function(theMatrixGeneData, theDataframeBatchData,
		theMinNumberOfGenes, theNumberOfPermutatedGenes, theNumberOfPermutations,
		theMinBatchSize=20, theAdjustedFlag=TRUE, theNumberOfThreads=1, theSeed=NULL,
		theContinueFunction=NULL, theContinueTestCount=0)
{
	logDebug("one2ManyBatchcorrForBatchType_calculate_calculate")
	logDebug("inside one2ManyBatchcorrForBatchType_calculate_calculate")
	stopifnotWithLogging("The number of columns of gene data must equal the number of rows of batch data.", ncol(theMatrixGeneData)==nrow(theDataframeBatchData))
	if(require(batchcorr, warn.conflicts=FALSE))
	{
		### for each batch type
		### loop for each batch type
		resultList <- lapply(c(2:length(theDataframeBatchData)), function(batchTypeIndex)
				{
					return(one2manyBatchcorr_calcForSingleIndex(batchTypeIndex, theMatrixGeneData, theDataframeBatchData,
									theMinNumberOfGenes, theNumberOfPermutatedGenes, theNumberOfPermutations,
									theMinBatchSize, theAdjustedFlag, theNumberOfThreads, theSeed,
									theContinueFunction, theContinueTestCount))
				}
		)
		names(resultList)<- names(theDataframeBatchData)[2:length(theDataframeBatchData)]
		return(resultList)
	}
	else
	{
		return(NULL)
	}
}

####################################################################
####################################################################

one2manyBatchcorr_calcForSingleIndex<-function(batchTypeIndex, theMatrixGeneData, theDataframeBatchData,
		theMinNumberOfGenes, theNumberOfPermutatedGenes, theNumberOfPermutations,
		theMinBatchSize, theAdjustedFlag, theNumberOfThreads, theSeed=NULL,
		theContinueFunction=NULL, theContinueTestCount=0)
{
	if(require(batchcorr, warn.conflicts=FALSE))
	{
		batchTypeName <- names(theDataframeBatchData)[batchTypeIndex]
		logDebug("one2manyBatchcorr_calcForSingleIndex batchTypeName ", batchTypeName)
		logDebug("one2manyBatchcorr_calcForSingleIndex theMinBatchSize ", theMinBatchSize)
		batchIdsForSamples <- as.vector(unlist(theDataframeBatchData[batchTypeIndex]))
		uniqueListOfBatchIds <- sort(unique(batchIdsForSamples))
		batchIdentifiersWithCount <- getListOfBatchNamesWithCountsAndTotal(batchIdsForSamples, uniqueListOfBatchIds)
		### loop for each unique batch id in this batch type
		resultList<-lapply(uniqueListOfBatchIds, function(batchId)
				{
					### duplicate theDataframeBatchData$Sample
					oneToManyDataframe <- data.frame(sample=theDataframeBatchData$Sample, stringsAsFactors=FALSE, check.names=FALSE)
					### duplicate theDataframeBatchData[batchTypeIndex] into theDataframeBatchData$<batchTypeName>
					###		replacing all but batchId with "Other Batches"
					oneToManyDataframe$newbatch<-sapply(batchIdsForSamples, function(currentBatchId)
							{
								if(currentBatchId==batchId)
								{
									return(currentBatchId)
								}
								else
								{
									return("Other Batches")
								}
							}, USE.NAMES=FALSE )
					names(oneToManyDataframe)[2] <- batchTypeName
					###
					stopifnotWithLogging("The number of columns of gene data must equal the number of rows of batch data.", ncol(theMatrixGeneData)==nrow(oneToManyDataframe))
					### list of results with names(tempResultList) being list of batch Types
					tempResultList<-batchcorrForBatchType_calculate(theMatrixGeneData, oneToManyDataframe,
							theMinNumberOfGenes, theNumberOfPermutatedGenes, theNumberOfPermutations,
							theMinBatchSize, theAdjustedFlag, theNumberOfThreads, theSeed,
							theContinueFunction, theContinueTestCount)
					tempResult<-tempResultList[[1]]
					czArray<-tempResult
					if (TRUE==is.array(tempResult))
					{
						czArray<-tempResult[1,2,,]
					}
					return(czArray)
			}
		)
		names(resultList)<-batchIdentifiersWithCount
		return(resultList)
	}
	else
	{
		return(NULL)
	}
}

####################################################################
####################################################################

getColors<-function(theOutputType)
{
	myBlue <-  rainbow(1, start=4/6, s=0.7)
	myColorRampPalette <- colorRampPalette(c("red", "white", myBlue))
	### this is a hard coded value but not necessarily related to any other value.
	### just cut into 20 colors
	myColors<-beaRainbow(20, s=0.7)
	if ("metric"==theOutputType)
	{
		myColors <- myColorRampPalette(100)
		###myBreaks <- seq(0, 1, length.out=101)
	}
	else
	{
		myColors <- c(c("red", "orange", "yellow", "green"), rainbow(10, start=0.50, end=0.67, s=0.7))
	}
	return(myColors)
}

getBreaks<-function(theOutputType)
{
	### this is a hard coded value but not necessarily related to any other value.
	### just cut into 20 colors
	myBreaks <- NULL
	if ("metric"==theOutputType)
	{
		myBreaks <- NULL
	}
	else
	{
		myBreaks <- c(c(0, 0.05, 0.1, 0.15, 0.2), seq(0.3, 1, length.out=10))
	}
	return(myBreaks)
}

writePvalueImage<-function(theBatchCorrMetricFile, theBatchCorrPValueFile,
		theBatchType, theOutputType, resultsBatch, theTitle)
{
	###logDebug("writePvalueImage theTitle ", theTitle)
	myColors<-getColors(theOutputType)
	myBreaks<-getBreaks(theOutputType)
	outputFile <- theBatchCorrPValueFile
	if ("metric"==theOutputType)
	{
		outputFile <- theBatchCorrMetricFile
	}
	tempTitle <- gsub("_", " ", theTitle, fixed=TRUE)
	###logDebug("writePvalueImage tempTitle ", tempTitle)
	###logDebug("writePvalueImage outputFile ", outputFile)
	writeBatchCoreFile(outputFile,
			theBatchType,
			tempTitle,
			resultsBatch,
			myColors,
			myBreaks)
}

writeBatchCoreFile<-function(theFile, theBatchType, theTitle, theData, theColors, theBreaks)
{
	###logDebug("writeBatchCoreFile file = ", theFile)
	CairoPNG(filename=theFile, pointsize=24, width = 1000, height = 1000)
	on.exit(dev.off())
	imageFromMatrix(data=theData, col=theColors, breaks=theBreaks,
			main=list(breakIntoTitle(theTitle, theOldChar = " ", theNewChar = " ", theWidth = 50)))
	mtext(paste("MBatch", packageDescription("MBatch")["Version"], "batchcorr", packageDescription("batchcorr")["Version"], sep=" "), side=1, adj=1, outer=FALSE, line=4, cex=0.50, bg="white")
}

imageFromMatrix<-function(data, col=beaRainbow(12), breaks=NULL, main=NULL)
{
	### data should be a matrix with column and row labels (names)
	###logDebug("imagewithkey image with key")
	omar<-par("mar")
	###logDebug("imagewithkey 1")
	evenBreaks <- t(seq(0,1, length.out=(length(col))))
	evenPlacement <- seq(0,1, length.out=(1+length(col)))
	axisLabelCex <- 0.7
	legendLabelCex <- 0.7
	chartValueCex <- 1.0
	###logDebug("imagewithkey 2")
	splitColNames <- sapply(colnames(data), function (mylabel)
			{
				return(breakIntoTitle(mylabel, theOldChar=" ", theNewChar=" ", theWidth=18))
			}, USE.NAMES=FALSE)
	###logDebug("imagewithkey 3")
	splitRowNames <- sapply(rownames(data), function (mylabel)
			{
				return(breakIntoTitle(mylabel, theOldChar=" ", theNewChar=" ", theWidth=18))
			}, USE.NAMES=FALSE)
	###logDebug("imagewithkey 4")
	if ((length(splitColNames)>=15)||(length(splitRowNames)>=15))
	{
		chartValueCex <- 0.5
		axisLabelCex <- 0.5
		legendLabelCex <- 0.5
	}
	else if ((length(splitColNames)>=10)||(length(splitRowNames)>=10))
	{
		chartValueCex <- 0.6
		axisLabelCex <- 0.5
		legendLabelCex <- 0.5
	}
	else if ((length(splitColNames)>=5)||(length(splitRowNames)>=5))
	{
		chartValueCex <- 0.75
		axisLabelCex <- 0.6
		legendLabelCex <- 0.6
	}
	###logDebug("imagewithkey 5")
	###labelWidth <- max(strwidth(rownames(data), cex=cex.label, units="inches"))
	###labelWidth <- max(strwidth(rownames(data), cex=3.6, units="inches"))
	labelWidth <- 2
	###c(bottom, left, top, right)
	par(mar=c(omar[1]+labelWidth,omar[2]+labelWidth,omar[3],1))
	### this writes the title, the axis labels (not the row and column labels), and colors the squares
	layout(mat=matrix(1:2, 1,2), c(10,1))
	###logDebug("imagewithkey 6")
	if (is.null(breaks))
	{
		image(data, col=col, axes=FALSE, zlim=c(0,1), main=main)
	}
	else
	{
		image(data, col=col, axes=FALSE, breaks=breaks, main=main)
	}
	### this adds the row and column labels
	###logDebug("imagewithkey 7")
	axis(2, at=seq(0, 1, length.out=ncol(data)), lwd=1, labels=splitColNames, cex.axis=axisLabelCex, las=2, lab=c(5,5,20) )
	axis(1, at=seq(0, 1, length.out=nrow(data)), lwd=1, labels=splitRowNames, cex.axis=axisLabelCex, las=2, lab=c(5,5,20) )
	### this populates the numbers in each cell
	###logDebug("imagewithkey 8")
	textLocationsX <- seq(0, 1, length.out=nrow(data))
	textLocationsY <- seq(0, 1, length.out=ncol(data))
	rectLocationsX <- seq(0, 1, length.out=nrow(data))
	widthX <- rectLocationsX[2]/2
	rectLocationsY <- seq(0, 1, length.out=ncol(data))
	widthY <- rectLocationsY[2]/2
	###logDebug("imagewithkey chartValueCex ", chartValueCex)
	for(i in 1:ncol(data))
	{
		for(j in 1:nrow(data))
		{
			text(x=textLocationsX[j], y=textLocationsY[i], labels=round(data[j,i], 3), cex=chartValueCex)
			rect(rectLocationsX[j]-widthX, rectLocationsY[i]-widthY, rectLocationsX[j]+widthX, rectLocationsY[i]+widthY, border="black")
		}
	}
	###logDebug("imagewithkey 10")
	### draw a box on the outside of the colored squares
	box()
	### this does the tall colored graph on the right
	###logDebug("imagewithkey 11")
	par(mar=c(omar[1]+labelWidth,0,omar[3],2.1))
	###logDebug("imagewithkey 12")
	image(evenBreaks, axes=FALSE, col=col)
	###logDebug("imagewithkey 13")
	if (!is.null(breaks))
	{
		mtext(side=4, at=evenPlacement, text=round(breaks,2), las=2, cex=legendLabelCex)
	}
	else
	{
		mtext(side=4, at=seq(0, 1, length=3), text=seq(0, 1, length=3), las=2, cex=legendLabelCex)
	}
	###logDebug("imagewithkey 14")
}

####################################################################
####################################################################

batchcorrForBatchType_calcForSingleIndex<-function(batchTypeIndex, theMatrixGeneData, theDataframeBatchData,
		theMinNumberOfGenes, theNumberOfPermutatedGenes, theNumberOfPermutations,
		theMinBatchSize, theAdjustedFlag, theNumberOfThreads, theSeed=NULL,
		theContinueFunction=NULL, theContinueTestCount=0)
{
	resultsBatch<-NULL
	###logDebug("batchcorrForBatchType_calcForSingleIndex - batch index ", batchTypeIndex)
	if (!is.null(theContinueFunction))
	{
	###logDebug("batchcorrForBatchType_calcForSingleIndex - continue function is not null")
	}
	else
	{
	###logDebug("batchcorrForBatchType_calcForSingleIndex - continue function found a null")
	}
	###logDebug("batchcorrForBatchType_calcForSingleIndex - theContinueTestCount = ", theContinueTestCount)
	if (0==theContinueTestCount)
	{
		###logDebug("batchcorrForBatchType_calcForSingleIndex - test count is 0, null out continue function")
		theContinueFunction <- NULL
	}
	### compile data and information for display
	sampleIdsAll <-  as.vector(unlist(theDataframeBatchData[1]))
	batchIdsAll <- as.vector(unlist(theDataframeBatchData[batchTypeIndex]))
	batchIdsUnique <- sort(unique(batchIdsAll))
	if(length(batchIdsUnique)<=1)
	{
		###logDebug("batchcorrForBatchType_calcForSingleIndex - More than one batch required (initial test)")
		resultsBatch<-"More than one batch required (initial test)"
	}
	else
	{
		keepSamples<-colSums(!is.na(theMatrixGeneData))!=0
		myData<-theMatrixGeneData[,keepSamples]
		keepSampleIds<-sampleIdsAll[keepSamples]
		keepBatchIds<-batchIdsAll[keepSamples]
		###OLD FILTER Mad<-apply(myData, 1, mad)
		###OLD FILTER uData<-myData[!is.na(Mad) & Mad> quantile(Mad, na.rm=T)[4],]
		if (nrow(myData)<theMinNumberOfGenes)
		{
			###logDebug("batchcorrForBatchType_calcForSingleIndex - Gene count is less than min number of genes ", nrow(myData), " < ", theMinNumberOfGenes)
			resultsBatch<-"Gene count is less than min number of genes"
		}
		else if (nrow(myData)<theNumberOfPermutatedGenes)
		{
			###logDebug("batchcorrForBatchType_calcForSingleIndex - Gene count is less than requested permutated genes ", nrow(myData), " < ", theNumberOfPermutatedGenes)
			resultsBatch<-"Gene count is less than requested permutated genes"
		}
		else
		{
			###
			validBatch<-names(table(keepBatchIds))[table(keepBatchIds)>=theMinBatchSize]
			nBatch<-length(validBatch)
			if(nBatch<=1)
			{
				###logDebug("batchcorrForBatchType_calcForSingleIndex - More than one batch required (second test)")
				resultsBatch<-"More than one batch required (second test)"
			}
			else
			{

				batchnamesWithCount <- sapply(validBatch, function (mylabel)
						{
							return(paste(mylabel, " (", sum(as.vector(batchIdsAll)==mylabel), ")", sep=""))
						}, USE.NAMES=FALSE)
				resultsBatch<-array(1, dim=c(nBatch, nBatch, 2, 5), dimnames=list(batchnamesWithCount, batchnamesWithCount, VALUE_TYPE_LIST, CZ_NAME_LIST))
				if(nBatch==1)
				{
					resultsBatch[,,,]<-1
				}
				else
				{
					###logDebug("batchcorr-time started")
					###OLD FILTER sampledData<-t(uData[sample(1:nrow(uData), min(theNumberOfPermutatedGenes, nrow(uData))),])
					if (0!=theSeed)
					{
						set.seed(theSeed)
					}
					sampledData<-t(myData[sample(1:nrow(myData), min(theNumberOfPermutatedGenes, nrow(myData))),])
					bcResults<-NULL
					###logDebug("batchcorr comparisons will be ", nBatch)
					###logDebug("theNumberOfThreads will be ", theNumberOfThreads)
					for(i in 1:(nBatch-1))
					{
						for(j in (i+1):nBatch)
						{
							###logDebug("doing batchcorr comparison ", i, j)
							stopifnotWithLogging(paste("The number of batches to process in batch ", i, " is less than the defined minimum batch size of ",theMinBatchSize),
										length(which(keepBatchIds==validBatch[i]))>=theMinBatchSize)
							stopifnotWithLogging(paste("The number of batches to process in batch ", j, " is less than the defined minimum batch size of ",theMinBatchSize),
										length(which(keepBatchIds==validBatch[j]))>=theMinBatchSize)
							###logDebug("batchcorrForBatchType_calcForSingleIndex -- get input matrix")
							input.mat<-asmatrixWithIssues(sampledData[c(which(keepBatchIds==validBatch[i]), which(keepBatchIds==validBatch[j])),])
							###logDebug("input.mat ", paste(dim(input.mat), collapse=", "))
							nbatch1<-sum(keepBatchIds==validBatch[i])
							###logDebug("batchcorrForBatchType_calcForSingleIndex -- nbatch1=", nbatch1)
							###we need the matrix of sample in row and gene in column
							###logDebug("nbatch1 ", nbatch1)
							###logDebug("theNumberOfPermutations ", theNumberOfPermutations)
							logDebug("batchcorrForBatchType_calcForSingleIndex -- theAdjustedFlag ", theAdjustedFlag)
							callBC <- get("batchcorr")
							bcResults<-callBC(input.mat, nbatch1, (theNumberOfPermutations+1), theAdjustedFlag,
										continue=theContinueFunction, minperm=theContinueTestCount, theThreads=theNumberOfThreads)
							###logDebug("batchcorrForBatchType_calcForSingleIndex -- bcResults=", paste(bcResults ,collapse=", "))
							###logDebug("batchcorrForBatchType_calcForSingleIndex -- bcResults[1,]=", paste(bcResults[1,],collapse=", "))
							resultsBatch[j,i,1,]<-resultsBatch[i,j,1,]<-bcResults[1,]
							###logDebug("batchcorrForBatchType_calcForSingleIndex -- get metric values=", paste(bcResults[1,],collapse=", "))
							nanList<-is.na(bcResults[1,])
							###logDebug("batchcorrForBatchType_calcForSingleIndex -- nanList=", paste(nanList,collapse=", "))
							numValidList<-1
							numLessList<-1
							if(nrow(bcResults)>1)
							{
								numValidList<-colSums(!is.na(bcResults[2:nrow(bcResults),]))
								numLessList<-colSums(bcResults[2:nrow(bcResults),]<=matrix(bcResults[1,], byrow=TRUE, ncol=ncol(bcResults), nrow=(nrow(bcResults)-1)), na.rm=TRUE)
							}
							else
							{
							###logDebug("batchcorrForBatchType_calcForSingleIndex -- found all NaNs")
							}
							###logDebug("batchcorrForBatchType_calcForSingleIndex -- numValidList=", paste(numValidList,collapse=", "))
							###logDebug("batchcorrForBatchType_calcForSingleIndex -- numLessList=", paste(numLessList,collapse=", "))
							resultsBatch[j,i,2,]<-resultsBatch[i,j,2,]<-sapply(c(1:length(nanList)), function(myIndex)
									{
										if (TRUE==nanList[[myIndex]])
										{
											###logDebug("batchcorrForBatchType_calcForSingleIndex -- get TRUE==nanList")
											return(NaN)
										}
										else
										{
											###logDebug("batchcorrForBatchType_calcForSingleIndex -- get FALSE==nanList")
											result<-numLessList[[myIndex]]/numValidList[[myIndex]]
											###logDebug("batchcorrForBatchType_calcForSingleIndex -- get after divide")
											if (is.na(result))
											{
												###logDebug("batchcorrForBatchType_calcForSingleIndex -- get is.na")
												result<- NaN
											}
											###logDebug("batchcorrForBatchType_calcForSingleIndex -- almost done")
											return(result)
										}
									}, USE.NAMES=FALSE)
							###logDebug("batchcorrForBatchType_calcForSingleIndex -- looping batchcorr")
						}
					}
					###logDebug("batchcorr finished")
				}
			}
		}
	}
	return(resultsBatch)
}

batchcorrForBatchType_writeImage<-function(theResultsBatchList, theTitle, theOutputDir, theBatchCorrBase)
{
	###logDebug("batchcorrForBatchType_writeImage")
	collateOrigValue<-Sys.getlocale("LC_COLLATE")
	on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
	Sys.setlocale("LC_COLLATE","C")
	logDebug("Changing LC_COLLATE to C for duration of run")
	checkPackageSettings()
	###logDebug("batchcorrForBatchType_writeImage theResultsBatchList ", theResultsBatchList)
	###logDebug("batchcorrForBatchType_writeImage theTitle ", theTitle)
	###logDebug("batchcorrForBatchType_writeImage theOutputDir ", theOutputDir)
	###logDebug("batchcorrForBatchType_writeImage theBatchCorrBase ", theBatchCorrBase)
	for (batchTypeName in names(theResultsBatchList))
	{
		###logDebug("batchcorrForBatchType_writeImage batchTypeName ", batchTypeName)
		resultsBatchArrayList<-theResultsBatchList[[batchTypeName]]
		if(TRUE==is.character(resultsBatchArrayList))
		{
			###logDebug("batchcorrForBatchType_writeImage character")
			### write string to file name
			outputDir <- checkCreateDir(checkCreateDir(checkCreateDir(theOutputDir, theBatchCorrBase), batchTypeName), "OneToMany")
			metricFile <- createDirPlusFilename(outputDir, theBatchCorrBase, "_Metric_Diagram.png")
			pvalueFile <- createDirPlusFilename(outputDir, theBatchCorrBase, "_Pvalue_Diagram.png")
			writeStringToImage(metricFile, breakIntoTitle(paste(theTitle, batchTypeName, resultsBatchArrayList, sep=" "), theOldChar=" ", theNewChar=" " ))
			writeStringToImage(pvalueFile, breakIntoTitle(paste(theTitle, batchTypeName, resultsBatchArrayList, sep=" "), theOldChar=" ", theNewChar=" " ))
		}
		else
		{
			###logDebug("batchcorrForBatchType_writeImage list")
			valueTypeList <- dimnames(resultsBatchArrayList)[[3]]
			czTypeList <- dimnames(resultsBatchArrayList)[[4]]
			for (valueType in valueTypeList)
			{
				###logDebug("batchcorrForBatchType_writeImage valueType ", valueType)
				for (czType in czTypeList)
				{
					###logDebug("batchcorrForBatchType_writeImage czType ", czType)
					resultBatch <- resultsBatchArrayList[,,valueType, czType]
					longBatchType <- paste(batchTypeName, valueType, czType, sep=" ")
					###logDebug("batchcorrForBatchType_writeImage array")
					stopifnotWithLogging("Batch result is not an array as expected.", is.array(resultBatch))
					###logDebug("batchcorrForBatchType_writeImage batch type ", longBatchType)
					myTitle <- paste(theTitle, longBatchType, sep=" ")
					### write string to file name
					outputDir <- checkCreateDir(checkCreateDir(checkCreateDir(theOutputDir, theBatchCorrBase), batchTypeName), "ManyToMany")
					###logDebug("batchcorrForBatchType_writeImage outputDir ", outputDir)
					metricFile <- createDirPlusFilename(outputDir, theBatchCorrBase, "_Metric", "_", czType, "_Diagram.png")
					pvalueFile <- createDirPlusFilename(outputDir, theBatchCorrBase, "_Pvalue", "_", czType, "_Diagram.png")
					###logDebug("batchcorrForBatchType_writeImage metricFile ", metricFile)
					###logDebug("batchcorrForBatchType_writeImage pvalueFile ", pvalueFile)
					writePvalueImage(metricFile, pvalueFile, batchTypeName, tolower(valueType), resultBatch, myTitle)
				}
			}
		}
	}
	return(NULL)
}

####################################################################
####################################################################
