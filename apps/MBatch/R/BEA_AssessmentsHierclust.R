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

####################################################################
####################################################################

createBatchEffectsOutput_hierclust<-function(theMatrixGeneData, theDataframeBatchData, theTitle,
                                             theDataVersion, theTestVersion,
																						 theHierClustOutputDir, theHierClustFileBase)
{
  # theHierClustOutputDir=getwd(), theHierClustFileBase="HierarchicalClustering")
	logDebug("hierclust_outputGraphics")
	collateOrigValue<-Sys.getlocale("LC_COLLATE")
	on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
	Sys.setlocale("LC_COLLATE","C")
	logDebug("Changing LC_COLLATE to C for duration of run")
	checkPackageSettings()
	stopifnotWithLogging("The number of columns of gene data must equal the number of rows of batch data.", ncol(theMatrixGeneData)==nrow(theDataframeBatchData))
	uDendSamples <- hierClust_calc(theMatrixGeneData)
	uDendFeatures <- hierClust_calc(t(theMatrixGeneData))
	hierClustOutputDir <- checkCreateDir(theHierClustOutputDir, theHierClustFileBase)
	hierClustOutputDir <- addVersionsIfNeeded(hierClustOutputDir, theDataVersion, theTestVersion)
	hierClustOutputDir <- checkDirForCreation(hierClustOutputDir)
	rdataFileSamples <- NULL
	rdataFileFeatures <- NULL
	newOutDir <- hierClustOutputDir
	if (is.null(uDendSamples))
	{
	  logDebug("createBatchEffectsOutput_hierclust hc_openAndWriteIssuesLogFile")
	  hc_openAndWriteIssuesLogFile(hierClustOutputDir)
	}
	else
	{
		logDebug("createBatchEffectsOutput_hierclust RData output")
		### RData output here, write uDend and theDataframeBatchData
		### D3 Hierarchical Clustering output here
	  rdataFileSamples <- writeHCDataTSVs(uDendSamples, hierClustOutputDir,
	                                      "HCData.tsv", "HCOrder.tsv", "uDend.RData")
	  rdataFileFeatures <-writeHCDataTSVs(uDendFeatures, hierClustOutputDir,
	                                      "HCData_feature.tsv", "HCOrder_feature.tsv", "uDend_feature.RData")
	  #logDebug("createBatchEffectsOutput_hierclust: theHierClustOutputDir=",theHierClustOutputDir, "theTitle=", theTitle, "theHierClustFileBase=", theHierClustFileBase)
		hierClust_draw(uDendSamples, theDataframeBatchData, theTitle, hierClustOutputDir, theHierClustFileBase)
	}
  list(newOutDir, rdataFileSamples, rdataFileFeatures)
}


hc_openAndWriteIssuesLogFile<-function(theOutputDir)
{
  myFile <- file(cleanFilePath(theOutputDir, "error.log"), "w+")
  on.exit(close(myFile))
  cat("HC not calculated -- too many nulls, NAs, infinities, or no deviation\n", file=myFile, append=TRUE)
}

####################################################################
####################################################################

hierClust_draw<-function(theUdend, theDataframeBatchData, theTitle,
												 theHierClustOutputDir=getwd(), theHierClustFileBase="HierarchicalClustering")
{
	logDebug("hierclust_draw")
	collateOrigValue<-Sys.getlocale("LC_COLLATE")
	on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
	Sys.setlocale("LC_COLLATE","C")
	logDebug("Changing LC_COLLATE to C for duration of run")
	checkPackageSettings()
	dataframeBatchDataWithMatchedColors <- matchColors(theDataframeBatchData)
	myColors <- getListOfColorsForHierClust(dataframeBatchDataWithMatchedColors, TRUE)
	matrixOfColors<-matrixWithIssues(character(0), nrow(dataframeBatchDataWithMatchedColors),0)
	for(batchName in c(2:length(dataframeBatchDataWithMatchedColors)))
	{
		selectedIndexes <- factor(dataframeBatchDataWithMatchedColors[, batchName])
		selectedColors <- myColors[selectedIndexes]
		matrixOfColors <- cbind(matrixOfColors, selectedColors)
	}
	#logDebug("hierClust_draw parameters: theHierClustOutputDir=", theHierClustOutputDir,"theTitle=",theTitle)
	writeClusteringImage(theHierClustOutputDir, theHierClustFileBase, theUdend, dataframeBatchDataWithMatchedColors, theTitle, matrixOfColors)
	writeClusteringLegends(theHierClustOutputDir, theHierClustFileBase, theUdend, dataframeBatchDataWithMatchedColors, theTitle, matrixOfColors)
	return(NULL)
}

hierClust_calc<-function(theMatrixGeneData)
{
	logDebug("hierclust_calc")
	collateOrigValue<-Sys.getlocale("LC_COLLATE")
	on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
	Sys.setlocale("LC_COLLATE","C")
	logDebug("Changing LC_COLLATE to C for duration of run")
	checkPackageSettings()
	###Do hierarchical clustering
	# do not need to check for NaN, as NaN is an NA
	if (is.null(theMatrixGeneData))
	{
		logWarn("Unable to calculate distanceMatrix--only one dimension")
		return(NULL)
	}
	else
	{
		subMatrix <- theMatrixGeneData[!is.na(rowSums(theMatrixGeneData)),]
		rowCount <- nrow(subMatrix)
		if ( (length(rowCount) == 1) && (!is.na(rowCount)) && (rowCount>0) )
		{
		  subMatrix <- subMatrix[!is.infinite(rowSums(subMatrix)),]
		  rowCount <- nrow(subMatrix)
		  if ( (length(rowCount) == 1) && (!is.na(rowCount)) && (rowCount>0) )
  		{
  		  uDend <- NULL
  		  tryCatch(
  		    {
  		      logDebug("calculating HC")
  		      checkIfTestError()
  			    d<-distanceMatrix(subMatrix, metric="pearson")
  			    #logDebug("checking distance matrix results")
  			    #subMatrix <- d[!is.na(rowSums(d)),]
  			    #subMatrix <- subMatrix[!is.infinite(rowSums(subMatrix)),]
  				  uDend<-hclust(d, method="ward.D2")
  				},
  				warning=function(e)
  				{
  				  logWarn("1 Unable to calculate hclust--too many NAs, Infinities or NaNs in data")
  				  uDend <- NULL
  				},
  				error=function(e)
  				{
  				  logWarn("2 Unable to calculate hclust--too many NAs, Infinities or NaNs in data")
  				  uDend <- NULL
  				})
  			return(uDend)
  		}
		  else
		  {
		    logWarn("3 Unable to calculate distanceMatrix--too many NAs, Infinities or NaNs in data")
		    return(NULL)
		  }
		}
		else
		{
			logWarn("4 Unable to calculate distanceMatrix--too many NAs, Infinities or NaNs in data")
			return(NULL)
		}
	}
}

####################################################################
####################################################################

getListOfColorsForHierClust<-function(theDataframeBatchData, theShuffle=FALSE)
{
	numberOfColors <- 2
	if (1==length(theDataframeBatchData))
	{
	  numberOfUniqueIds <- length(unique(as.vector(unlist(theDataframeBatchData[1]))))
	  if (numberOfUniqueIds>numberOfColors)
	  {
	    numberOfColors <- numberOfUniqueIds
	  }
	}
	else
	{
  	for(batchTypeIndex in c(2:length(theDataframeBatchData)))
  	{
  		numberOfUniqueIds <- length(unique(as.vector(unlist(theDataframeBatchData[batchTypeIndex]))))
  		if (numberOfUniqueIds>numberOfColors)
  		{
  			numberOfColors <- numberOfUniqueIds
  		}
  	}
	}
	myColors <- beaRainbow(numberOfColors, v=1.0, shuffle=theShuffle)
	return(myColors)
}

writeClusteringImage<-function(theHierClustOutputDir, theHierClustFileBase, uDend, theDataframeBatchData, theTitle, theMatrixOfColors)
{
	filename<-makeHCFileName_PNG(theHierClustOutputDir, "Diagram")
	logDebug("writeClusteringImage", filename)
	CairoPNG(filename=filename, width = 1000, height = 1000, pointsize=24)
	###defaultPar <- par(no.readonly = TRUE)
	on.exit(dev.off(), add = TRUE)
	###foo<-layout(matrix(c(1,2), nrow=2), heights=c(98, 2))
	par(cex=0.55)
	par(cex.lab=0.55)
	par(cex.axis=0.55)
	labelCex <- 1.0
	labelColumns <- as.vector(unlist(names(theDataframeBatchData)[2:length(theDataframeBatchData)]))
	labelWidths <- sapply(labelColumns, function (individualLabel)
	{
		return(nchar(individualLabel))
	}, USE.NAMES=FALSE)
	maxWidth <- max(labelWidths)
	#logDebug("maxWidth ", maxWidth)
	if (maxWidth>30)
	{
		###logDebug("maxWidth>30")
		par(cex=0.35)
		par(cex.lab=0.35)
		par(cex.axis=0.35)
	}
	else if (maxWidth>15)
	{
		###logDebug("maxWidth>15")
		par(cex=0.45)
		par(cex.lab=0.45)
		par(cex.axis=0.45)
	}
	#changed 'main= theTitle' to main=hierClusteringImageTitle
	dendromat(uDend, theMatrixOfColors, main=theTitle, labCol=labelColumns, labRow="",
						cex.lab=labelCex, cex.axis=labelCex, cex=labelCex)
	###par(defaultPar)
	###mtext(paste("MBatch", packageDescription("MBatch")["Version"], sep=" "), side=1, adj=1, outer=FALSE, line=4, cex=0.50, bg="white")
}

writeClusteringLegends<-function(theHierClustOutputDir, theHierClustFileBase, h, theDataframeBatchData, theTitle, theMatrixOfColors)
{
	legendList <- writeSeparateLegends(theHierClustOutputDir, theHierClustFileBase, h, theDataframeBatchData, theTitle, theMatrixOfColors)
	writeCombinedLegendHC(theHierClustOutputDir, legendList, theTitle)
}

writeSeparateLegends<-function(theHierClustOutputDir, theHierClustFileBase, h, theDataframeBatchData, theMainTitle, theMatrixOfColors)
{
	listOfFiles <- NULL
	for(batchTypeIndex in c(2:length(theDataframeBatchData)))
	{
		###logDebug("plot new")
		listOfFiles <- c(listOfFiles,
										 writeIndividualLegend(theHierClustOutputDir, theHierClustFileBase, theDataframeBatchData, theMatrixOfColors, batchTypeIndex))
	}
	return(listOfFiles)
}

writeIndividualLegend<-function(theHierClustOutputDir, theHierClustFileBase,
																theDataframeBatchData, theMatrixOfColors, batchTypeIndex)
{
	batchTypeName <- names(theDataframeBatchData)[batchTypeIndex]
	filename<-makeHCFileName_PNG(theHierClustOutputDir, "Legend", batchTypeName)
	logDebug("writeIndividualLegend", filename)
	version <- paste("MBatch", packageDescription("MBatch")["Version"], sep=" ")
	###draw LegendOn ExistingPlot(theDataframeBatchData, theMatrixOfColors, batchTypeIndex)
	batchTypeName <- names(theDataframeBatchData)[batchTypeIndex]
	batchIdsListAll <- as.vector(unlist(theDataframeBatchData[batchTypeIndex]))
	batchIdsListUnique <- unique(batchIdsListAll)
	batchIdColors <- unique(as.vector(unlist(theMatrixOfColors[,batchTypeIndex-1])))
	sortedBatchIds <- sort(batchIdsListUnique)
	sortedColors <- as.vector(unlist(sortXBasedOnOrderOfY(batchIdColors, batchIdsListUnique)))
	stopifnotWithLogging("The number of sorted colors and number of sorted batch ids must match.", length(sortedBatchIds)==length(sortedColors))
	sortedBatchIdsWithCount <- as.vector(unlist(sapply(sortedBatchIds, function (batchId)
	{
		return(paste(batchId," (",sum(batchIdsListAll==batchId), ")", sep=""))
	}, USE.NAMES=FALSE)))
	##leg end("center", legend=sortedBatchIdsWithCount, fill=sortedColors, title=batchTypeName)
	mbatchStandardLegend(batchTypeName, version, sortedBatchIdsWithCount, sortedColors, NULL, filename)
	return(filename)
}

writeCombinedLegendHC<-function(theHierClustOutputDir, theLegendList, theMainTitle)
{
	filename<-makeHCFileName_PNG(theHierClustOutputDir, "Legend", "ALL")
	logDebug("writeCombinedLegendHC", filename)
	mbatchStandardCombineLegends(theMainTitle, filename, theLegendList)
}

matchColors<-function(theDataframeBatchData)
{
	theDataframeBatchData<-data.frame(theDataframeBatchData, check.names=FALSE)
	for(i in 1:ncol(theDataframeBatchData))
	{
		theDataframeBatchData[,i]<-factor(theDataframeBatchData[,i])
	}
	nl<-sapply(theDataframeBatchData, nlevels)
	ref.ind<-which(nl==max(nl))[1]
	Lev<-levels(theDataframeBatchData[,ref.ind])
	for(i in 1:ncol(theDataframeBatchData))
	{
		ori.lev<-levels(theDataframeBatchData[,i])
		temp.tab<-table(theDataframeBatchData[,ref.ind], theDataframeBatchData[,i])
		indx<-unlist(sapply(Lev, function(x) which(temp.tab[x,]==max(temp.tab[x,]))[1]))
		if(any(duplicated(indx)))
		{
			dup<-unique(indx[duplicated(indx)])
			for(ii in dup)
			{
				keep<-which(temp.tab[,ii]==max(temp.tab[,ii]))[1]
				indx[indx==ii]<-NA
				indx[keep]<-ii
			}
		}
		new.lev<-ori.lev[indx]
		nfill<-sum(is.na(new.lev))
		not.in<-ori.lev[!ori.lev %in% new.lev]
		if(nfill>0)
		{
			new.lev[is.na(new.lev)]<-paste('NAfill', 1:nfill, sep="")
		}
		new.lev<-c(new.lev, not.in)
		theDataframeBatchData[,i]<-factor(as.vector(theDataframeBatchData[,i]), levels=new.lev)
	}
	return(theDataframeBatchData)
}

####################################################################
####################################################################

makeHCFileName_PNG<-function(theDir, theDiagramOrLegend, theLegendType="")
{
	###stopifnotWithLogging("theCentroidOrPlain should be PCA-Plus or PcaPlain", (("PCA-Plus"==theCentroidOrPlain)||("PcaPlain"==theCentroidOrPlain)))
	###stopifnotWithLogging("thePlotType should be Many2Many|One2Many|DualBatch", (("Many2Many"==thePlotType)||("One2Many"==thePlotType)||("DualBatch"==thePlotType)))
	stopifnotWithLogging("theDiagramOrLegend should be Diagram|Legend", (("Diagram"==theDiagramOrLegend)||("Legend"==theDiagramOrLegend)))
	if ("Legend"!=theDiagramOrLegend)
	{
		if (""!=theLegendType)
		{
			stopWithLogging(paste("theLegendType is ", theLegendType, " which must be empty for theDiagramOrLegend=", theDiagramOrLegend, sep=""))
		}
	}
	if ("Legend"==theDiagramOrLegend)
	{
		createDirPlusFilename(theDir, "HierarchicalClustering_", theDiagramOrLegend, "-", theLegendType, ".png")
	}
	else
	{
		createDirPlusFilename(theDir, "HierarchicalClustering_", theDiagramOrLegend, ".png")
	}
}

####################################################################
####################################################################

writeHCDataTSVs<-function(uDend, theHierClustOutputDir, theOutputHCDataFileName, theOutputHCOrderFileName,
                          theUdendRData)
{
  rdataFile <- cleanFilePath(theHierClustOutputDir,theUdendRData)
  logInfo("writeHCDataTSVs rdataFile=", rdataFile)
	if (!is.null(uDend))
	{
	  data<-cbind(uDend$merge, uDend$height, deparse.level=0)
	  colnames(data)<-c("A", "B", "Height")
  	###Write out the data as a Tab separated file to the specified location
  	write.table(data, file = cleanFilePath(theHierClustOutputDir,theOutputHCDataFileName), append = FALSE, quote = FALSE, sep = "\t", row.names=FALSE)

  	data<-cbind(uDend$labels, uDend$order, deparse.level=0)
  	colnames(data)<-c("Id", "Order")
  	###Write out the order data as a Tab separated file to the specified location (1 more row than data file)
  	write.table(data, file = cleanFilePath(theHierClustOutputDir,theOutputHCOrderFileName), append = FALSE, quote = FALSE, sep = "\t", row.names=FALSE)
    # write udend RData file
  	save(uDend, file=rdataFile)
	}
	rdataFile
}

####################################################################
####################################################################
