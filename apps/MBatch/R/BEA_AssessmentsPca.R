# MBatch Copyright (c) 2011-2024 University of Texas MD Anderson Cancer Center
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
### top level functions
####################################################################

createBatchEffectsOutput_pca<-function(theMatrixGeneData, theDataframeBatchData, theTitle,
		theDoPlainMtoMFlag, theDoCentroidsMtoMFlag, theDoDSCFlag, theDoDscPermsFileFlag, theDoSampleLocatorFlag,
		theIsPcaTrendFunction, theListOfComponentsToPlot,
		theDSCPermutations, theDSCThreads, theMinBatchSize,
		theDataVersion, theTestVersion,
		theOutputDir, thePcaCentroidsBase="PCA-Plus", thePcaPlainBase="PCAPlain",
		theSeed=NULL, theGeneLimit=0)
		###theIsPcaTrendFunction=NULL, theListOfComponentsToPlot=c(1, 2),
		###theDSCPermutations=2000, theDSCThreads=1, theMinBatchSize=2,
{
	logDebug("createBatchEffectsOutput_pca")
  collateOrigValue<-Sys.getlocale("LC_COLLATE")
	on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
	Sys.setlocale("LC_COLLATE","C")
	logDebug("Changing LC_COLLATE to C for duration of run")
	checkPackageSettings()
	### code assumes that theMatrixGeneData and theDataframeBatchData contain
	### the same samples and the same number of samples and in the same order
	stopifnotWithLogging("The number of columns of gene data must equal the number of rows of batch data.", ncol(theMatrixGeneData)==nrow(theDataframeBatchData))
	stopifnotWithLogging("createBatchEffectsOutput_pca - must have > 0 rows in theMatrixGeneData.", nrow(theMatrixGeneData)>0)
	stopifnotWithLogging("createBatchEffectsOutput_pca - must have > 0 cols in theMatrixGeneData.", ncol(theMatrixGeneData)>0)
	dscOutputDir <- NULL
	vectorOfOutputDirs <- c()
	if (TRUE==theDoDSCFlag)
	{
	  # make directory for later use
	  #dscOutputDir <- addVersionsIfNeeded(cleanFilePath(theOutputDir, "DSC"), theDataVersion, theTestVersion)
	  #checkDirForCreation(dscOutputDir)
	  #file.create(file.path(dscOutputDir, "DSCOverview.tsv"))
	  #logDebug("make future dir", dscOutputDir)
	  # write complist file
	  # NOTE - MUST BE AT TOP LEVEL FOR DSC OVERVIEW FUNCTION TO WORK
	  logDebug("write ALL__CompListDSC.RData", theOutputDir)
		saveCompListDscData(cleanFilePath(theOutputDir, "ALL__CompListDSC.RData"), theListOfComponentsToPlot)
	}
	###logDebug("createBatchEffectsOutput_pca - rows in theMatrixGeneData=", nrow(theMatrixGeneData))
	###logDebug("createBatchEffectsOutput_pca - cols in theMatrixGeneData=", ncol(theMatrixGeneData))
	for(batchTypeIndex in c(2:length(theDataframeBatchData)))
	{
		### compile data and information for display
		batchTypeName <- names(theDataframeBatchData)[batchTypeIndex]
		batchIdsForSamples <- as.character(as.vector(unlist(theDataframeBatchData[batchTypeIndex])))
		###logDebug("createBatchEffectsOutput_pca - before sync length(batchIdsForSamples)=", length(batchIdsForSamples))
		pca <- doSamplePcaCall(theMatrixGeneData, theMinBatchSize, batchIdsForSamples, theListOfComponentsToPlot, theSeed, theGeneLimit)
		dscOutputDir <- checkCreateDir(checkCreateDir(theOutputDir, batchTypeName), "ManyToMany")
		dscOutputDir <- addVersionsIfNeeded(dscOutputDir, theDataVersion, theTestVersion)
		if (is.null(pca))
		{
		  checkDirForCreation(dscOutputDir)
		  openAndWriteIssuesLogFile(dscOutputDir)
		}
		else
		{
			###logDebug("createBatchEffectsOutput_pca - rows in pca@scores=", nrow(pca@scores))
			###logDebug("createBatchEffectsOutput_pca - cols in pca@scores=", ncol(pca@scores))
			batchIdsForSamples <- syncBatchIdsAndPcaResults(batchIdsForSamples, rownames(pca@scores), theDataframeBatchData$Sample)
			###logDebug("createBatchEffectsOutput_pca - after sync length(rownames(pca@scores))=", length(rownames(pca@scores)))
			###logDebug("createBatchEffectsOutput_pca - after sync length(batchIdsForSamples)=", length(batchIdsForSamples))
			###logDebug("createBatchEffectsOutput_pca - after sync rows in pca@scores=", nrow(pca@scores))
			###		logDebug("createBatchEffectsOutput_pca - after sync cols in pca@scores=", ncol(pca@scores))
			stopifnotWithLogging("createBatchEffectsOutput_pca - Number of rownames for pca score (sample ids) should match length of batch ids for samples", (length(rownames(pca@scores))==length(batchIdsForSamples)))
			sortedListOfBatchIds <- sortBatchesBasedOnSize(unique(batchIdsForSamples), batchIdsForSamples)
			###logDebug("sortBatchesBasedOnSize(sortedListOfBatchIds) = ", paste(sortedListOfBatchIds, collapse=","))
			batchTypeOutputDir <- checkCreateDir(checkCreateDir(theOutputDir, batchTypeName), "ManyToMany")
			batchTypeOutputDir <- addVersionsIfNeeded(batchTypeOutputDir, theDataVersion, theTestVersion)
			vectorOfOutputDirs <- c(vectorOfOutputDirs, batchTypeOutputDir)
			checkDirForCreation(batchTypeOutputDir)
			centroidsOutputDir <- ""
			if((TRUE==theDoCentroidsMtoMFlag)||(TRUE==theDoSampleLocatorFlag))
			{
				centroidsOutputDir <- checkCreateDir(batchTypeOutputDir, thePcaCentroidsBase)
				checkDirForCreation(centroidsOutputDir)
			}
			plainOutputDir <- ""
			if(TRUE==theDoPlainMtoMFlag)
			{
			  plainOutputDir <- checkCreateDir(batchTypeOutputDir, thePcaPlainBase)
				checkDirForCreation(plainOutputDir)
			}
			dscAllResults <- NULL
			if(TRUE==theDoDSCFlag)
			{
				logDebug("call pvalueDSC")
				dscAllResults <- pvalueDSC(pca, batchIdsForSamples, theDSCPermutations, 0, 0, theDSCThreads, theSeed)
				###logDebug("openAndWriteDscAllFile")
				openAndWriteDscAllFile(pca, dscAllResults, batchTypeOutputDir, "ANY", "NA", "NA")
				if (TRUE==theDoDscPermsFileFlag)
				{
					openAndWriteDscPermsFile(dscAllResults, batchTypeOutputDir, "ANY", "NA", "NA")
				}
				###logDebug("after openAndWriteDscAllFile")
			}
			logDebug("Write PCA Files 1 ", batchTypeOutputDir)
			mytitle <- paste(theTitle, "/", "PCA+", batchTypeName, sep=" ")
			writePcaDataFilesForDataset(batchTypeOutputDir, theDSCPermutations, pca, dscAllResults, theListOfComponentsToPlot, mytitle)
			####
			logDebug("write writeSharedFveWeightScoresFiles")
			writeSharedFveWeightScoresFiles(theDoSampleLocatorFlag, batchTypeOutputDir, pca, theSampleIds=rownames(theMatrixGeneData), theGeneIds=colnames(theMatrixGeneData))
			logDebug("createBatchEffectsOutput_pca loop")
			###logDebug("createBatchEffectsOutput_pca():", centroidsOutputDir, " ", plainOutputDir)
			for(i in seq(from=1, to=(length(theListOfComponentsToPlot)-1), by=2 ))
			{
				componentA <- theListOfComponentsToPlot[[i]]
				componentB <- theListOfComponentsToPlot[[i+1]]
				###
				dscTxtFile <- makePcaFileName_TXT(batchTypeOutputDir, "ANY", componentA, componentB, "DSC")
				errorFile <- file.path(batchTypeOutputDir, "error.log")
				###
				tempCentroids <- makePcaFileName_PNG(centroidsOutputDir, "ALL", componentA, componentB, "Diagram")
				centrLegendAll <- makePcaFileName_PNG(centroidsOutputDir, "ALL", componentA, componentB, "Legend", "ALL")
				centrLegendDSC <- makePcaFileName_PNG(centroidsOutputDir, "ALL", componentA, componentB, "Legend", "DSC")
				centrLegendRays <- makePcaFileName_PNG(centroidsOutputDir, "ALL", componentA, componentB, "Legend", "Rays")
				centrLegendPoints <- makePcaFileName_PNG(centroidsOutputDir, "ALL", componentA, componentB, "Legend", "Points")
				centroidsImageTitle <- paste(theTitle, "ALL", componentA, componentB)
				###
				tempPlain <- makePcaFileName_PNG(plainOutputDir, "ALL", componentA, componentB, "Diagram")
				plainLegendAll <- makePcaFileName_PNG(plainOutputDir, "ALL", componentA, componentB, "Legend", "ALL")
				plainLegendDSC <- makePcaFileName_PNG(plainOutputDir, "ALL", componentA, componentB, "Legend", "DSC")
				plainLegendRays <- makePcaFileName_PNG(plainOutputDir, "ALL", componentA, componentB, "Legend", "Rays")
				plainLegendPoints <- makePcaFileName_PNG(plainOutputDir, "ALL", componentA, componentB, "Legend", "Points")
				plainImageTitle <- paste(theTitle, "ALL", componentA, componentB)
				###
				### TODO: MAny2Many DSC objects returned from doInternalPCa need to be writen to Calc.RData file
				doInternalPca(pca, dscAllResults, componentA, componentB, theMatrixGeneData,
				    theDoPlainMtoMFlag, theDoCentroidsMtoMFlag, theDoDSCFlag,
						theDoDscPermsFileFlag, tempCentroids, tempPlain, dscTxtFile,
						centrLegendAll, centrLegendDSC, centrLegendRays, centrLegendPoints,
						plainLegendAll, plainLegendDSC, plainLegendRays, plainLegendPoints,
						centroidsImageTitle, plainImageTitle,
						batchIdsForSamples, NULL, sortedListOfBatchIds, NULL, batchTypeName, NULL,
						theIsPcaTrendFunction, theDSCPermutations, theDSCThreads, theMinBatchSize, colnames(theMatrixGeneData),
						errorFile, theSeed, dscOutputDir)
			}
		}
	}
	vectorOfOutputDirs
}

even <- function(theVal)
{
  return ((theVal%%2) == 0)
}

createBatchEffectsOutput_pca_dualBatch<-function(theMatrixGeneData, theDataframeBatchData,
    theListForDoCentroidDualBatchType,
		theTitle, theDoDSCFlag, theDoDscPermsFileFlag, theDoSampleLocatorFlag,
		theIsPcaTrendFunction, theListOfComponentsToPlot,
		theDSCPermutations, theDSCThreads, theMinBatchSize,
		theOutputDir, theDataVersion, theTestVersion,
		theFileBase="PCA-Plus",
		theSeed=NULL, theGeneLimit=0)
		###theIsPcaTrendFunction=NULL, theListOfComponentsToPlot=c(1, 2),
		###theDSCPermutations=2000, theDSCThreads=1, theMinBatchSize=2,
{
	stopifnotWithLogging("createBatchEffectsOutput_pca_dualBatch - The number of columns of gene data must equal the number of rows of batch data.", ncol(theMatrixGeneData)==nrow(theDataframeBatchData))
	stopifnotWithLogging("createBatchEffectsOutput_pca_dualBatch - Must have more than one batch type to plot.", length(theListForDoCentroidDualBatchType)>1)
	stopifnotWithLogging("createBatchEffectsOutput_pca_dualBatch - Must have an even number of batch types to plot.", even(length(theListForDoCentroidDualBatchType)))
	stopifnotWithLogging("createBatchEffectsOutput_pca_dualBatch - must have > 0 rows in theMatrixGeneData.", nrow(theMatrixGeneData)>0)
	stopifnotWithLogging("createBatchEffectsOutput_pca_dualBatch - must have > 0 cols in theMatrixGeneData.", ncol(theMatrixGeneData)>0)
	includedBatchTypes <- as.vector(unlist(names(theDataframeBatchData)[2:length(theDataframeBatchData)]))
	for (xx in 1:length(theListForDoCentroidDualBatchType))
	{
		stopifnotWithLogging(paste("createBatchEffectsOutput_pca_dualBatch - All batch types in theListForDoCentroidDualBatchType must match types available in the data.",  theListForDoCentroidDualBatchType[[xx]], " is invalid.", sep=" "), (theListForDoCentroidDualBatchType[[xx]] %in% includedBatchTypes)>0)
	}
	collateOrigValue<-Sys.getlocale("LC_COLLATE")
	on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
	Sys.setlocale("LC_COLLATE","C")
	logDebug("Changing LC_COLLATE to C for duration of run")
	checkPackageSettings()
	dscOutputDir <- NULL
	vectorOfOutputDirs <- c()
	if (TRUE==theDoDSCFlag)
	{
	  # make directory for later user
	  dscOutputDir <- addVersionsIfNeeded(cleanFilePath(theOutputDir, "DSC"), theDataVersion, theTestVersion)
	  checkDirForCreation(dscOutputDir)
	  file.create(file.path(dscOutputDir, "DSCOverview.tsv"))
	  logDebug("make future dir", dscOutputDir)
	  # write complist file
	  # NOTE - MUST BE AT TOP LEVEL FOR DSC OVERVIEW FUNCTION TO WORK
	  logDebug("write ALL__CompListDSC.RData", theOutputDir)
	  saveCompListDscData(cleanFilePath(theOutputDir, "ALL__CompListDSC.RData"), theListOfComponentsToPlot)
	}
	logDebug("createBatchEffectsOutput_pca_dualBatch - start loop")
	for(x in seq(from=1, to=(length(theListForDoCentroidDualBatchType)-1), by=2 ))
	{
		batchTypeA <- theListForDoCentroidDualBatchType[[x]]
		batchTypeB <- theListForDoCentroidDualBatchType[[x+1]]
		### compile data and information for display
		batchIdsForSamplesA <- as.vector(unlist(theDataframeBatchData[batchTypeA]))
		batchIdsForSamplesB <- as.vector(unlist(theDataframeBatchData[batchTypeB]))
		###logDebug("createBatchEffectsOutput_pca_dualBatch - loop doSamplePcaCall")
		pca <- doSamplePcaCall(theMatrixGeneData, theMinBatchSize, batchIdsForSamplesA, theListOfComponentsToPlot, theSeed, theGeneLimit)
		if (!is.null(pca))
		{
			batchIdsForSamplesA <- syncBatchIdsAndPcaResults(batchIdsForSamplesA, rownames(pca@scores), theDataframeBatchData$Sample)
			batchIdsForSamplesB <- syncBatchIdsAndPcaResults(batchIdsForSamplesB, rownames(pca@scores), theDataframeBatchData$Sample)
			stopifnotWithLogging("createBatchEffectsOutput_pca_dualBatch A - Number of rownames for pca score (sample ids) should match length of batch ids for samples", (length(rownames(pca@scores))==length(batchIdsForSamplesA)))
			stopifnotWithLogging("createBatchEffectsOutput_pca_dualBatch B - Number of rownames for pca score (sample ids) should match length of batch ids for samples", (length(rownames(pca@scores))==length(batchIdsForSamplesB)))
			sortedListOfBatchIdsA <- sort(unique(batchIdsForSamplesA))
			sortedListOfBatchIdsB <- sort(unique(batchIdsForSamplesB))
			###logDebug("createBatchEffectsOutput_pca_dualBatch - loop pvalueDSC")
			tempBatchTypes <- paste(batchTypeA, "with", batchTypeB, sep="")
			outputDir <- checkCreateDir(checkCreateDir(theOutputDir, tempBatchTypes), "DualBatch")
			outputDir <- addVersionsIfNeeded(outputDir, theDataVersion, theTestVersion)
			outputDir <- checkCreateDir(outputDir, theFileBase)
			checkDirForCreation(outputDir)
			vectorOfOutputDirs <- c(vectorOfOutputDirs, outputDir)
			datafilesOutputDir <- checkCreateDir(checkCreateDir(theOutputDir, tempBatchTypes), "DualBatch")
			datafilesOutputDir <- addVersionsIfNeeded(datafilesOutputDir, theDataVersion, theTestVersion)
			checkDirForCreation(datafilesOutputDir)
			dscAllResults <- NULL
			if(TRUE==theDoDSCFlag)
			{
				dscAllResults <- pvalueDSC(pca, batchIdsForSamplesA, theDSCPermutations, 0, 0, theDSCThreads, theSeed)
				openAndWriteDscAllFile(pca, dscAllResults, datafilesOutputDir, "ANY", "NA", "NA")
				if(TRUE==theDoDscPermsFileFlag)
				{
					openAndWriteDscPermsFile(dscAllResults, datafilesOutputDir, "ANY", "NA", "NA")
				}
			}
			logDebug("Write PCA Files 2 ", datafilesOutputDir)
			mytitle <- paste(theTitle, "/", "PCA+", batchTypeA, "vs", batchTypeB, sep=" ")
			writePcaDataFilesForDataset(datafilesOutputDir, theDSCPermutations, pca, dscAllResults, theListOfComponentsToPlot, mytitle)
			####
			writeSharedFveWeightScoresFiles(theDoSampleLocatorFlag, datafilesOutputDir, pca, theSampleIds=rownames(theMatrixGeneData), theGeneIds=colnames(theMatrixGeneData))
			for(i in seq(from=1, to=(length(theListOfComponentsToPlot)-1), by=2 ))
			{
				componentA <- theListOfComponentsToPlot[[i]]
				componentB <- theListOfComponentsToPlot[[i+1]]
				###
				dscTxtFile <- makePcaFileName_TXT(datafilesOutputDir, "ANY", componentA, componentB, "DSC")
				errorFile <- file.path(outputDir, "error.log")
				###
				tempCentroids <- makePcaFileName_PNG(outputDir, tempBatchTypes, componentA, componentB, "Diagram")
				centrLegendAll <- makePcaFileName_PNG(outputDir, tempBatchTypes, componentA, componentB, "Legend", "ALL")
				centrLegendDSC <- makePcaFileName_PNG(outputDir, tempBatchTypes, componentA, componentB, "Legend", "DSC")
				centrLegendRays <- makePcaFileName_PNG(outputDir, tempBatchTypes, componentA, componentB, "Legend", "Rays")
				centrLegendPoints <- makePcaFileName_PNG(outputDir, tempBatchTypes, componentA, componentB, "Legend", "Points")
				centroidsImageTitle <- paste(theTitle, tempBatchTypes, componentA, componentB)
				###
				errorFilePlain <- ""
				tempPlain <- ""
				plainLegendAll <- ""
				plainLegendDSC <- ""
				plainLegendRays <- ""
				plainLegendPoints <- ""
				plainImageTitle <- ""
				###
				### TODO: DualBatch DSC objects returned from doInternalPCa need to be writen to Calc.RData file
				doInternalPca(pca, dscAllResults, componentA, componentB, theMatrixGeneData, FALSE, TRUE, theDoDSCFlag, theDoDscPermsFileFlag,
						tempCentroids, tempPlain, dscTxtFile,
						centrLegendAll, centrLegendDSC, centrLegendRays, centrLegendPoints,
						plainLegendAll, plainLegendDSC, plainLegendRays, plainLegendPoints,
						centroidsImageTitle, plainImageTitle,
						batchIdsForSamplesA, batchIdsForSamplesB, sortedListOfBatchIdsA, sortedListOfBatchIdsB, batchTypeA, batchTypeB,
						theIsPcaTrendFunction, theDSCPermutations, theDSCThreads, theMinBatchSize, colnames(theMatrixGeneData),
						errorFile, theSeed, dscOutputDir)
			}
		}
	}
	return(NULL)
}

createBatchEffectsOutput_pca_one2many<-function(theMatrixGeneData, theDataframeBatchData,
		thePcaCentroidsTitle, thePcaPlainTitle,
		theDoPlainOtoMFlag, theDoCentroidsOtoMFlag, theDoDSCFlag, theDoDscPermsFileFlag, theDoSampleLocatorFlag,
		theIsPcaTrendFunction, theListOfComponentsToPlot,
		theDSCPermutations, theDSCThreads, theMinBatchSize,
		theOutputDir, theDataVersion, theTestVersion,
		thePcaCentroidsBase="PCA-Plus", thePcaPlainBase="PCAPlain",
		theSeed=NULL, theGeneLimit=0)
{
  ###theIsPcaTrendFunction=NULL, theListOfComponentsToPlot=c(1, 2),
  ###theDSCPermutations=2000, theDSCThreads=1, theMinBatchSize=2,
	logDebug("createBatchEffectsOutput_pca_one2many")
  collateOrigValue<-Sys.getlocale("LC_COLLATE")
	on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
	Sys.setlocale("LC_COLLATE","C")
	logDebug("Changing LC_COLLATE to C for duration of run")
	checkPackageSettings()
	stopifnotWithLogging("The number of columns of gene data must equal the number of rows of batch data.", ncol(theMatrixGeneData)==nrow(theDataframeBatchData))
	stopifnotWithLogging("createBatchEffectsOutput_pca_one2many - must have > 0 rows in theMatrixGeneData.", nrow(theMatrixGeneData)>0)
	stopifnotWithLogging("createBatchEffectsOutput_pca_one2many - must have > 0 cols in theMatrixGeneData.", ncol(theMatrixGeneData)>0)
	dscOutputDir <- NULL
	if (TRUE==theDoDSCFlag)
	{
	  # make directory for later user
	  dscOutputDir <- addVersionsIfNeeded(cleanFilePath(theOutputDir, "OneToMany-DSC"), theDataVersion, theTestVersion)
	  checkDirForCreation(dscOutputDir)
	  file.create(file.path(dscOutputDir, "DSCOverview.tsv"))
	  logDebug("make future dir", dscOutputDir)
	  # write complist file
	  # NOTE - MUST BE AT TOP LEVEL FOR DSC OVERVIEW FUNCTION TO WORK
	  logDebug("write ALL__CompListDSC.RData", theOutputDir)
	  saveCompListDscData(cleanFilePath(theOutputDir, "ALL__CompListDSC.RData"), theListOfComponentsToPlot)
	}
	for(batchTypeIndex in c(2:length(theDataframeBatchData)))
	{
		batchTypeName <- names(theDataframeBatchData)[batchTypeIndex]
		batchIdsForSamples <- as.vector(unlist(theDataframeBatchData[batchTypeIndex]))
		uniqueListOfBatchIds <- unique(batchIdsForSamples)
		for(batchId in uniqueListOfBatchIds)
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
			o2mBatchIdsForSamples <- as.vector(unlist(oneToManyDataframe[2]))
			pca <- doSamplePcaCall(theMatrixGeneData, theMinBatchSize, o2mBatchIdsForSamples, theListOfComponentsToPlot, theSeed, theGeneLimit)
			if (!is.null(pca))
			{
				o2mBatchIdsForSamples <- syncBatchIdsAndPcaResults(o2mBatchIdsForSamples, rownames(pca@scores), theDataframeBatchData$Sample)
				stopifnotWithLogging("createBatchEffectsOutput_pca_one2many - Number of rownames for pca score (sample ids) should match length of batch ids for samples", (length(rownames(pca@scores))==length(o2mBatchIdsForSamples)))
				o2mUniqueListOfBatchIds <- sortBatchesBasedOnSize(unique(o2mBatchIdsForSamples), o2mBatchIdsForSamples)
				batchIdentifiersWithCount <- getListOfBatchNamesWithCounts(o2mBatchIdsForSamples, c(batchId, "Other Batches"))
				batchIdFilename <- compressIntoFilename(batchId)
				batchTypeOutputDir <- cleanFilePath(cleanFilePath(theOutputDir, batchTypeName), paste("OneToMany", batchId, sep="-"))
				batchTypeOutputDir <- addVersionsIfNeeded(batchTypeOutputDir, theDataVersion, theTestVersion)
				checkDirForCreation(batchTypeOutputDir)
				centroidsOutputDir <- ""
				if((TRUE==theDoCentroidsOtoMFlag)||(TRUE==theDoSampleLocatorFlag))
				{
				  centroidsOutputDir <- cleanFilePath(cleanFilePath(cleanFilePath(theOutputDir, batchTypeName), paste("OneToMany", batchId, sep="-")), thePcaCentroidsBase)
				  centroidsOutputDir <- addVersionsIfNeeded(batchTypeOutputDir, theDataVersion, theTestVersion)
				  checkDirForCreation(centroidsOutputDir)
				}
				plainOutputDir <- ""
				if(TRUE==theDoPlainOtoMFlag)
				{
					plainOutputDir <- cleanFilePath(cleanFilePath(cleanFilePath(theOutputDir, batchTypeName), paste("OneToMany", batchId, sep="-")), thePcaPlainBase)
					plainOutputDir <- addVersionsIfNeeded(batchTypeOutputDir, theDataVersion, theTestVersion)
					checkDirForCreation(plainOutputDir)
				}
				pcaTitle <- paste(batchIdentifiersWithCount, sep="", collapse=" versus ")
				dscAllResults <- NULL
				if(TRUE==theDoDSCFlag)
				{
					dscAllResults <- pvalueDSC(pca, o2mBatchIdsForSamples, theDSCPermutations, 0, 0, theDSCThreads, theSeed)
					openAndWriteDscAllFile(pca, dscAllResults, batchTypeOutputDir, "ANY", "NA", "NA")
					if(TRUE==theDoDscPermsFileFlag)
					{
						openAndWriteDscPermsFile(dscAllResults, batchTypeOutputDir, "ANY", "NA", "NA")
					}
				}
				### RData output here for one2many, write pca (scores has Principle Component and Sample ids), dscAllResults and batchTypeName and batchIdsForSamples
				logDebug("Write PCA Files 3 ", batchTypeOutputDir)
				mytitle <- paste(thePcaPlainTitle, "/", "PCA+", batchTypeName, sep=" ")
				writePcaDataFilesForDataset(batchTypeOutputDir, theDSCPermutations, pca, dscAllResults, theListOfComponentsToPlot)
				####
				writeSharedFveWeightScoresFiles(theDoSampleLocatorFlag, batchTypeOutputDir, pca, theSampleIds=rownames(theMatrixGeneData), theGeneIds=colnames(theMatrixGeneData))
				for(i in seq(from=1, to=(length(theListOfComponentsToPlot)-1), by=2 ))
				{
					componentA <- theListOfComponentsToPlot[[i]]
					componentB <- theListOfComponentsToPlot[[i+1]]
					###
					dscTxtFile <- makePcaFileName_TXT(batchTypeOutputDir, "ANY", componentA, componentB, "DSC")
					errorFile <- file.path(batchTypeOutputDir, "error.log")
					###
					tempCentroids <- makePcaFileName_PNG(centroidsOutputDir, batchIdFilename, componentA, componentB, "Diagram")
					centrLegendAll <- makePcaFileName_PNG(centroidsOutputDir, batchIdFilename, componentA, componentB, "Legend", "ALL")
					centrLegendDSC <- makePcaFileName_PNG(centroidsOutputDir, batchIdFilename, componentA, componentB, "Legend", "DSC")
					centrLegendRays <- makePcaFileName_PNG(centroidsOutputDir, batchIdFilename, componentA, componentB, "Legend", "Rays")
					centrLegendPoints <- makePcaFileName_PNG(centroidsOutputDir, batchIdFilename, componentA, componentB, "Legend", "Points")
					centroidsImageTitle <- paste(thePcaCentroidsTitle, batchIdFilename, componentA, componentB)
					###
					tempPlain <- makePcaFileName_PNG(plainOutputDir, batchIdFilename, componentA, componentB, "Diagram")
					plainLegendAll <- makePcaFileName_PNG(plainOutputDir, batchIdFilename, componentA, componentB, "Legend", "ALL")
					plainLegendDSC <- makePcaFileName_PNG(plainOutputDir, batchIdFilename, componentA, componentB, "Legend", "DSC")
					plainLegendRays <- makePcaFileName_PNG(plainOutputDir, batchIdFilename, componentA, componentB, "Legend", "Rays")
					plainLegendPoints <- makePcaFileName_PNG(plainOutputDir, batchIdFilename, componentA, componentB, "Legend", "Points")
					plainImageTitle <- paste(thePcaPlainTitle, batchIdFilename, componentA, componentB)
					###
					### TODO: One2Many DSC objects returned from doInternalPCa need to be writen to Calc.RData file
					doInternalPca(pca, dscAllResults, componentA, componentB, theMatrixGeneData, theDoPlainOtoMFlag, theDoCentroidsOtoMFlag, theDoDSCFlag,
							theDoDscPermsFileFlag, tempCentroids, tempPlain, dscTxtFile,
							centrLegendAll, centrLegendDSC, centrLegendRays, centrLegendPoints,
							plainLegendAll, plainLegendDSC, plainLegendRays, plainLegendPoints,
							centroidsImageTitle, plainImageTitle,
							o2mBatchIdsForSamples, NULL, o2mUniqueListOfBatchIds, NULL, batchTypeName, NULL,
							theIsPcaTrendFunction, theDSCPermutations, theDSCThreads, theMinBatchSize, colnames(theMatrixGeneData),
							errorFile, theSeed, dscOutputDir)
				}
			}
		}
	}
	return(NULL)
}

####################################################################
### shared functions
####################################################################

syncBatchIdsAndPcaResults<-function(theBatchIdsForSamples, thePcaSamples, theOriginalSamples)
{
	if (length(thePcaSamples)==length(theOriginalSamples))
	{
		return(theBatchIdsForSamples)
	}
	else
	{
		return(theBatchIdsForSamples[is.element(theOriginalSamples, thePcaSamples)])
	}
}

getPchList<-function(theLength=0)
{
	###return(c(0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15:19))
	###pchList <- 0:17
	pchList <- 0:14
	if(0!=theLength)
	{
		pchList <- rep(pchList, length.out=theLength)
	}
	return(pchList)
}

calculatePch<-function(theIndex, thePchList)
{
	retVal <- 1
	if (theIndex<(length(thePchList)+1))
	{
		retVal <- theIndex
	}
	else
	{
		retVal <- ((theIndex) %% (length(thePchList)))
	}
	if (0==retVal)
	{
		retVal <- length(thePchList)
	}
	return(retVal)
}

### calculate batch id color list - list with names set to batch ids
calculateBatchIdColors<-function(theSortedListOfBatchIds)
{
	color <- beaRainbow(length(theSortedListOfBatchIds), v=0.7)
	names(color) <- theSortedListOfBatchIds
	return(color)
}

### calculate batch id symbol list - list with names set to batch ids
calculateBatchIdSymbols<-function(theSortedListOfBatchIds)
{
	pchValues <- getPchList(length(theSortedListOfBatchIds))
	names(pchValues) <- theSortedListOfBatchIds
	return(pchValues)
}

doSamplePcaCall<-function(theMatrixGeneData, theMinBatchSize, theBatchIds, theCompPairList, theSeed=NULL, theGeneLimit=0)
{
	logDebug("doSamplePcaCall - pre-filter matrix size=", paste(dim(theMatrixGeneData),collapse=","))
	myMatrixGeneData<-filterTheGeneDataOnBatchSize(theMatrixGeneData, theMinBatchSize, theBatchIds)
	myMatrixGeneData<-filterBasedOnGeneLimit(myMatrixGeneData, theSeed, theGeneLimit)
	logDebug("doSamplePcaCall - pre-filter matrix size=", paste(dim(theMatrixGeneData),collapse=","))
	logDebug("doSamplePcaCall - no-NA matrix length=", sum(!is.na(rowSums(myMatrixGeneData))))
	pca <- NULL
	if (sum(!is.na(rowSums(myMatrixGeneData)))<2)
	{
		logWarn("All rows contain NA values")
	}
	else
	{
		### calculate PCA
		errorsContinued <- FALSE
		##pca<-SamplePCA(theMatrixGeneData[!is.na(rowSums(theMatrixGeneData)),], usecor=FALSE, center=TRUE),
		##to fix the routine error code 1 in Lapack 9/1/11 by Nianxiang
		pca <- try( SamplePCA(myMatrixGeneData[!is.na(rowSums(myMatrixGeneData)),], usecor=FALSE, center=TRUE) )
		if(is(pca,'try-error'))
		{
			logDebug("SamplePCA first call threw a problem -- dropping outlier and retrying")
			myNoMissingMatrixGeneData<-myMatrixGeneData[!is.na(rowSums(myMatrixGeneData)),]
			myGeneVar<-apply(myNoMissingMatrixGeneData, 1, var)
			sortedGeneVar<-sort(myGeneVar, decreasing=TRUE)
			myOutlierNumber<-1
			###logDebug("pca =", pca)
			logDebug("myOutlierNumber =", myOutlierNumber)
			###logDebug("nrow(sortedGeneVar) =", length(sortedGeneVar))
			while((is(pca,'try-error'))&&(myOutlierNumber<length(sortedGeneVar)))
			{
				logDebug("SamplePCA subsequent call generated a failed (not necessarily a problem yet) ", myOutlierNumber)
				pca<-try(SamplePCA(myNoMissingMatrixGeneData[myGeneVar<sortedGeneVar[myOutlierNumber],], usecor=FALSE, center=TRUE))
				myOutlierNumber<-myOutlierNumber+1
			}
			if (myOutlierNumber>=length(sortedGeneVar))
			{
				logWarn("SamplePCA continued generating failures (data set cannot be used) ", myOutlierNumber)
				errorsContinued <- TRUE
			}
		}
		if (FALSE==errorsContinued)
		{
			stopifnotWithLogging("SamplePCA returned a null", (!is.null(pca)))
			tempCols <- colnames(pca@scores)
			tempRows <- rownames(pca@scores)
			if ((0==length(tempCols))||(0==length(tempRows)))
			{
				stopWithLogging("SamplePCA in ClassDiscovery is returning empty rows and/or columns. You need a newer version of ClassDiscovery")
			}
		}
		else
		{
			pca <- NULL
		}
	}
	if(is.null(pca))
	{
		logDebug("doSamplePcaCall - pca is null")
	}
	else if (ncol(pca@scores)<max(theCompPairList))
	{
	  logDebug("doSamplePcaCall - max(theCompPairList)=", max(theCompPairList))
	  logDebug("doSamplePcaCall - ncol(pca@scores)=", ncol(pca@scores))
	  logWarn("doSamplePcaCall - pca set to null, since number of components calculated is less than components requested")
		pca <- NULL
	}
	else
	{
	  logDebug("doSamplePcaCall - max(theCompPairList)=", max(theCompPairList))
	  logDebug("doSamplePcaCall - ncol(pca@scores)=", ncol(pca@scores))
	  logDebug("doSamplePcaCall - pca scores size=", paste(dim(pca@scores),collapse=","))
	}
	return(pca)
}

pcaErrorNormal <- function(theOutputFile, theDscResultsDir)
{
  logDebug("pcaErrorNormal - theOutputFile=", theOutputFile)
  myFile <- file(theOutputFile, "w+")
  on.exit(close(myFile))
  cat("Unable to Generate PCA Results\n", file=myFile, append=TRUE)
}

pcaErrorDSC <- function(theOutputFile, theDscResultsDir)
{
  logDebug("pcaErrorDSC - theOutputFile=", theOutputFile)
  logDebug("pcaErrorDSC - myFile=", file.path(theDscResultsDir, "error.log"))
  unlink(file.path(theDscResultsDir, "DSCOverview.tsv"))
  myFile <- file(file.path(theDscResultsDir, "error.log"), "w+")
  on.exit(close(myFile))
  cat("Unable to Generate PCA Results\n", file=myFile, append=TRUE)
}


openAndWriteIssuesLogFilePCA<-function(theOutputFile, theDscResultsDir)
{
  pcaErrorNormal(theOutputFile, theDscResultsDir)
  pcaErrorDSC(theOutputFile, theDscResultsDir)
}

doInternalPca<-function(thePca, theDscAllResults, theComponentA, theComponentB, theMatrixGeneData,
                        theDoPlainMtoMFlag, theDoCentroidsMtoMFlag, theDoDSCFlag, theDoDscPermsFileFlag,
                        theCentroidFile, thePlainFile, theDSCFile,
                        theCentrLegendAll, theCentrLegendDSC, theCentrLegendRays, theCentrLegendPoints,
                        thePlainLegendAll, thePlainLegendDSC, thePlainLegendRays, thePlainLegendPoints,
                        theCentroidTitle, thePlainTitle,
                        batchIdsForSamplesA, batchIdsForSamplesB, sortedListOfBatchIdsA, sortedListOfBatchIdsB, batchTypeNameA, batchTypeNameB,
                        theIsPcaTrendFunction, theDSCPermutations, theDSCThreads, theMinBatchSize, theSampleIds,
                        theErrorFile, theSeed, theDscResultsDir)
{
  returnValue <- NULL
  tryCatch({
    checkIfTestError()
    doInternalPcaThrowable(thePca, theDscAllResults, theComponentA, theComponentB, theMatrixGeneData,
                           theDoPlainMtoMFlag, theDoCentroidsMtoMFlag, theDoDSCFlag, theDoDscPermsFileFlag,
                           theCentroidFile, thePlainFile, theDSCFile,
                           theCentrLegendAll, theCentrLegendDSC, theCentrLegendRays, theCentrLegendPoints,
                           thePlainLegendAll, thePlainLegendDSC, thePlainLegendRays, thePlainLegendPoints,
                           theCentroidTitle, thePlainTitle,
                           batchIdsForSamplesA, batchIdsForSamplesB, sortedListOfBatchIdsA, sortedListOfBatchIdsB, batchTypeNameA, batchTypeNameB,
                           theIsPcaTrendFunction, theDSCPermutations, theDSCThreads, theMinBatchSize, theSampleIds,
                           theSeed)
  },
  error=function(e)
  {
    logWarn("2 Unable to generate PCA")
    logDebug("2 write error to ", theErrorFile)
    # unlink(dirname(theErrorFile), recursive = TRUE, force = TRUE)
    # dir.create(dirname(theErrorFile), showWarnings=FALSE, recursive=TRUE)
    openAndWriteIssuesLogFilePCA(theErrorFile, theDscResultsDir)
  })
}

doInternalPcaThrowable<-function(thePca, theDscAllResults, theComponentA, theComponentB, theMatrixGeneData,
                                 theDoPlainMtoMFlag, theDoCentroidsMtoMFlag, theDoDSCFlag, theDoDscPermsFileFlag,
                                 theCentroidFile, thePlainFile, theDSCFile,
                                 theCentrLegendAll, theCentrLegendDSC, theCentrLegendRays, theCentrLegendPoints,
                                 thePlainLegendAll, thePlainLegendDSC, thePlainLegendRays, thePlainLegendPoints,
                                 theCentroidTitle, thePlainTitle,
                                 batchIdsForSamplesA, batchIdsForSamplesB, sortedListOfBatchIdsA, sortedListOfBatchIdsB, batchTypeNameA, batchTypeNameB,
                                 theIsPcaTrendFunction, theDSCPermutations, theDSCThreads, theMinBatchSize, theSampleIds,
                                 theSeed)
{
	doTrending<-FALSE
	if(!is.null(theIsPcaTrendFunction))
	{
		doTrending<-theIsPcaTrendFunction(batchTypeNameA, unique(sortedListOfBatchIdsA))
	}
	logDebug("doInternalPca before filter - Number of rownames for pca score (sample ids) ", length(rownames(thePca@scores)))
	logDebug("doInternalPca before filter - length of batch ids for samples ", length(batchIdsForSamplesA))
	### filter by batch size
	if (theMinBatchSize>=1)
	{
		logDebug("doInternalPca - filtering batch size of ", theMinBatchSize)
		###logDebug("BEFORE ",
		###		" length(batchIdsForSamplesA)=", length(batchIdsForSamplesA),
		###		" length(theSampleIds)=", length(theSampleIds),
		###		" ncol(theMatrixGeneData)=", ncol(theMatrixGeneData),
		###		" nrow(theMatrixGeneData)=", nrow(theMatrixGeneData),
		###		" length(sortedListOfBatchIdsA)=", length(sortedListOfBatchIdsA))
		threshHoldList<-batchTreshold(batchIdsForSamplesA, theMinBatchSize)
		batchIdsForSamplesA <- batchIdsForSamplesA[which(threshHoldList)]
		batchIdsForSamplesB <- batchIdsForSamplesB[which(threshHoldList)]
		theSampleIds <- theSampleIds[which(threshHoldList)]
		theMatrixGeneData <- theMatrixGeneData[,which(threshHoldList)]
		intA <- intersect(sortedListOfBatchIdsA, batchIdsForSamplesA)
		sortedListOfBatchIdsA <- NULL
		if (!is.null(intA))
		{
			sortedListOfBatchIdsA <- sort(intA)
		}
		intB <- intersect(sortedListOfBatchIdsB, batchIdsForSamplesB)
		sortedListOfBatchIdsB <- NULL
		if (!is.null(intB))
		{
			sortedListOfBatchIdsB <- sort(intB)
		}
		###logDebug("AFTER  ",
		###		" length(batchIdsForSamplesA)=", length(batchIdsForSamplesA),
		###		" length(theSampleIds)=", length(theSampleIds),
		###		" ncol(theMatrixGeneData)=", ncol(theMatrixGeneData),
		###		" nrow(theMatrixGeneData)=", nrow(theMatrixGeneData),
		###		" length(sortedListOfBatchIdsA)=", length(sortedListOfBatchIdsA),
		###		" length(threshHoldList)=", length(threshHoldList),
		###		" sum(threshHoldList)=", sum(threshHoldList),
		###		" ncol(pca@scores)=", ncol(pca@scores),
		###		" nrow(thePca@scores)=", nrow(thePca@scores))
	}
	logDebug("doInternalPca after filter - Number of rownames for pca score (sample ids) ", length(rownames(thePca@scores)))
	logDebug("doInternalPca after filter - length of batch ids for samples ", length(batchIdsForSamplesA))
	###else
	###{
	###	logDebug("doInternalPca - NOT filtering for ", theMinBatchSize)
	###}
	extraLegend <- NULL
	diagramLabelDsc <- NULL
	diagramLabelDscXY <- NULL
	results <- NULL
	if (TRUE==theDoDSCFlag)
	{
		logDebug("doInternalPca - pvalueDSC")
		### TODO: make calls to doInternalPca return the list of DSC objects created here
		results <- pvalueDSC(thePca, batchIdsForSamplesA, theDSCPermutations, theComponentA, theComponentB, theDSCThreads, theSeed)
		if(TRUE==theDoDscPermsFileFlag)
		{
			openAndWriteDscPermsFile(results, dirname(theDSCFile), "ANY", theComponentA, theComponentB)
		}
		diagramLabelDsc <- theDscAllResults@mDSC
		diagramLabelDscXY <- results@mDSC
		DSCLabels <- c(
				"First PCA Component:",
				"First Component FVE (%):",
				"Second PCA Component:",
				"Second Component FVE (%):",
				"",
				paste("Disp. Sep. Crit. (DSC) (", theComponentA, ",", theComponentB, "):", sep=""),
				paste("Disp. within groups (Dw) (", theComponentA, ",", theComponentB, "):", sep=""),
				paste("Disp. between groups (Db) (", theComponentA, ",", theComponentB, "):", sep=""),
				paste("DSC pvalue(", theComponentA, ",", theComponentB, "):", sep=""),
				"",
				"Disp. Sep. Crit. (DSC):",
				"Disp. within groups (Dw):",
				"Disp. between groups (Db):",
				"DSC pvalue:")
		DSCvalues <- c(
				theComponentA,
				calculateFVEForSingleComponent(thePca, theComponentA)*100,
				theComponentB,
				calculateFVEForSingleComponent(thePca, theComponentB)*100,
				"",
				results@mDSC,
				results@mDW,
				results@mDB,
				results@mPvalue,
				"",
				theDscAllResults@mDSC,
				theDscAllResults@mDW,
				theDscAllResults@mDB,
				theDscAllResults@mPvalue)
		### combine DSC and DSC pvalue data
		extraLegend <- combineStringPairs(DSCLabels, DSCvalues, " ", FALSE, FALSE, FALSE)
		combineForFile <- combineStringPairs(DSCLabels, DSCvalues, "\t", FALSE, FALSE, FALSE)
		### save value objects for later reuse in making DSC summary / passing in PCA-PVALUE-DSC / file is *__CompDSC.RData
		saveCompDscData(paste(theDSCFile, "__CompDSC.RData", sep=""), results, theComponentA, theComponentB)
	}
	logDebug("write PCA file direct 1 ", dirname(dirname(theCentroidFile)))
	writePCAAnnotations(dirname(dirname(theCentroidFile)), theDSCPermutations, thePca, NULL, results, theComponentA, theComponentB )
	### calculate batch id color list - list with names set to batch ids
	batchIdColorsA <- calculateBatchIdColors(sortedListOfBatchIdsA)
	batchIdColorsB <- NULL
	if (!is.null(sortedListOfBatchIdsB))
	{
		batchIdColorsB <- calculateBatchIdColors(sortedListOfBatchIdsB)
	}
	### calculate batch id symbol list - list with names set to batch ids
	batchIdSymbols <- calculateBatchIdSymbols(sortedListOfBatchIdsA)
	if (!is.null(sortedListOfBatchIdsB))
	{
		batchIdSymbols <- calculateBatchIdSymbols(sortedListOfBatchIdsB)
	}
	### create output
	if(TRUE==theDoCentroidsMtoMFlag)
	{
		openAndWritePcaAnalysisImage(thePca, theComponentA, theComponentB, theCentroidFile, theCentroidTitle,
				theCentrLegendAll, theCentrLegendDSC, theCentrLegendRays, theCentrLegendPoints,
				batchIdsForSamplesA, batchIdsForSamplesB,
				sortedListOfBatchIdsA, sortedListOfBatchIdsB,
				batchTypeNameA, batchTypeNameB,
				batchIdColorsA, batchIdColorsB,
				batchIdSymbols, TRUE, doTrending, extraLegend, diagramLabelDsc, diagramLabelDscXY )
	}
	if(TRUE==theDoPlainMtoMFlag)
	{
		openAndWritePcaAnalysisImage(thePca, theComponentA, theComponentB, thePlainFile, thePlainTitle,
				thePlainLegendAll, thePlainLegendDSC, thePlainLegendRays, thePlainLegendPoints,
				batchIdsForSamplesA, batchIdsForSamplesB,
				sortedListOfBatchIdsA, sortedListOfBatchIdsB,
				batchTypeNameA, batchTypeNameB,
				batchIdColorsA, batchIdColorsB,
				batchIdSymbols, FALSE, doTrending, extraLegend, diagramLabelDsc, diagramLabelDscXY )
	}
	if(TRUE==theDoDSCFlag)
	{
		openAndWriteDSCFile(combineForFile, theDSCFile)
	}
}

calculateFVE<-function(thePca, theComponentA, theComponentB)
{
	return(calculateFVEForSingleComponent(thePca, theComponentA)+calculateFVEForSingleComponent(thePca, theComponentB))
}

calculateFVEForSingleComponent<-function(thePca, theComponent)
{
	varianceValue <- (thePca@variances[[theComponent]])/sum(thePca@variances)
	return(varianceValue)
}

writePcaAnalysis_Image<-function(thePca, theComponentA, theComponentB, theTitle,
		theBatchIdsForSamplesRays, theBatchIdsForSamplesPoints,
		theSortedListOfBatchIdsRays, theSortedListOfBatchIdsPoints,
		theBatchTypeRays, theBatchTypePoints,
		theBatchIdColorsRays, theBatchIdColorsPoints,
		theBatchIdSymbols, theDoCentroids, theDoTrending, theDiagramLabelDsc, theDiagramLabelDscXY,
		theCex=0.5, theTitleCex=1.0, theSubCex=1.0)
{
	center1<-NULL
	center2<-NULL
	if (TRUE==theDoCentroids)
	{
		###logDebug("c1")
		center1<-tapply(thePca@scores[,theComponentA], factor(theBatchIdsForSamplesRays), mean)
		###logDebug("c2")
		center2<-tapply(thePca@scores[,theComponentB], factor(theBatchIdsForSamplesRays), mean)
	}
	xlim <- range(thePca@scores[,theComponentA], finite=TRUE )
	ylim <- range(thePca@scores[,theComponentB], finite=TRUE )
	###logDebug("xlim[", theComponentA, "]=", paste(xlim, collapse=", "))
	###logDebug("ylim[", theComponentB, "]=", paste(ylim, collapse=", "))
	###logDebug("max=", max(thePca@scores))
	###logDebug("min=", min(thePca@scores))
	###logDebug("max[", theComponentA, "]=", max(thePca@scores[,theComponentA]))
	###logDebug("min[", theComponentA, "]=", min(thePca@scores[,theComponentA]))
	###logDebug("max[", theComponentB, "]=", max(thePca@scores[,theComponentB]))
	###logDebug("min[", theComponentB, "]=", min(thePca@scores[,theComponentB]))
	###logDebug("for batch")
	###logDebug("writePcaAnalysis_Image theSortedListOfBatchIdsRays ", paste(theSortedListOfBatchIdsRays, collapse=", "))
	###logDebug("writePcaAnalysis_Image theBatchIdsForSamplesRays ", paste(theBatchIdsForSamplesRays, collapse=", "))
	###logDebug("writePcaAnalysis_Image theBatchIdColorsRays ", paste(theBatchIdColorsRays, collapse=", "))
	###logDebug("writePcaAnalysis_Image theBatchTypeRays ", paste(theBatchTypeRays, collapse=", "))
	previousCentroidXY<-NULL
	if (is.null(theBatchIdsForSamplesPoints))
	{
		###logDebug("copying rays to points")
		theSortedListOfBatchIdsPoints <- theSortedListOfBatchIdsRays
		theBatchIdsForSamplesPoints <- theBatchIdsForSamplesRays
		theBatchIdColorsPoints <- theBatchIdColorsRays
		theBatchTypePoints <- theBatchTypeRays
	}
	#######################################################
	### do points
	#######################################################
	firstRun <- TRUE
	stopifnotWithLogging("Number of rownames for pca score (sample ids) should match length of batch ids for samples", (length(rownames(thePca@scores))==length(theBatchIdsForSamplesPoints)))
	logDebug("writePcaAnalysis_Image - Number of rownames for pca score (sample ids) ", length(rownames(thePca@scores)))
	logDebug("writePcaAnalysis_Image - length of batch ids for samples ", length(theBatchIdsForSamplesPoints))
	plotOrderForBatchIds <- sortBatchesBasedOnSize(theSortedListOfBatchIdsPoints, theBatchIdsForSamplesPoints)
	for (batchId in plotOrderForBatchIds)
	{
		###logDebug("batchId ", batchId)
		color <- theBatchIdColorsPoints[[batchId]]
		###logDebug("writePcaAnalysis_Image - Do points, get batch id symbol for ", batchId, " from ", paste(names(theBatchIdSymbols), collapse=", "))
		symbolPchValue <- theBatchIdSymbols[[batchId]]
		###logDebug("BatchId=", batchId, " color=", color, " symbol=", symbolPchValue)
		groupi<-thePca@scores[factor(theBatchIdsForSamplesPoints)==batchId, c(theComponentA,theComponentB)]
		if (FALSE==is.array(groupi))
		{
			### we have a one dimensional array, recast to get correct internal format
			groupi <- array(groupi,dim=c(1,2))
		}
		if (TRUE==firstRun)
		{
			firstRun <- FALSE
			### create plot for first occurance
			###logDebug("plot")
			varianceValue <- calculateFVE(thePca, theComponentA, theComponentB)
			subTitle <- paste("FVE(", theComponentA, ",", theComponentB, "): ", round(varianceValue*100,3) ,"%", sep="")
			if ((!is.null(theDiagramLabelDscXY))&&(!is.null(theDiagramLabelDsc)))
			{
				subTitle <- paste("FVE(", theComponentA, ",", theComponentB, "): ", round(varianceValue*100,3) ,"%",
													"   DSC(", theComponentA, ",", theComponentB, "): ", round(theDiagramLabelDscXY,3),
													"   DSC: ", round(theDiagramLabelDsc,3),
													sep="")
			}
			plot(groupi[,1], groupi[,2], cex=1.5, xlab=paste("Comp. ", theComponentA, sep=""), ylab=paste("Comp. ", theComponentB, sep=""),
					col=color, xlim=xlim, ylim=ylim, xaxs="r", yaxs="r",
					main=theTitle, cex.main=theTitleCex, cex.sub=theSubCex, pch=symbolPchValue,
					sub=subTitle
			)
		}
		else
		{
			### plot just the points
			###logDebug("points 1")
			points(groupi[,1], groupi[,2], cex=1.5, col=color, pch=symbolPchValue)
		}
	}
	#######################################################
	### do lines
	#######################################################
	for (batchId in theSortedListOfBatchIdsRays)
	{
		color <- theBatchIdColorsRays[[batchId]]
		groupi<-thePca@scores[factor(theBatchIdsForSamplesRays)==batchId, c(theComponentA,theComponentB)]
		if (FALSE==is.array(groupi))
		{
			### we have a one dimensional array, recast to get correct internal format
			groupi <- array(groupi,dim=c(1,2))
		}
		if (TRUE==theDoCentroids)
		{
			for (j in (1:nrow(groupi)))
			{
				###logDebug("segments")
				segments( groupi[j,1], groupi[j,2], center1[[batchId]], center2[[batchId]], col=color , lwd=0.5)
			}
		}
	}
	#######################################################
	### do centroids
	#######################################################
	if (TRUE==theDoCentroids)
	{
		for (batchId in theSortedListOfBatchIdsRays)
		{
			color <- theBatchIdColorsRays[[batchId]]
			points(center1[[batchId]], center2[[batchId]], lwd=1.5, col=color, bg=color, pch=21, cex=1.5)
			previousCentroidXY<-c(previousCentroidXY, center1[[batchId]], center2[[batchId]])
			if(TRUE==all.equal(theBatchIdsForSamplesPoints, theBatchIdsForSamplesRays))
			{
				symbolPchValue <- theBatchIdSymbols[[batchId]]
				points(center1[[batchId]], center2[[batchId]], lwd=1.5, col="black", pch=symbolPchValue, cex=1.5)
			}
			else
			{
				points(center1[[batchId]], center2[[batchId]], lwd=1.5, col="black", pch=21, cex=1.5)
			}
		}
	}
	###logDebug("abline")
	abline(h=0, v=0, col="brown", lty=2)
	###logDebug("theDoTrending")
	###logDebug("theDoTrending -- ", theDoTrending)
	if (TRUE==theDoTrending)
	{
		if (length(previousCentroidXY)>2)
		{
			for(i in seq(from=1, to=(length(previousCentroidXY)-1), by=2 ))
			{
				if ((i+3)<=length(previousCentroidXY))
				{
					###logDebug("Centroid points (", previousCentroidXY[[i]], ") , (", previousCentroidXY[[i+1]], ") , (", previousCentroidXY[[i+2]], ") , (", previousCentroidXY[[i+3]], ")")
					###logDebug("Centroid distances =", )
					distPoints <- sqrt((previousCentroidXY[[i+2]]-previousCentroidXY[[i]])^2+(previousCentroidXY[[i+3]]-previousCentroidXY[[i+1]])^2)
					if (!is.na(distPoints))
					{
						if (distPoints>0.000001)
						{
							arrows( previousCentroidXY[[i]], previousCentroidXY[[i+1]], previousCentroidXY[[i+2]], previousCentroidXY[[i+3]], col="black" , lwd=2.0)
						}
					}
				}
			}
		}
	}
}

openAndWritePcaAnalysisImage<-function(thePca, theComponentA, theComponentB, thePcaFile, theTitle,
		theLegendAll, theLegendDSC, theLegendRays, theLegendPoints,
		theBatchIdsForSamplesRays, theBatchIdsForSamplesPoints,
		theSortedListOfBatchIdsRays, theSortedListOfBatchIdsPoints,
		theBatchTypeRays, theBatchTypePoints,
		theBatchIdColorsRays, theBatchIdColorsPoints,
		theBatchIdSymbols, theDoCentroids, theDoTrending,
		theExtraLegend, theDiagramLabelDsc, theDiagramLabelDscXY)
{
		###logDebug("openAndWritePcaAnalysisImage -- theSortedListOfBatchIdsRays = ", paste(theSortedListOfBatchIdsRays, collapse=","))
	###logDebug("openAndWritePcaAnalysisImage -- theSortedListOfBatchIdsPoints = ", paste(theSortedListOfBatchIdsPoints, collapse=","))
	openAndWritePcaAnalysisImage_Diagram(thePcaFile, thePca, theComponentA, theComponentB, theTitle,
			theBatchIdsForSamplesRays, theBatchIdsForSamplesPoints,
			theSortedListOfBatchIdsRays, theSortedListOfBatchIdsPoints,
			theBatchTypeRays, theBatchTypePoints,
			theBatchIdColorsRays, theBatchIdColorsPoints,
			theBatchIdSymbols, theDoCentroids, theDoTrending, theDiagramLabelDsc, theDiagramLabelDscXY)
	plotSize <- nrow(thePca@scores)
	listOfFiles <- NULL
	if (!is.null(theBatchIdsForSamplesPoints))
	{
		listOfFiles <- c(listOfFiles, openAndWritePcaAnalysisImage_LegendRays(theLegendRays, theBatchIdsForSamplesRays, theSortedListOfBatchIdsRays, theBatchTypeRays,
			theBatchIdColorsRays, theBatchIdSymbols, plotSize))
		listOfFiles <- c(listOfFiles, openAndWritePcaAnalysisImage_LegendPoints(theLegendPoints, theBatchIdsForSamplesPoints, theSortedListOfBatchIdsPoints, theBatchTypePoints,
				theBatchIdColorsPoints, theBatchIdSymbols, plotSize))
	}
	else
	{
		listOfFiles <- c(listOfFiles, openAndWritePcaAnalysisImage_LegendPoints(theLegendPoints, theBatchIdsForSamplesRays, theSortedListOfBatchIdsRays, theBatchTypeRays,
				theBatchIdColorsRays, theBatchIdSymbols, plotSize))
	}
	listOfFiles <- c(openAndWritePcaAnalysisImage_LegendDSC(theLegendDSC, theExtraLegend), listOfFiles)
	#openAndWritePcaAnalysisImage_LegendAll(theLegendAll, theTitle, theBatchIdsForSamplesRays, theSortedListOfBatchIdsRays, theBatchTypeRays, theBatchIdColorsRays,
	#																			 theBatchIdsForSamplesPoints, theSortedListOfBatchIdsPoints, theBatchTypePoints, theBatchIdColorsPoints,
	#																			 theBatchIdSymbols, theExtraLegend, plotSize)
	mbatchStandardCombineLegends(theTitle, theLegendAll, listOfFiles)
}

openAndWritePcaAnalysisImage_LegendRays<-function(theLegendRays, theBatchIdsForSamplesRays,
																									theSortedListOfBatchIdsRays, theBatchTypeRays,
																									theBatchIdColorsRays, theBatchIdSymbols, theSize)
{
	sortedBatch <- getListOfBatchNamesWithCounts(theBatchIdsForSamplesRays, theSortedListOfBatchIdsRays)
	mbatchStandardLegend(paste("Rays: ", theBatchTypeRays, " (", theSize, ")", sep=""),
											 paste("MBatch", packageDescription("MBatch")["Version"], sep=" "),
											 sortedBatch, as.vector(theBatchIdColorsRays), NULL, theLegendRays)
#	Cairo PNG(filename=theLegendRays, width = 600, height = 1000, pointsize=24, bg = "transparent")
#	on.exit(dev.off(), add = TRUE)
#	plot.new()
#	if (!is.null(theBatchIdSymbols))
#	{
#		writePcaAnalysis_Legend("center", theBatchIdsForSamplesRays, theSortedListOfBatchIdsRays, paste("Rays: ", theBatchTypeRays, " (", theSize, ")", sep=""), theBatchIdColorsRays, theBatchIdSymbols, 0.5)
#	}
#	else
#	{
#		writePcaAnalysis_Legend("center", theBatchIdsForSamplesRays, theSortedListOfBatchIdsRays, paste("Rays: ", theBatchTypeRays, " (", theSize, ")", sep=""), theBatchIdColorsRays, NULL, 0.5)
#	}
	return(theLegendRays)
}

openAndWritePcaAnalysisImage_LegendPoints<-function(theLegendPoints, theBatchIdsForSamplesPoints,
																										theSortedListOfBatchIdsPoints, theBatchTypePoints,
																										theBatchIdColorsPoints, theBatchIdSymbols, theSize)
{
	sortedBatch <- getListOfBatchNamesWithCounts(theBatchIdsForSamplesPoints, theSortedListOfBatchIdsPoints)
	mbatchStandardLegend(paste("Points: ", theBatchTypePoints, " (", theSize, ")", sep=""),
											 paste("MBatch", packageDescription("MBatch")["Version"], sep=" "),
											 sortedBatch, as.vector(theBatchIdColorsPoints), as.vector(theBatchIdSymbols), theLegendPoints)
#	if (!is.null(theBatchIdsForSamplesPoints))
#	{
#		Cairo PNG(filename=theLegendPoints, width = 600, height = 1000, pointsize=24, bg = "transparent")
#		on.exit(dev.off(), add = TRUE)
#		plot.new()
#		writePcaAnalysis_Legend("center", theBatchIdsForSamplesPoints, theSortedListOfBatchIdsPoints, paste("Points: ", theBatchTypePoints, " (", theSize, ")", sep=""), theBatchIdColorsPoints, theBatchIdSymbols)
#	}
	return(theLegendPoints)
}

openAndWritePcaAnalysisImage_LegendDSC<-function(theLegendDSC, theExtraLegend)
{
	mbatchStandardLegend("Dispersion Metrics",
											 paste("MBatch", packageDescription("MBatch")["Version"], sep=" "),
											 theExtraLegend, NULL, NULL, theLegendDSC)
#	if (!is.null(theExtraLegend))
#	{
#		Cairo PNG(filename=theLegendDSC, width = 600, height = 1000, pointsize=24, bg = "transparent")
#		on.exit(dev.off(), add = TRUE)
#		plot.new()
#		leg end("center", legend=theExtraLegend, title="Dispersion Metrics", cex=0.70)
#	}
	return(theLegendDSC)
}

openAndWritePcaAnalysisImage_Diagram<-function(thePcaFile, thePca, theComponentA, theComponentB, theTitle,
		theBatchIdsForSamplesRays, theBatchIdsForSamplesPoints,
		theSortedListOfBatchIdsRays, theSortedListOfBatchIdsPoints,
		theBatchTypeRays, theBatchTypePoints,
		theBatchIdColorsRays, theBatchIdColorsPoints,
		theBatchIdSymbols, theDoCentroids, theDoTrending, theDiagramLabelDsc, theDiagramLabelDscXY)
{
	###logDebug("Write", thePcaFile)
	CairoPNG(filename=thePcaFile, width = 1000, height = 1000, pointsize=24, bg = "transparent")
	on.exit(dev.off(), add = TRUE)
	writePcaAnalysis_Image(thePca, theComponentA, theComponentB, theTitle,
			theBatchIdsForSamplesRays, theBatchIdsForSamplesPoints,
			theSortedListOfBatchIdsRays, theSortedListOfBatchIdsPoints,
			theBatchTypeRays, theBatchTypePoints,
			theBatchIdColorsRays, theBatchIdColorsPoints,
			theBatchIdSymbols, theDoCentroids, theDoTrending, theDiagramLabelDsc, theDiagramLabelDscXY)
	mtext(paste("MBatch", packageDescription("MBatch")["Version"], sep=" "), side=1, adj=1, outer=FALSE, line=4, cex=0.50)
}

####################################################################
###
####################################################################


openAndWriteDSCFile <- function(theExtraLegend, theDSCFile)
{
  #logDebug("openAndWriteDSCFile -- start")
  #logDebug(paste(theDSCFile, sep=" + ", collapse=" | "))
  #logDebug(theExtraLegend)
  if (!file.exists(theDSCFile))
	{
		myFile <- file(theDSCFile, "w+")
		on.exit(close(myFile))
		for (line in theExtraLegend)
		{
			cat(line, file=myFile)
			cat("\n", file=myFile)
		}
	}
	#logDebug("openAndWriteDSCFile -- end")
}

openAndWriteFveFile <- function(thePca, theOutputDir)
{
	###logDebug("openAndWriteFveFile -- start")
	fveFile <- makePcaFileName_TXT(theOutputDir, "ANY", "NA", "NA", "FVE")
	if (!file.exists(fveFile))
	{
		myFile <- file(fveFile, "w+")
		on.exit(close(myFile))
		cat("Component Number\tFVE (%)\tCumulative FVE (%)\n", file=myFile, append=TRUE)
		runningTotal <- 0
		for (myIndex in 1:length(thePca@variances))
		{
			fve <- calculateFVEForSingleComponent(thePca, myIndex)
			runningTotal <- runningTotal + fve
			cat(myIndex, file=myFile, append=TRUE)
			cat("\t", file=myFile, append=TRUE)
			cat(round(fve*100, 3), file=myFile, append=TRUE)
			cat("\t", file=myFile, append=TRUE)
			cat(round(runningTotal*100, 3), file=myFile, append=TRUE)
			cat("\n", file=myFile, append=TRUE)
		}
	}
	###logDebug("openAndWriteFveFile -- end")
}

openAndWriteWeightFile <- function(thePca, theGeneList, theOutputDir)
{
	### TODO: limit this to the components from the output
	### components <- order(unique(theCompPairList))
	### myDataFrame <- thePca@scores[,components]
	###logDebug("openAndWriteWeightFile -- start")
	weightFile <- makePcaFileName_TXT(theOutputDir, "ANY", "NA", "NA", "Weights")
	if (!file.exists(weightFile))
	{
		### thePca@components
		myFile <- file(weightFile, "w+")
		on.exit(close(myFile))
		for (myIndex in 1:ncol(thePca@components))
		{
			cat(paste("\tPCA ", myIndex, sep=""), file=myFile, append=TRUE)
		}
		cat("\n", file=myFile, append=TRUE)
		for (myRow in 1:nrow(thePca@components))
		{
			myGene <- theGeneList[[myRow]]
			myCol <- thePca@components[myRow,]
			cat(paste(c(myGene, myCol), sep="\t", collapse="\t"),file=myFile, append=TRUE)
			cat("\n", file=myFile, append=TRUE)
		}
	}
	###logDebug("openAndWriteWeightFile -- end")
}

openAndWriteDscAllFile <- function(thePca, thePvalueDscObj, theOutputDir, theNameOverrideBatch=NULL, theNameOverrideCompA=NULL, theNameOverrideCompB=NULL)
{
	###logDebug("openAndWriteDscAllFile -- start")
	dscAllFile <- makePcaFileName_TXT(theOutputDir, "ALL", "NA", "NA", "DSC")
	if (!is.null(theNameOverrideBatch))
	{
		dscAllFile <- makePcaFileName_TXT(theOutputDir, theNameOverrideBatch, theNameOverrideCompA, theNameOverrideCompB, "DSC")
		### save value objects for later reuse in making DSC summary / passing in PCA-PVALUE-DSC / file is *__OverallDSC.RData
		saveOverallDscData(paste(dscAllFile, "__OverallDSC.RData", sep=""), thePvalueDscObj)
	}
	if (!file.exists(dscAllFile))
	{
		myFile <- file(dscAllFile, "w+")
		on.exit(close(myFile))
		cat("Component Number\tFVE\tDSC\tDB\tDW\tpValue\n", file=myFile, append=TRUE)
		for (myIndex in 1:length(thePvalueDscObj@mListOfGeneDSC))
		{
			cat(myIndex, file=myFile, append=TRUE)
			cat("\t", file=myFile, append=TRUE)
			cat(round(calculateFVEForSingleComponent(thePca, myIndex)*100, 8), file=myFile, append=TRUE)
			cat("\t", file=myFile, append=TRUE)
			cat(round(thePvalueDscObj@mListOfGeneDSC[[myIndex]], 8), file=myFile, append=TRUE)
			cat("\t", file=myFile, append=TRUE)
			cat(round(thePvalueDscObj@mListOfGeneDB[[myIndex]], 8), file=myFile, append=TRUE)
			cat("\t", file=myFile, append=TRUE)
			cat(round(thePvalueDscObj@mListOfGeneDW[[myIndex]], 8), file=myFile, append=TRUE)
			cat("\t", file=myFile, append=TRUE)
			cat(round(thePvalueDscObj@mListOfGenePvalue[[myIndex]], 8), file=myFile, append=TRUE)
			cat("\n", file=myFile, append=TRUE)
		}
	}
}

openAndWriteDscPermsFile <- function(thePvalueDscObj, theOutputDir, theNameOverrideBatch=NULL, theNameOverrideCompA=NULL, theNameOverrideCompB=NULL)
{
  logDebug("openAndWriteDscPermsFile -- start")
	dscPermsFile <- makePcaFileName_TXT(theOutputDir, "ALL", "NA", "NA", "DSCPerms")
	if (!is.null(theNameOverrideBatch))
	{
	  logDebug("openAndWriteDscPermsFile -- redo name")
	  dscPermsFile <- makePcaFileName_TXT(theOutputDir, theNameOverrideBatch, theNameOverrideCompA, theNameOverrideCompB, "DSCPerms")
	}
	logDebug("openAndWriteDscPermsFile -- check file")
	if (!file.exists(dscPermsFile))
	{
	  logDebug("openAndWriteDscPermsFile -- no file, write")
		counter <- 1
		if (length(thePvalueDscObj@mListOfResults)>1)
		{
		  myFile <- file(dscPermsFile, "w+")
		  on.exit(close(myFile))
		  cat("perm", file=myFile, append=TRUE)
			cat("\t", file=myFile, append=TRUE)
			cat("DSC", file=myFile, append=TRUE)
			cat("\t", file=myFile, append=TRUE)
			cat("Db", file=myFile, append=TRUE)
			cat("\t", file=myFile, append=TRUE)
			cat("Dw", file=myFile, append=TRUE)
			for(x in 1:10)
			{
				cat("\t", file=myFile, append=TRUE)
				cat("DSC[", file=myFile, append=TRUE)
				cat(x, file=myFile, append=TRUE)
				cat("]", file=myFile, append=TRUE)
			}
			for(x in 1:10)
			{
				cat("\t", file=myFile, append=TRUE)
				cat("Db[", file=myFile, append=TRUE)
				cat(x, file=myFile, append=TRUE)
				cat("]", file=myFile, append=TRUE)
			}
			for(x in 1:10)
			{
				cat("\t", file=myFile, append=TRUE)
				cat("Dw[", file=myFile, append=TRUE)
				cat(x, file=myFile, append=TRUE)
				cat("]", file=myFile, append=TRUE)
			}
			cat("\n", file=myFile, append=TRUE)
			for (permResult in thePvalueDscObj@mListOfResults)
			{
				######
				cat(counter, file=myFile, append=TRUE)
				cat("\t", file=myFile, append=TRUE)
				cat(round(permResult@mDSC, 8), file=myFile, append=TRUE)
				######
				cat("\t", file=myFile, append=TRUE)
				cat(round(permResult@mDB, 8), file=myFile, append=TRUE)
				######
				cat("\t", file=myFile, append=TRUE)
				cat(round(permResult@mDW, 8), file=myFile, append=TRUE)
				######
				for(value in permResult@mListOfDSCbyGene[1:10])
				{
					cat("\t", file=myFile, append=TRUE)
					cat(round(value, 8), file=myFile, append=TRUE)
				}
				######
				for(value in permResult@mListOfDBbyGene[1:10])
				{
					cat("\t", file=myFile, append=TRUE)
					cat(round(value, 8), file=myFile, append=TRUE)
				}
				######
				for(value in permResult@mListOfDWbyGene[1:10])
				{
					cat("\t", file=myFile, append=TRUE)
					cat(round(value, 8), file=myFile, append=TRUE)
				}
				######
				cat("\n", file=myFile, append=TRUE)
				counter <- counter + 1
			}
		}
		else
		{
		  logDebug("openAndWriteDscPermsFile -- nothing to write")
		}
	}
	else
	{
	  logDebug("openAndWriteDscPermsFile -- file exists")
	}
}

openAndWriteScoresFile <- function(theOutputDir, thePca, theGeneIds)
{
	###logDebug("openAndWriteScoresFile -- start")
	sampleLocatorFile <- makePcaFileName_TXT(theOutputDir, "ANY", "NA", "NA", "Scores")
	if (!file.exists(sampleLocatorFile))
	{
		myFile <- file(sampleLocatorFile, "w+")
		on.exit(close(myFile))
		cat(paste(c("Component Number", theGeneIds), collapse="\t"), file=myFile, append=TRUE)
		cat("\n", file=myFile, append=TRUE)
		for (myIndex in 1:ncol(thePca@scores))
		{
			pcaId <- c(paste("PCA ", myIndex, sep=""))
			scoreList <- paste(thePca@scores[,myIndex], sep="\t", collapse="\t")
			dataList <- paste(pcaId, scoreList, sep="\t", collapse="\t")
			cat(dataList, file=myFile, append=TRUE)
			cat("\n", file=myFile, append=TRUE)
		}
	}
	###logDebug("openAndWriteScoresFile -- end")
}

writeSharedFveWeightScoresFiles <- function(theDoSampleLocatorFlag, theOutputDir, thePca, theSampleIds, theGeneIds)
{
	###logDebug("writeSharedDscFveWeightScoresFiles -- start")
	openAndWriteFveFile(thePca, theOutputDir)
	openAndWriteWeightFile(thePca, theSampleIds, theOutputDir)
	if (TRUE==theDoSampleLocatorFlag)
	{
		openAndWriteScoresFile(theOutputDir, thePca, theGeneIds)
	}
	###logDebug("writeSharedDscFveWeightScoresFiles -- end")
}

####################################################################
###
####################################################################

makePcaFileName_PNG<-function(theDir, theAllOrBatch, theComponentA, theComponentB, theDiagramOrLegend, theLegendType="")
{
	###stopifnotWithLogging("theCentroidOrPlain should be PCA-Plus or PcaPlain", (("PCA-Plus"==theCentroidOrPlain)||("PCAPlain"==theCentroidOrPlain)))
	###stopifnotWithLogging("thePlotType should be Many2Many|One2Many|DualBatch", (("Many2Many"==thePlotType)||("One2Many"==thePlotType)||("DualBatch"==thePlotType)))
	stopifnotWithLogging("theDiagramOrLegend should be Diagram|Legend|Complete", (("Diagram"==theDiagramOrLegend)||("Legend"==theDiagramOrLegend)))
	if ("Legend"!=theDiagramOrLegend)
	{
		if (""!=theLegendType)
		{
			stopWithLogging(paste("theLegendType is ", theLegendType, " which must be empty for theDiagramOrLegend=", theDiagramOrLegend, sep=""))
		}
	}
	componentSection <- ""
	if (("NA"!=theComponentA)&&("NA"!=theComponentB))
	{
		componentSection <- paste("_Comp", theComponentA, "_Comp", theComponentB, sep="")
	}
	if ("Legend"==theDiagramOrLegend)
	{
		createDirPlusFilename(theDir, theAllOrBatch, componentSection, "_", theDiagramOrLegend, "-", theLegendType, ".png")
	}
	else
	{
		createDirPlusFilename(theDir, theAllOrBatch, componentSection, "_", theDiagramOrLegend, ".png")
	}
}

makePcaFileName_TXT<-function(theDir, theAllOrBatch, theComponentA, theComponentB, theDataOrError)
{
	stopifnotWithLogging("theDataOrError should be DSC|DSCPerms|Scores|FVE|ERROR",
			(("DSC"==theDataOrError)||("DSCPerms"==theDataOrError)||("Scores"==theDataOrError)||("FVE"==theDataOrError)||("ERROR"==theDataOrError)||("Weights"==theDataOrError)))
	componentSection <- ""
	if (("NA"!=theComponentA)&&("NA"!=theComponentB))
	{
		componentSection <- paste("_Comp", theComponentA, "_Comp", theComponentB, sep="")
	}
	createDirPlusFilename(theDir, theAllOrBatch, componentSection, "_", theDataOrError, ".txt")
}

####################################################################
###
####################################################################

writePcaDataFilesForDataset<-function(theOutputDir, theDSCPermutations, thePca, theDscAllResults, theCompPairList, theTitle)
{
	logDebug("write writePcaDataFilesForDataset")
  checkDirForCreation(theOutputDir)
	writePCAValuesTSV(theOutputDir, thePca, theCompPairList, theTitle)
	writePCAAnnotations(theOutputDir, theDSCPermutations, thePca, theDscAllResults, NULL, NULL, NULL )
}

convertDataframeToMatrix<-function(theDataFrame, theHeader)
{
	myVector <- c(theHeader, colnames(theDataFrame))
	for (rowCounter in 1:(nrow(theDataFrame)))
	{
		myVector <- c(myVector, rownames(theDataFrame)[rowCounter], theDataFrame[rowCounter,] )
	}
	foo <- matrixWithIssues(data=myVector, nrow=nrow(theDataFrame)+1, ncol=ncol(theDataFrame)+1, byrow=TRUE)
	return(foo)
}

convertDataframeToMatrixAddingNonHeaderRow<-function(theDataFrame, theHeader, theNonHeader)
{
	myVector <- c(theHeader, colnames(theDataFrame))
	myVector <- c(myVector, theNonHeader)
	for (rowCounter in 1:(nrow(theDataFrame)))
	{
		myVector <- c(myVector, rownames(theDataFrame)[rowCounter], theDataFrame[rowCounter,] )
	}
	foo <- matrixWithIssues(data=myVector, nrow=nrow(theDataFrame)+2, ncol=ncol(theDataFrame)+1, byrow=TRUE)
	return(foo)
}

writePCAValuesTSV <- function(theOutputDir, thePca, theCompPairList, theTitle)
{
	logDebug("writePCAValuesTSV")
	### File for PCA Fractional Variances for components in PCA plot.
	### PCAFractionalVariances.tsv is a tab-delimited file containing two rows.
	### The top row is PC1	PC2	PC3 PC$ and so on. The second row is the fractional variance value.
	### PC1	PC2	PC3	PC4
  logDebug("writePCAValuesTSV title file")
  titleFile <- cleanFilePath(theOutputDir, "PCA_Title.txt")
  writeTitleFile(theTitle, titleFile)
  logDebug("writePCAValuesTSV outputFile")
  outputFile <- cleanFilePath(theOutputDir, "PCAValues.tsv")
  logDebug("writePCAValuesTSV components")
	components <- order(unique(theCompPairList))
	logDebug("writePCAValuesTSV myDataFrame")
	myDataFrame <- thePca@scores[,components]
	# values should go into first row and have row label of FVE
	logDebug("writePCAValuesTSV - values")
	values <- c("FVE", lapply(components, function(comp) {
		calculateFVEForSingleComponent(thePca, comp)
	}))
	logDebug("writePCAValuesTSV - matrix")
	myMatrix <- convertDataframeToMatrixAddingNonHeaderRow(myDataFrame, "Id", values)
	logDebug("writePCAValuesTSV - write")
	write.table(myMatrix, file=outputFile, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	logDebug("writePCAValuesTSV finished ", outputFile)
}

writePCAAnnotations <- function(theOutputDir, thePerms, thePca, theDscAllResults, theDscResults, theComponentA, theComponentB )
{
	logDebug("writePCAAnnotations")
	outputFile <- cleanFilePath(theOutputDir, "PCAAnnotations.tsv")
	##############################################
	### File for annotations related to a particular diagram. (DSC values for the diagram, that is, overall DSC values can go here.)
	### PCADiagramAnnotations.tsv is a tab-delimited file containing two columns, Annotation and Value.
	### Annotation is a display string for a value associated with a diagram and Value is the value to go with the display string.
	### Annotation	Value
	if ((!is.null(theDscAllResults))&&(!file.exists(outputFile)))
	{
		logDebug("writePCAAnnotations - header and all data")
		version <- getMBatchVersion()
		#packageDescription("MBatch")["Version"]
		if (!is.null(installed.packages()[ ,"Version"]["MBatch"]))
		{
			version <- as.character(unlist(installed.packages()[ ,"Version"]["MBatch"]))
		}
		myDataFrame <- data.frame(
			Type="Run",
			SubType="-",
			Annotation="MBatch Version",
			Value=version,
			stringsAsFactors=FALSE, check.names=FALSE)
		write.table(myDataFrame, file=outputFile, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
		myDataFrame <- data.frame(
			Type="Run",
			SubType="-",
			Annotation="DSC Permutations",
			Value=thePerms,
			stringsAsFactors=FALSE, check.names=FALSE)
		write.table(myDataFrame, file=outputFile, append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
		myDataFrame <- data.frame(
			Type="Diagram",
			SubType="-",
			Annotation=c("Disp. Sep. Crit. (DSC)", "Disp. within groups (Dw)", "Disp. between groups (Db)", "DSC pvalue"),
			Value=c(theDscAllResults@mDSC, theDscAllResults@mDW, theDscAllResults@mDB, theDscAllResults@mPvalue),
			stringsAsFactors=FALSE, check.names=FALSE)
		write.table(myDataFrame, file=outputFile, append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	}
	##############################################
	### File for annotations related to a PCA components or component pairs. (DSC values for the components or component pairs can go here.)
	### PCAComponentAnnotations.tsv is a tab-delimited file containing three columns, Components, Annotation and Value.
	### Components contains either a single component number or a pair of component numbers separated by a comma.
	### Annotation is a display string for a value associated with a diagram and Value is the value to go with the display string.
	### Components	Annotation	Value
	if (!is.null(theDscResults))
	{
		compAcommaB <- paste("PC", theComponentA, ",", "PC", theComponentB, sep="")
		logDebug("writePCAAnnotations - data ", compAcommaB)
		myDataFrame <- data.frame(
				Type="Component",
				SubType=c(compAcommaB, compAcommaB, compAcommaB, compAcommaB),
				Annotation=c(paste("Disp. Sep. Crit. (DSC) (", theComponentA, ",", theComponentB, ")", sep=""),
						paste("Disp. within groups (Dw) (", theComponentA, ",", theComponentB, ")", sep=""),
						paste("Disp. between groups (Db) (", theComponentA, ",", theComponentB, ")", sep=""),
						paste("DSC pvalue(", theComponentA, ",", theComponentB, ")", sep="")),
				Value=c(theDscResults@mDSC, theDscResults@mDW, theDscResults@mDB, theDscResults@mPvalue),
				stringsAsFactors=FALSE)
		write.table(myDataFrame, file=outputFile, append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	}
}

openAndWriteIssuesLogFile<-function(theOutputDir)
{
	myFile <- file(cleanFilePath(theOutputDir, "error.log"), "w+")
	on.exit(close(myFile))
	cat("PCA not calculated or number of components calculated less than number requested\n", file=myFile, append=TRUE)
}

####################################################################
###
####################################################################
