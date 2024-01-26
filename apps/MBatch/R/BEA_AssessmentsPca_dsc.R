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

setClass("PCA-DSC", representation(mListOfDSCbyGene="vector", mListOfDWbyGene="vector",mListOfDBbyGene="vector",mDSC="numeric",mDB="numeric",mDW="numeric"))


clearDSCOverviewFiles<-function(theStartDir)
{
	logDebug("clearDSCOverviewFiles - start")
	logDebug("clearDSCOverviewFiles - find all the *__CompListDSC.RData files")
	myListOfCompListDscFiles <- dir(path=theStartDir, pattern=glob2rx("*__CompListDSC.RData"), full.names = TRUE, recursive = TRUE)
	### for each *__CompListDSC.RData
	for (myCompListDscFile in myListOfCompListDscFiles)
	{
		###		find all the *__OverallDSC.RData files in the directories below the CompListDSC file
		currentDir <- dirname(myCompListDscFile)
		logDebug("clearDSCOverviewFiles - find all the *__OverallDSC.RData files")
		myOverallDscFileList <- dir(path=currentDir, pattern=glob2rx("*__OverallDSC.RData"), full.names = TRUE, recursive = TRUE)
		### 	For each *_OverallDSC.RData files
		for (myOverallDscFile in myOverallDscFileList)
		{
			###			find all the *__CompDSC.RData files in the same directory as the OverallDSC file
			logDebug("clearDSCOverviewFiles - find all the *__CompDSC.RData files")
			myCompDscFileList <- dir(path=dirname(myOverallDscFile), pattern=glob2rx("*__CompDSC.RData"), full.names = TRUE, recursive = FALSE)
			###			For each of the *__CompDSC.RData files
			for(myCompDscFile in myCompDscFileList)
			{
				###				delete the *__CompDSC.RData files
				unlink(myCompDscFile)
			}
			###			delete the *__OverallDSC.RData file
			unlink(myOverallDscFile)
		}
		###		delete the *__CompListDSC.RData file
		unlink(myCompListDscFile)
	}
	logDebug("clearDSCOverviewFiles - done")
}

getDSCSubPath<-function(theBaseDir, theCurrentDir)
{
  gsub(theBaseDir, "", theCurrentDir, fixed=TRUE)
}

dirAboveDSC <- function(theStartingDir)
{
  newDir <- dirname(theStartingDir)
  oldDir <- basename(theStartingDir)
  if ("DSC" != oldDir)
  {
    newDir <- dirAboveDSC(newDir)
  }
  return(newDir)
}

buildDSCOverviewFile<-function(theSearchDir, theOutputDir, theBaseDir, theOutputFile)
{
	### buildDSCOverviewFile
	########################################################################
	### first find all the *__CompListDSC.RData file (aka PCA directories)
	### for each *__CompListDSC.RData
	###		if structure for storing DSC Overview data is null
	### 		load the *__CompListDSC.RData file which has the list of components
	### 		get number of entries for overall DSC structure
	###			build structure for storing DSC Overview data
	###		delete the *__CompListDSC.RData file
	###		find all the *_OverallDSC.RData files in the directories below the CompListDSC file
	### 	For each *_OverallDSC.RData files
	###			Load the *__OverallDSC.RData file which has the overall DSC object
	###			add the new Overall DSC data into structure for storing DSC Overview data
	###			delete the *__OverallDSC.RData file
	###			find all the *__CompDSC.RData files in the same directory as the OverallDSC file
	###			For each of the *__CompDSC.RData files
	###				load the file
	###				add the component data into exising entry in structure for storing DSC Overview data
	###				delete the *__CompDSC.RData files
	### sort the structure for storing DSC Overview data
	### write DSC Overview file
	########################################################################
  logDebug("buildDSCOverviewFile theSearchDir=", theSearchDir)
  logDebug("buildDSCOverviewFile theOutputDir=", theOutputDir)
  logDebug("buildDSCOverviewFile theBaseDir=", theBaseDir)
  logDebug("buildDSCOverviewFile theOutputFile=", theOutputFile)
  #### flag for actually getting data to write
  foundDataFlag <- FALSE
  ####
	myDscOverviewStruct <- NULL
	### first find all the *__CompListDSC.RData file (aka PCA directories)
	myListOfCompListDscFiles <- dir(path=theSearchDir, pattern=glob2rx("*__CompListDSC.RData"), full.names = TRUE, recursive = TRUE)
	### for each *__CompListDSC.RData
	rowCounter <- 0
	logDebug("myCompListDscFile")
	for (myCompListDscFile in myListOfCompListDscFiles)
	{
	  logDebug("loop myCompListDscFile=", myCompListDscFile)
	  ###logDebug("loop myCompListDscFile")
		###		if structure for storing DSC Overview data is null
		if (is.null(myDscOverviewStruct))
		{
			logDebug("set struct")
			### "define" the variables being loaded by load
			theList <- NULL
			theOverallObject <- NULL
			theComponentA <- NULL
			theComponentB <- NULL
			theObject <- NULL
			### 		load the *__CompListDSC.RData file which has theList of components
			load(myCompListDscFile)
			logDebug("myCompListDscFile", myCompListDscFile)
			### theList has list of components
			columnNames <- c("dataset", "Overall-DSC", "Overall-DSC-pvalue")
			for(i in seq(from=1, to=(length(theList)-1), by=2 ))
			{
				componentA <- theList[[i]]
				componentB <- theList[[i+1]]
				columnNames <- c(columnNames, paste("DSC(", componentA, ",", componentB, ")", sep=""))
				columnNames <- c(columnNames, paste("DSC-pvalue(", componentA, ",", componentB, ")", sep=""))
				logDebug("comp A ", componentA, ", comp B", componentB)
			}
			### 		get number of entries for overall DSC structure
			numberOfRows <- length(dir(path=theSearchDir, pattern=glob2rx("*__OverallDSC.RData"), full.names = FALSE, recursive = TRUE))
			numberOfColumns <- length(columnNames)
			###			build structure for storing DSC Overview data
			### this has the columnNames across the top and numberOfEntries number of row entries, prepopulated with NA values
			if ((0!=numberOfRows)&&(0!=numberOfColumns))
			{
				myDscOverviewStruct <- matrixWithIssues(data = buildList(NA, (numberOfRows*numberOfColumns)), nrow=numberOfRows, ncol=numberOfColumns)
				colnames(myDscOverviewStruct) <- columnNames
			}
		}
		###		find all the *__OverallDSC.RData files in the directories below the CompListDSC file
		currentDir <- dirname(myCompListDscFile)
		myOverallDscFileList <- dir(path=currentDir, pattern=glob2rx("*__OverallDSC.RData"), full.names = TRUE, recursive = TRUE)
		### 	For each *_OverallDSC.RData files
		for (myOverallDscFile in myOverallDscFileList)
		{
			logDebug("myOverallDscFile ", myOverallDscFile)
			rowCounter <- rowCounter + 1
			logDebug("rowCounter ", rowCounter)
			###			Load the *__OverallDSC.RData file which has the overall DSC object
			### theOverallObject PCA-PVALUE-DSC mListOfGenePvalue, mListOfGeneDSC, mListOfGeneDB, mListOfGeneDW, mPvalue, mDSC, mDB, mDW
			load(myOverallDscFile)
			###			add the new Overall DSC data into structure for storing DSC Overview data
			### info from path
			dirpathForData <- dirname(myOverallDscFile)
			datasetDesc <- getDSCSubPath(theBaseDir, dirpathForData)
			logDebug("datasetDesc ", datasetDesc)
			#if (!is.null(theAltDatasetName))
			#{
  		#	  datasetDesc <- paste(theAltDatasetName, datasetDesc, sep="")
			#}
			### info from file
			overallDSC <- theOverallObject@mDSC
			overallPvalue <- theOverallObject@mPvalue
			###
			myDscOverviewStruct[rowCounter, "dataset"] <- datasetDesc
			myDscOverviewStruct[rowCounter, "Overall-DSC"] <- overallDSC
			myDscOverviewStruct[rowCounter, "Overall-DSC-pvalue"] <- overallPvalue
			###			find all the *__CompDSC.RData files in the same directory as the OverallDSC file
			myCompDscFileList <- dir(path=dirpathForData, pattern=glob2rx("*__CompDSC.RData"), full.names = TRUE, recursive = FALSE)
			###			For each of the *__CompDSC.RData files
			logDebug("length(myCompDscFileList) ", length(myCompDscFileList))
			if (length(myCompDscFileList)>0)
			{
  			for(myCompDscFile in myCompDscFileList)
  			{
  			  foundDataFlag <- TRUE
  				logDebug("myCompDscFile ", myCompDscFile)
  				###				load the file
  				load(myCompDscFile)
  				### theObject, theComponentA, theComponentB
  				dscColId <- paste("DSC(", theComponentA, ",", theComponentB, ")", sep="")
  				pvaColId <- paste("DSC-pvalue(", theComponentA, ",", theComponentB, ")", sep="")
  				###				add the component data into exising entry in structure for storing DSC Overview data
  				myDscOverviewStruct[rowCounter, dscColId] <- theObject@mDSC
  				myDscOverviewStruct[rowCounter, pvaColId] <- theObject@mPvalue
  				###logDebug("myCompDscFile values (", theComponentA, ",", theComponentB, ") DSC=", theObject@mDSC, " pvalue=", theObject@mPvalue)
  				###				delete the *__CompDSC.RData files
  				#--unlink(myCompDscFile)
  			}
			}
			###			delete the *__OverallDSC.RData file
			#--unlink(myOverallDscFile)
		}
		###		delete the *__CompListDSC.RData file
		#--unlink(myCompListDscFile)
	}
	logDebug("dim(myDscOverviewStruct) = ", dim(myDscOverviewStruct))
	logDebug("foundDataFlag = ", foundDataFlag)
	if ((!is.null(myDscOverviewStruct)) && (isTRUE(foundDataFlag)))
	{
	  outDscTSVpath <- cleanFilePath(theOutputDir, theOutputFile)
	  logDebug("overall DSC ", outDscTSVpath)
		### sort the structure for storing DSC Overview data
		oldColNames <- colnames(myDscOverviewStruct)
		myDscOverviewStruct <- data.frame(myDscOverviewStruct, check.names=FALSE)
		### myDscOverviewStruct <- myDscOverviewStruct[do.call("order", list(quote(myDscOverviewStruct))),]
		myDscOverviewStruct <- myDscOverviewStruct[do.call(order, myDscOverviewStruct),]
		myDscOverviewStruct <- as.matrix(myDscOverviewStruct)
		colnames(myDscOverviewStruct) <- oldColNames
		### write DSC Overview file
		checkDirForCreation(theOutputDir)
		write.table(myDscOverviewStruct, file=outDscTSVpath, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
	}
	return(myDscOverviewStruct)
}

getDSC <- function(thePca, theBatchIdsForSamples, theFirstComponent, theSecondComponent)
{
	### pca@scores[sample,componentId]
	pcaDataExcerpt <- t(thePca@scores[,((1:ncol(thePca@scores))==theFirstComponent)|((1:ncol(thePca@scores))==theSecondComponent)])
	results <- getDSCwithExcerpt(pcaDataExcerpt, theBatchIdsForSamples)
	###logDebug(" DSCij=", results[[1]], " Dwij=", results[[2]], " Dbij=", results[[3]])
	return(results)
}

getDSCwithExcerpt <- function(thePcaDataExcerpt, theBatchIdsForSamples)
{
  logDebug("nrow(thePcaDataExcerpt)=", nrow(thePcaDataExcerpt))
  logDebug("ncol(thePcaDataExcerpt)=", ncol(thePcaDataExcerpt))
  logDebug("length(theBatchIdsForSamples)=", length(theBatchIdsForSamples))
  stopifnotWithLogging("Number of batch ids should match number of samples in pca data",
			length(theBatchIdsForSamples)==length(thePcaDataExcerpt[1,]),length(theBatchIdsForSamples)==length(thePcaDataExcerpt[2,]))
	# title not used anymore
	logDebug("getDSCwithExcerpt before Python")
	# logDebug("getDSCwithExcerpt - use_condaenv = ", getGlobalMBatchEnv())
	# use_condaenv(getGlobalMBatchEnv())
	logDebug("getDSCwithExcerpt - import(mbatch.dsc.dsc_calc)")
	calc <- import("mbatch.dsc.dsc_calc")
	# samples need to be cast as character for np_array to work
	batchIdsForSamples <- as.character(as.vector(unlist(theBatchIdsForSamples)))
	names(batchIdsForSamples) <- colnames(thePcaDataExcerpt)
	# matrix needs to be data.frame for Python/Reticulate
	# same for theBatchIdsForSamples
	calcInfo <- calc$dsc_calc(r_to_py(thePcaDataExcerpt), np_array(batchIdsForSamples), TRUE)
	logDebug("getDSCwithExcerpt after Python")
	results <- new("PCA-DSC")
	results@mListOfDSCbyGene <- calcInfo$m_list_of_feature_dsc
	results@mListOfDWbyGene <-  calcInfo$m_list_of_feature_dw
	results@mListOfDBbyGene <-  calcInfo$m_list_of_feature_db
	results@mDSC <- calcInfo$m_dsc
	results@mDB <-  calcInfo$m_db
	results@mDW <-  calcInfo$m_dw
	return(results)
}

####################################################################
###
####################################################################
