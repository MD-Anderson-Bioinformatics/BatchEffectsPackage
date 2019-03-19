#MBatch Copyright ? 2011, 2012, 2013, 2014, 2015, 2016, 2017 University of Texas MD Anderson Cancer Center
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
	
buildDSCOverviewFile<-function(theStartDir, theOutputFile)
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
	
	myDscOverviewStruct <- NULL
	### first find all the *__CompListDSC.RData file (aka PCA directories)
	myListOfCompListDscFiles <- dir(path=theStartDir, pattern=glob2rx("*__CompListDSC.RData"), full.names = TRUE, recursive = TRUE)
	### for each *__CompListDSC.RData
	rowCounter <- 0
	logDebug("myCompListDscFile")
	for (myCompListDscFile in myListOfCompListDscFiles)
	{
		###logDebug("loop myCompListDscFile")
		###		if structure for storing DSC Overview data is null
		if (is.null(myDscOverviewStruct))
		{
			###logDebug("set struct")
			### "define" the variables being loaded by load
			theList <- NULL
			theOverallObject <- NULL
			theComponentA <- NULL
			theComponentB <- NULL
			theObject <- NULL
			### 		load the *__CompListDSC.RData file which has theList of components
			load(myCompListDscFile)
			###logDebug("myCompListDscFile", myCompListDscFile)
			### theList has list of components
			columnNames <- c("run-date", "disease-type", "data-type", "platform", "data-level", "correction-type", "PCA", "batch-type", "diagram-batch", "Overall-DSC", "Overall-DSC-pvalue")
			for(i in seq(from=1, to=(length(theList)-1), by=2 ))
			{
				componentA <- theList[[i]]
				componentB <- theList[[i+1]]
				columnNames <- c(columnNames, paste("DSC(", componentA, ",", componentB, ")", sep="")) 
				columnNames <- c(columnNames, paste("DSC-pvalue(", componentA, ",", componentB, ")", sep=""))
				###logDebug("comp A ", componentA, ", comp B", componentB)
			}
			### 		get number of entries for overall DSC structure
			numberOfRows <- length(dir(path=theStartDir, pattern=glob2rx("*__OverallDSC.RData"), full.names = FALSE, recursive = TRUE))
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
			###logDebug("myOverallDscFile ", myOverallDscFile)
			rowCounter <- rowCounter + 1
			###logDebug("rowCounter ", rowCounter)
			###			Load the *__OverallDSC.RData file which has the overall DSC object
			### theOverallObject PCA-PVALUE-DSC mListOfGenePvalue, mListOfGeneDSC, mListOfGeneDB, mListOfGeneDW, mPvalue, mDSC, mDB, mDW
			load(myOverallDscFile)
			###			add the new Overall DSC data into structure for storing DSC Overview data
			### info from path
			dirpathForData <- dirname(myOverallDscFile)
			diagramBatch <- basename(dirpathForData)
			dirpathForData <- dirname(dirpathForData)
			batchtype <- basename(dirpathForData)
			dirpathForData <- dirname(dirpathForData)
			PCA <- basename(dirpathForData)
			dirpathForData <- dirname(dirpathForData)
			correctiontype <- basename(dirpathForData)
			dirpathForData <- dirname(dirpathForData)
			datalevel <- basename(dirpathForData)
			dirpathForData <- dirname(dirpathForData)
			platform <- basename(dirpathForData)
			dirpathForData <- dirname(dirpathForData)
			datatype <- basename(dirpathForData)
			dirpathForData <- dirname(dirpathForData)
			diseasetype <- basename(dirpathForData)
			dirpathForData <- dirname(dirpathForData)
			rundate <- basename(dirpathForData)
			### info from file
			overallDSC <- theOverallObject@mDSC
			overallPvalue <- theOverallObject@mPvalue
			###
			myDscOverviewStruct[rowCounter, "run-date"] <- rundate  
			myDscOverviewStruct[rowCounter, "disease-type"] <- diseasetype  
			myDscOverviewStruct[rowCounter, "data-type"] <- datatype  
			myDscOverviewStruct[rowCounter, "platform"] <- platform  
			myDscOverviewStruct[rowCounter, "data-level"] <- datalevel  
			myDscOverviewStruct[rowCounter, "correction-type"] <- correctiontype  
			myDscOverviewStruct[rowCounter, "PCA"] <- PCA  
			myDscOverviewStruct[rowCounter, "batch-type"] <- batchtype 
			myDscOverviewStruct[rowCounter, "diagram-batch"] <- diagramBatch  
			myDscOverviewStruct[rowCounter, "Overall-DSC"] <- overallDSC  
			myDscOverviewStruct[rowCounter, "Overall-DSC-pvalue"] <- overallPvalue  
			###			find all the *__CompDSC.RData files in the same directory as the OverallDSC file 
			myCompDscFileList <- dir(path=dirname(myOverallDscFile), pattern=glob2rx("*__CompDSC.RData"), full.names = TRUE, recursive = FALSE)
			###			For each of the *__CompDSC.RData files
			for(myCompDscFile in myCompDscFileList)
			{
				###logDebug("myCompDscFile ", myCompDscFile)
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
			###			delete the *__OverallDSC.RData file 
			#--unlink(myOverallDscFile)
		}
		###		delete the *__CompListDSC.RData file
		#--unlink(myCompListDscFile)
	}
	if (!is.null(myDscOverviewStruct))
	{
		### sort the structure for storing DSC Overview data
		oldColNames <- colnames(myDscOverviewStruct)
		myDscOverviewStruct <- data.frame(myDscOverviewStruct, check.names=FALSE)
		### myDscOverviewStruct <- myDscOverviewStruct[do.call("order", list(quote(myDscOverviewStruct))),]
		myDscOverviewStruct <- myDscOverviewStruct[do.call(order, myDscOverviewStruct),]
		myDscOverviewStruct <- as.matrix(myDscOverviewStruct)
		colnames(myDscOverviewStruct) <- oldColNames
		### write DSC Overview file
		write.table(myDscOverviewStruct, file=file.path(theStartDir, theOutputFile), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
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
	stopifnotWithLogging("Number of batch ids should match number of samples in pca data",
			length(theBatchIdsForSamples)==length(thePcaDataExcerpt[1,]),length(theBatchIdsForSamples)==length(thePcaDataExcerpt[2,]))
	logDebug("nrow(thePcaDataExcerpt)=", nrow(thePcaDataExcerpt))
	logDebug("ncol(thePcaDataExcerpt)=", ncol(thePcaDataExcerpt))
	logDebug("length(theBatchIdsForSamples)=", length(theBatchIdsForSamples))
	dscJavaObj <- .jnew("org/mda/dscjava/DscJava")
	logDebug("getDSCwithExcerpt before java")
	javaPcaDscObj <- .jcall(dscJavaObj, "Lorg/mda/dscjava/PcaDsc;", "getDSCwithExcerpt", 
													as.vector(thePcaDataExcerpt), 
													dim(thePcaDataExcerpt), 
													as.vector(theBatchIdsForSamples))
	logDebug("getDSCwithExcerpt after java")
	results <- new("PCA-DSC")
	results@mListOfDSCbyGene <-  .jcall(javaPcaDscObj, "[D", "getmListOfGeneDSC")
	results@mListOfDWbyGene <-  .jcall(javaPcaDscObj, "[D", "getmListOfGeneDW")
	results@mListOfDBbyGene <-  .jcall(javaPcaDscObj, "[D", "getmListOfGeneDB")
	results@mDSC <- .jcall(javaPcaDscObj, "D", "getmDSC")
	results@mDB <- .jcall(javaPcaDscObj, "D", "getmDB")
	results@mDW <- .jcall(javaPcaDscObj, "D", "getmDW")
	return(results)
}

####################################################################
### 
####################################################################
