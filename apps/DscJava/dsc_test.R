#MBatch Copyright 2011, 2012 University of Texas MD Anderson Cancer Center
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(snowfall, warn.conflicts=FALSE, verbose=FALSE)

baseDir <-  "/tmp/parallel_tests"

myScores <- file.path(baseDir, "PCAScores.tsv")
myBatches <- file.path(baseDir, "BatchData.tsv")
mysource <- file.path(file.path(baseDir, "dsc-r"), "dsc_test.R")

################################################################################

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

################################################################################

buildList<-function(theValue, theLength)
{
	### looks weird so that it works with NA and NAN
	return(as.vector(unlist(lapply(c(1:theLength), function(x) {return(theValue)}))))
}

################################################################################

####################################################################
### 
####################################################################

stopifnotWithLogging<-function(msg="", ... )
{
	if (sum(...)!=length(c(...))) 
	{
		message(msg)
		stop(msg, call.=FALSE)
	}
}

####################################################################
### 
####################################################################

setClass("PCA-PVALUE-DSC", representation(
				mListOfGenePvalue="vector", 
				mListOfGeneDSC="vector",
				mListOfGeneDB="vector",
				mListOfGeneDW="vector",
				mPvalue="numeric",
				mDSC="numeric",
				mDB="numeric",
				mDW="numeric"))

pvalueDSC <- function(thePcaScores, theBatchIdsForSamples, thePermutations, theComponentA, theComponentB, theThreads)
{
	###message("pvalueDSC")
	pcaDataExcerpt <- NULL
	if (theComponentA < 1)
	{
		pcaDataExcerpt <- t(thePcaScores)
	}
	else
	{
		pcaDataExcerpt <- t(thePcaScores[,((1:ncol(thePcaScores))==theComponentA)|((1:ncol(thePcaScores))==theComponentB)])
	}
	results <- pvalueDSCwithExcerpt(pcaDataExcerpt, theBatchIdsForSamples, thePermutations, theThreads)
	###message(" DSC=", results[[1]], " dbValue=", results[[2]], " dwValue=", results[[3]], " pvalue=", results[[4]])
	return(results)
}

doDscPerms <- function(theIgnoredIndex)
{					
	### permute the data randomly, keeping sizes of batches the same
	permData <- t(apply(thePcaDataExcerpt, 1, function(x) 
					{
						###message("old= ", paste(x, collapse=", "))
						sx <- sample(x, length(x))
						###message("new= ", paste(x, collapse=", "))
						return(sx)
					}))
	DSCPermObj <- getDSCwithExcerpt(permData, theBatchIdsForSamples)
}

pvalueDSCwithExcerpt <- function(thePcaDataExcerpt, theBatchIdsForSamples, thePermutations, theThreads)
{
	###message("pvalueDSCwithExcerpt")
	### get initial DSC values
	DSCValues <- getDSCwithExcerpt(thePcaDataExcerpt, theBatchIdsForSamples)
	originalDscByGeneList <- DSCValues@mListOfDSCbyGene
	DSC <- DSCValues@mDSC
	dbValue <- DSCValues@mDB
	dwValue <- DSCValues@mDW
	###message("DSC=", DSC, " dbValue=", dbValue, " dwValue=", dwValue)
	genePValues <- buildList(1, length(originalDscByGeneList))
	pValue <- 1
	###message("theThreads=", theThreads, " thePermutations=", thePermutations)
	sfInit(parallel=TRUE, cpus=theThreads)
	#message("on.exit")
	on.exit(sfStop(), add=TRUE)
	#message("mysource=", mysource)
	sfSource(mysource)
	#message("thePcaDataExcerpt")
	sfExport("thePcaDataExcerpt")
	#message("theBatchIdsForSamples")
	sfExport("theBatchIdsForSamples")
	#message("if")
	if (length(unique(theBatchIdsForSamples))>1)
	{
		numGreater <- 0
		numPerms <- 0
		geneNumGreater <- buildList(0, length(originalDscByGeneList))
		permCount <- buildList(0, length(originalDscByGeneList))
		### this part can be multi-threaded if needed
		startTime <- Sys.time()
		permResults <- sfLapply(c(1:thePermutations), function(x) { doDscPerms() })
		for(DSCPermObj in permResults)
		{
			###message("Perm results.mDSC =" , DSCPermObj@mDSC, " results.mDB =" , DSCPermObj@mDB, " results.mDW =" , DSCPermObj@mDW);
			DSCPerm <- DSCPermObj@mDSC
			geneDSCPerm <- DSCPermObj@mListOfDSCbyGene
			for (myIndex in 1:length(geneNumGreater))
			{
				if (!is.nan(geneDSCPerm[[myIndex]]))
				{
					if (!is.nan(originalDscByGeneList[[myIndex]]))
					{
						permCount[myIndex] <- permCount[[myIndex]] + 1
						if (geneDSCPerm[[myIndex]] >= originalDscByGeneList[[myIndex]])
						{
							geneNumGreater[myIndex] <- geneNumGreater[[myIndex]] + 1
						}
					}
				}
			}
			### Trace_SbPerm and Trace_SwPerm are not used
			###message(" DSCPerm=", DSCPerm, " >= DSC=", DSC)
			
			if(!is.nan(DSCPerm))
			{
				if(!is.nan(DSC))
				{
					numPerms <- numPerms + 1 
					if(DSCPerm >= DSC)
					{
						numGreater <- numGreater + 1
					}
				}
			}
		}
		endTime <- Sys.time()
		###message("p-value run time = ", (endTime - startTime))
		### basic pvalue
		pValue <- 0
		if(numPerms > 0)
		{
			pValue <- numGreater / numPerms
		}
		### gene pvalues
		genePValues <- buildList(0, length(geneNumGreater))
		if(thePermutations > 0)
		{
			for (myIndex in 1:length(genePValues))
			{
				if (is.nan(originalDscByGeneList[[myIndex]]))
				{
					genePValues[myIndex] <- NaN
				}
				else
				{
					genePValues[myIndex] <- geneNumGreater[[myIndex]]/permCount[[myIndex]]
				}
			}
		}
		###message("p-value is ", pValue)
	}
	results <- new("PCA-PVALUE-DSC")
	results@mListOfGenePvalue <-  genePValues
	results@mListOfGeneDSC <- originalDscByGeneList
	results@mListOfGeneDB <- DSCValues@mListOfDBbyGene
	results@mListOfGeneDW <- DSCValues@mListOfDWbyGene
	results@mPvalue <-  pValue
	results@mDSC <- DSC
	results@mDB <- dbValue
	results@mDW <- dwValue
	return(results)
}

####################################################################
### 
####################################################################


setClass("PCA-DSC", representation(
				mListOfDSCbyGene="vector", 
				mListOfDWbyGene="vector", 
				mListOfDBbyGene="vector", 
				mDSC="numeric",
				mDB="numeric",
				mDW="numeric"))

getDSC <- function(thePca, theBatchIdsForSamples, theFirstComponent, theSecondComponent)
{
	### pca@scores[sample,componentId]
	pcaDataExcerpt <- t(thePca@scores[,((1:ncol(thePca@scores))==theFirstComponent)|((1:ncol(thePca@scores))==theSecondComponent)]) 
	results <- getDSCwithExcerpt(pcaDataExcerpt, theBatchIdsForSamples)
	###message(" DSCij=", results[[1]], " Dwij=", results[[2]], " Dbij=", results[[3]])
	return(results)
}

getDSCwithExcerpt <- function(thePcaDataExcerpt, theBatchIdsForSamples)
{
	###message("length(theBatchIdsForSamples)=",length(theBatchIdsForSamples))
	###message("length(thePcaDataExcerpt[1,])=",length(thePcaDataExcerpt[1,]))
	###message("length(thePcaDataExcerpt[2,])=",length(thePcaDataExcerpt[2,]))
	stopifnotWithLogging("Number of batch ids should match number of samples in pca data",
			length(theBatchIdsForSamples)==length(thePcaDataExcerpt[1,]),length(theBatchIdsForSamples)==length(thePcaDataExcerpt[2,]))
	### theBatchIdsForSamples - this is part of the @scores attribute from the pca object which now contains only the two components specified for this diagram
	### columns are samples - rows are components
	
	### mean of across all the samples in the dataset
	###message("mean of across all the samples in the dataset")
	meanListAcrossSamples <- rowMeans(thePcaDataExcerpt, na.rm = TRUE)
	###message("meanListAcrossSamples=", paste(meanListAcrossSamples, collapse=", ", sep=", "))
	###message("number of samples")
	numberOfSamples <- length(theBatchIdsForSamples)
	###message("numberOfSamples=", numberOfSamples)
	dwValue <- 0
	dbValue <- 0
	
	###message("batch id list")
	batchIdList <- sort(unique(sort(theBatchIdsForSamples)))
	###message("batchIdList=", batchIdList)
	gDSC <- buildList(NA, nrow(thePcaDataExcerpt))
	gDW <- buildList(NA, nrow(thePcaDataExcerpt))
	gDB <- buildList(NA, nrow(thePcaDataExcerpt))
	### for each gene
	for (compRow in 1:nrow(thePcaDataExcerpt))
	{
		geneDw <- 0
		geneDb <- 0
		### for each batch
		for (batchName in batchIdList)
		{
			###message("batchName=", batchName)
			### nb <- number of samples in batch
			###message("number of samples in batch")
			numSamplesInBatch <- sum(theBatchIdsForSamples==batchName)
			###message("batchName=", batchName, "numSamplesInBatch=", numSamplesInBatch)
			###message("values for batch")
			valuesForBatch <- thePcaDataExcerpt[compRow,theBatchIdsForSamples==batchName]
			###message("valuesForBatch=", paste(valuesForBatch, collapse=", "))
			### mean across samples in batch b only
			###message("mean across samples in batch b only")
			batchMean <- mean(valuesForBatch, na.rm = TRUE)
			### calculate sample variance across samples 
			###message("calculate sample variance across samples")
			sampleVariance <- var(valuesForBatch, na.rm = TRUE)
			###message("batchMean=", batchMean, "sampleVariance=", sampleVariance)
			### use sample variance to get population variance
			###message("use sample variance to get population variance")
			geneDw <- geneDw + ((numSamplesInBatch-1) / numberOfSamples * sampleVariance)
			geneDb <- geneDb + (
						(numSamplesInBatch / numberOfSamples) * ( batchMean - meanListAcrossSamples[[compRow]] ) * ( batchMean - meanListAcrossSamples[[compRow]] )
						)
			options(digits=22)
			###message("b=", batchName, " nsb=", numSamplesInBatch, " ns=", numberOfSamples, " bm=", batchMean, " sv=", sampleVariance, " geneDw=", geneDw, " geneDb=", geneDb, " m=", meanListAcrossSamples[[compRow]])
		}
		dwValue <- dwValue + geneDw
		dbValue <- dbValue + geneDb
		### take square root
		geneDw <- sqrt(geneDw)
		geneDb <- sqrt(geneDb)
		#
		geneDSC <- NaN
		if ((!is.na(geneDw))&&(geneDw > 0))
		{
			geneDSC <- geneDb/geneDw
		}
		gDSC[compRow] <- geneDSC
		gDW[compRow] <- geneDw
		gDB[compRow] <- geneDb
	}
	### take square root
	dwValue <- sqrt(dwValue)
	dbValue <- sqrt(dbValue)
	### avoid divide by zero error
	###message("avoid divide by zero error")
	###message("dwValue ", dwValue)
	DSC <- NaN
	if ((!is.na(dwValue))&&(dwValue > 0))
	{
		DSC <- dbValue/dwValue
	}
	if (length(unique(theBatchIdsForSamples))<2)
	{
		gDSC <- buildList(0, nrow(thePcaDataExcerpt))
		gDB <- buildList(0, nrow(thePcaDataExcerpt))
		DSC <- 0
		dbValue <- 0
	}
	results <- new("PCA-DSC")
	results@mListOfDSCbyGene <-  lapply(gDSC, epsilonZeroCheck)
	results@mListOfDWbyGene <-  gDW
	results@mListOfDBbyGene <-  gDB
	results@mDSC <- epsilonZeroCheck(DSC)
	results@mDB <- dbValue
	results@mDW <- dwValue
	###message("RESULTS results.mDSC =" , results@mDSC);
	###message("RESULTS results.mDB =" , results@mDB);
	###message("RESULTS results.mDW =" , results@mDW);
	return(results)
}

####################################################################
### 
####################################################################

doRTest <- function(theThreads)
{
	pcascores <- t(read.table(myScores, row.names=1, sep="\t", na.strings="NA", header=TRUE, stringsAsFactors=FALSE, quote="", check.names=FALSE))
	batchInfo <- read.table(myBatches, sep="\t", na.strings="NA", header=TRUE, stringsAsFactors=FALSE, quote="", colClasses=c("character"), as.is=TRUE, check.names=FALSE)
	indexList <- !is.na(match(batchInfo$Sample, rownames(pcascores)))
	batchIdsForSamples <- batchInfo$BatchId[indexList]
	myDSCPermutations <- 1000
	myDSCThreads <- theThreads
	dscAllResults <- NULL
	time <- system.time(dscAllResults <- pvalueDSC(pcascores, batchIdsForSamples, myDSCPermutations, 0, 0, myDSCThreads))
	###dscCompResults <- pvalueDSC(pcascores, batchIdsForSamples, myDSCPermutations, 1, 2, myDSCThreads)
	message("All")
	#message("mListOfGenePvalue=",paste(dscAllResults@mListOfGenePvalue, collapse=","))
	#message("mListOfGeneDSC=",paste(dscAllResults@mListOfGeneDSC, collapse=","))
	#message("mListOfGeneDB=",paste(dscAllResults@mListOfGeneDB, collapse=","))
	#message("mListOfGeneDW=",paste(dscAllResults@mListOfGeneDW, collapse=","))
	message("mPvalue=",dscAllResults@mPvalue, " should be ", "0")
	message("mDSC=",dscAllResults@mDSC, " should be ", "0.523")
	message("mDB=",dscAllResults@mDB, " should be ", "33.837")
	message("mDW=",dscAllResults@mDW, " should be ", "64.659")
	#message("Comp1, Comp2")
	#message("mListOfGenePvalue=",paste(dscCompResults@mListOfGenePvalue, collapse=","))
	#message("mListOfGeneDSC=",paste(dscCompResults@mListOfGeneDSC, collapse=","))
	#message("mListOfGeneDB=",paste(dscCompResults@mListOfGeneDB, collapse=","))
	#message("mListOfGeneDW=",paste(dscCompResults@mListOfGeneDW, collapse=","))
	###message("mPvalue=",dscCompResults@mPvalue, " should be ", "0")
	###message("mDSC=",dscCompResults@mDSC, " should be ", "0.922")
	###message("mDB=",dscCompResults@mDB, " should be ", "25.619")
	###message("mDW=",dscCompResults@mDW, " should be ", "27.793")
	time
}



