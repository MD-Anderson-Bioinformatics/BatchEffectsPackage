#MBatch Copyright ? 2011, 2012, 2013, 2014, 2015, 2016, 2017 University of Texas MD Anderson Cancer Center
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

setClass("PCA-PVALUE-DSC", representation(
				mListOfGenePvalue="vector",
				mListOfGeneDSC="vector",
				mListOfGeneDB="vector",
				mListOfGeneDW="vector",
				mListOfResults="vector",
				mPvalue="character",
				mDSC="numeric",
				mDB="numeric",
				mDW="numeric"))

pvalueDSC <- function(thePca, theBatchIdsForSamples, thePermutations, theComponentA, theComponentB, theThreads)
{
	logDebug("pvalueDSC start")
	pcaDataExcerpt <- NULL
	if (theComponentA < 1)
	{
	  logDebug("pvalueDSC component A < 1")
	  pcaDataExcerpt <- t(thePca@scores)
	}
	else
	{
	  logDebug("pvalueDSC transpose")
	  pcaDataExcerpt <- t(thePca@scores[,((1:ncol(thePca@scores))==theComponentA)|((1:ncol(thePca@scores))==theComponentB)])
	}
	logDebug("pvalueDSC call pvalueDSCwithExcerpt")
	results <- pvalueDSCwithExcerpt(pcaDataExcerpt, theBatchIdsForSamples, thePermutations, theThreads)
	###logDebug(" DSC=", results[[1]], " dbValue=", results[[2]], " dwValue=", results[[3]], " pvalue=", results[[4]])
	return(results)
}

doDscPerms <- function(thePcaDataExcerpt, theBatchIdsForSamples, thePermutations, theThreads)
{
	dscJavaObj <- .jnew("org/mda/dscjava/DscJava")
	logDebug("doDscPerms before java")
	javaPcaDscObjList <- .jcall(dscJavaObj, "[Lorg/mda/dscjava/PcaDsc;", "doDscPerms", as.vector(thePcaDataExcerpt), dim(thePcaDataExcerpt), as.vector(theBatchIdsForSamples), as.integer(thePermutations), as.integer(theThreads))
	logDebug("doDscPerms after java")
	resultsList <- lapply(javaPcaDscObjList, function(javaPcaDscObj)
			{
				results <- new("PCA-DSC")
				results@mListOfDSCbyGene <-  .jcall(javaPcaDscObj, "[D", "getmListOfGeneDSC")
				results@mListOfDWbyGene <-  .jcall(javaPcaDscObj, "[D", "getmListOfGeneDW")
				results@mListOfDBbyGene <-  .jcall(javaPcaDscObj, "[D", "getmListOfGeneDB")
				results@mDSC <- .jcall(javaPcaDscObj, "D", "getmDSC")
				results@mDB <- .jcall(javaPcaDscObj, "D", "getmDB")
				results@mDW <- .jcall(javaPcaDscObj, "D", "getmDW")
				results
			})
	logDebug("doDscPerms length(resultsList)=", length(resultsList))
	return(resultsList)
}

pvalueDSCwithExcerpt <- function(thePcaDataExcerpt, theBatchIdsForSamples, thePermutations, theThreads)
{
  logDebug("pvalueDSCwithExcerpt start")
	### get initial DSC values
	DSCValues <- getDSCwithExcerpt(thePcaDataExcerpt, theBatchIdsForSamples)
	logDebug("pvalueDSCwithExcerpt after getDSCwithExcerpt")
	originalDscByGeneList <- DSCValues@mListOfDSCbyGene
	DSC <- DSCValues@mDSC
	dbValue <- DSCValues@mDB
	dwValue <- DSCValues@mDW
	###logDebug("DSC=", DSC, " dbValue=", dbValue, " dwValue=", dwValue)
	genePValues <- buildList(1, length(originalDscByGeneList))
	pValue <- 1
	permResults <- c("")
	logDebug("pvalueDSCwithExcerpt length(unique(theBatchIdsForSamples))=", length(unique(theBatchIdsForSamples)))
	if (length(unique(theBatchIdsForSamples))>1)
	{
		numGreater <- 0
		numPerms <- 0
		geneNumGreater <- buildList(0, length(originalDscByGeneList))
		permCount <- buildList(0, length(originalDscByGeneList))
		### this part can be multi-threaded if needed
		startTime <- Sys.time()
		logDebug("pvalueDSCwithExcerpt call doDscPerms")
		permResults <- doDscPerms(thePcaDataExcerpt, theBatchIdsForSamples, thePermutations, theThreads)
		logDebug("pvalueDSCwithExcerpt length(permResults)=", length(permResults))
		for(DSCPermObj in permResults)
		{
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
			###logDebug(" DSCPerm=", DSCPerm, " >= DSC=", DSC)

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
		logDebug("p-value run time = ", (endTime - startTime))
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
		logDebug("p-value is ", pValue)
	}
	logDebug("pvalueDSCwithExcerpt 2222 length(permResults)=", length(permResults))
	logDebug("pvalueDSCwithExcerpt PCA-PVALUE-DSC")
	results <- new("PCA-PVALUE-DSC")
	results@mListOfGenePvalue <-  genePValues
	results@mListOfGeneDSC <- originalDscByGeneList
	results@mListOfGeneDB <- DSCValues@mListOfDBbyGene
	results@mListOfGeneDW <- DSCValues@mListOfDWbyGene
	results@mListOfResults <- permResults
	results@mPvalue <-  pvalueZeroCheck(pValue, thePermutations)
	results@mDSC <- DSC
	results@mDB <- dbValue
	results@mDW <- dwValue
	return(results)
}

saveCompListDscData<-function(theFile, theList)
{
	### save value objects for later reuse in making DSC summary / passing in list of components / file is *__CompListDSC.RData
	save(theList, file=theFile)
}

saveCompDscData<-function(theFile, theObject, theComponentA, theComponentB)
{
	### save value objects for later reuse in making DSC summary / passing in PCA-PVALUE-DSC / file is *__CompDSC.RData
	save(theObject, theComponentA, theComponentB, file=theFile)
}

saveOverallDscData<-function(theFile, theOverallObject)
{
	### save value objects for later reuse in making DSC summary / passing in PCA-PVALUE-DSC / file is *__OverallDSC.RData
	save(theOverallObject, file=theFile)
}

####################################################################
###
####################################################################

