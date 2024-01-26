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

pvalueDSC_scores <- function(thePcaScores, theBatchIdsForSamples, thePermutations, theComponentA, theComponentB, theThreads, theSeed)
{
  logDebug("pvalueDSC_scores start")
  stopifnotWithLogging("Component A should be a number", is.numeric(theComponentA))
  stopifnotWithLogging("Component B should be a number", is.numeric(theComponentB))
  pcaDataExcerpt <- NULL
  if (theComponentA < 1)
  {
    logDebug("pvalueDSC_scores component A < 1")
    pcaDataExcerpt <- t(thePcaScores)
  }
  else
  {
    logDebug("pvalueDSC_scores transpose")
    pcaDataExcerpt <- t(thePcaScores[,((1:ncol(thePcaScores))==theComponentA)|((1:ncol(thePcaScores))==theComponentB)])
  }
  logDebug("pvalueDSC_scores call pvalueDSCwithExcerpt")
  results <- pvalueDSCwithExcerpt(pcaDataExcerpt, theBatchIdsForSamples, theSeed, thePermutations, theThreads)
  ###logDebug(" DSC=", results[[1]], " dbValue=", results[[2]], " dwValue=", results[[3]], " pvalue=", results[[4]])
  return(results)
}

pvalueDSC <- function(thePca, theBatchIdsForSamples, thePermutations, theComponentA, theComponentB, theThreads, theSeed)
{
	logDebug("pvalueDSC start")
  stopifnotWithLogging("Component A should be a number", is.numeric(theComponentA))
  stopifnotWithLogging("Component B should be a number", is.numeric(theComponentB))
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
	results <- pvalueDSCwithExcerpt(pcaDataExcerpt, theBatchIdsForSamples, theSeed, thePermutations, theThreads)
	###logDebug(" DSC=", results[[1]], " dbValue=", results[[2]], " dwValue=", results[[3]], " pvalue=", results[[4]])
	return(results)
}

doDscPerms <- function(thePcaDataExcerpt, theBatchIdsForSamples, theSeed, thePermutations, theThreads)
{
  # set names for list of batch ids
  # special build of list for python
  python_vector <- as.character(as.vector(unlist(theBatchIdsForSamples)))
  names(python_vector) <- colnames(thePcaDataExcerpt)
  # use_condaenv(getGlobalMBatchEnv())
  logDebug("doDscPerms - import(mbatch.dsc.dsc_perm)")
  calc <- import("mbatch.dsc.dsc_perm")
  # the_df: pandas.DataFrame, the_batches: pandas.Series, the_seed: int, the_threads: int
  logDebug("doDscPerms thePcaDataExcerpt dim")
  logDebug(dim(thePcaDataExcerpt))
  logDebug("doDscPerms python_vector length")
  logDebug(length(python_vector))
  logDebug("doDscPerms theSeed")
  logDebug(theSeed)
  logDebug("doDscPerms thePermutations")
  logDebug(thePermutations)
  logDebug("doDscPerms theThreads")
  logDebug(theThreads)
  # as.data.frame needed for proper R to Python conversion
  # same with as.integer
  logDebug("doDscPerms before Python")
  calcInfoList <- calc$dsc_perm_calc_count(r_to_py(thePcaDataExcerpt),
                                           np_array(python_vector),
                                           as.integer(theSeed),
                                           as.integer(thePermutations),
                                           as.integer(theThreads))
  logDebug("doDscPerms after Python")
	resultsList <- lapply(calcInfoList, function(pythonDscObj)
	  {
  	  results <- new("PCA-DSC")
  	  results@mListOfDSCbyGene <- pythonDscObj$m_list_of_feature_dsc
  	  results@mListOfDWbyGene <-  pythonDscObj$m_list_of_feature_dw
  	  results@mListOfDBbyGene <-  pythonDscObj$m_list_of_feature_db
  	  results@mDSC <- pythonDscObj$m_dsc
  	  results@mDB <-  pythonDscObj$m_db
  	  results@mDW <-  pythonDscObj$m_dw
			results
		})
	logDebug("doDscPerms length(resultsList)=", length(resultsList))
	return(resultsList)
}

pvalueDSCwithExcerpt <- function(thePcaDataExcerpt, theBatchIdsForSamples, theSeed, thePermutations, theThreads)
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
		permResults <- doDscPerms(thePcaDataExcerpt, theBatchIdsForSamples, theSeed, thePermutations, theThreads)
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

