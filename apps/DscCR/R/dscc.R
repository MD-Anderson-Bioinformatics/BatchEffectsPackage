#DscCR Copyright 2023 University of Texas MD Anderson Cancer Center
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

################################################################################


####################################################################
###
####################################################################

getTestInputDirDscCR <- function()
{
  value <- Sys.getenv("DSCCR_TEST_INPUT")
  if (!isTRUE(file.exists(value)))
  {
    value <- "/BatchEffectsPackage_data/testing_static/DSCCR"
  }
  value
}

getTestOutputDirDscCR <- function()
{
  value <- Sys.getenv("DSCCR_TEST_OUTPUT")
  if (!isTRUE(file.exists(value)))
  {
    value <- "/BatchEffectsPackage_data/testing_dynamic/DSCCR"
  }
  value
}

getTestCompareDirDscCR <- function()
{
  value <- Sys.getenv("DSCCR_TEST_COMPARE")
  if (!isTRUE(file.exists(value)))
  {
    value <- "/BatchEffectsPackage_data/testing_static/DSCCR"
  }
  value
}

cleanFilePath <- function(theDir, theNewDir)
{
  #cleanDir <- iconv(theNewDir, from = "UTF-8", to = "ASCII", sub = "byte")
  #cleanDir <- gsub("><", "_", cleanDir, fixed=TRUE)
  #cleanDir <- gsub("<", "_", cleanDir, fixed=TRUE)
  #cleanDir <- gsub(">", "_", cleanDir, fixed=TRUE)
  #file.path(theDir, cleanDir)
  tmpPath <- file.path(theDir, theNewDir)
  # clean double // from the path
  tmpPath <- gsub("//", "/", tmpPath)
  tmpPath
}

####################################################################
###
####################################################################

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

pvalueDSC_R <- function(thePcaScores, theBatchIdsForSamples, thePermutations, theComponentA, theComponentB, theThreads, theRandomSeed)
{
  pvalueDSC_scores(thePcaScores, theBatchIdsForSamples, thePermutations, theComponentA, theComponentB, theThreads, theRandomSeed)
}

pvalueDSC_C <- function(thePcaScores, theBatchIdsForSamples, thePermutations, theComponentA, theComponentB, theThreads)
{
  # TODO: add passing in SEED
  ###message("pvalueDSC_C")
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

doDscPerms <- function(thePcaDataExcerpt, theBatchIdsForSamples, thePermutations, theThreads)
{
  dsccDoDscPerms(thePcaDataExcerpt, dim(thePcaDataExcerpt), theBatchIdsForSamples, thePermutations, theThreads)
  #	resultsList <- lapply(javaPcaDscObjList, function(javaPcaDscObj)
  #	{
  #		results <- new("PCA-DSC")
  #		results@mListOfDSCbyGene <-  .jcall(javaPcaDscObj, "[D", "getmListOfGeneDSC")
  #		results@mListOfDWbyGene <-  .jcall(javaPcaDscObj, "[D", "getmListOfGeneDW")
  #		results@mListOfDBbyGene <-  .jcall(javaPcaDscObj, "[D", "getmListOfGeneDB")
  #		results@mDSC <- .jcall(javaPcaDscObj, "D", "getmDSC")
  #		results@mDB <- .jcall(javaPcaDscObj, "D", "getmDB")
  #		results@mDW <- .jcall(javaPcaDscObj, "D", "getmDW")
  #		return(results)
  #	})
  #	return(resultsList);
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
  if (length(unique(theBatchIdsForSamples))>1)
  {
    numGreater <- 0
    numPerms <- 0
    geneNumGreater <- buildList(0, length(originalDscByGeneList))
    permCount <- buildList(0, length(originalDscByGeneList))
    ### this part can be multi-threaded if needed
    startTime <- Sys.time()
    permResults <- doDscPerms(thePcaDataExcerpt, theBatchIdsForSamples, thePermutations, theThreads)
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
  ###message(paste(dim(thePcaDataExcerpt), collapse=","))
  results <- dsccGetDSCwithExcerpt(thePcaDataExcerpt, dim(thePcaDataExcerpt), theBatchIdsForSamples)
  #	results <- new("PCA-DSC")
  #	results@mListOfDSCbyGene <-  .jcall(javaPcaDscObj, "[D", "getmListOfGeneDSC")
  #	results@mListOfDWbyGene <-  .jcall(javaPcaDscObj, "[D", "getmListOfGeneDW")
  #	results@mListOfDBbyGene <-  .jcall(javaPcaDscObj, "[D", "getmListOfGeneDB")
  #	results@mDSC <- .jcall(javaPcaDscObj, "D", "getmDSC")
  #	results@mDB <- .jcall(javaPcaDscObj, "D", "getmDB")
  #	results@mDW <- .jcall(javaPcaDscObj, "D", "getmDW")
  #	message("RESULTS results.mDSC =" , results@mDSC);
  #	message("RESULTS results.mDB =" , results@mDB);
  #	message("RESULTS results.mDW =" , results@mDW);
  return(results)
}
