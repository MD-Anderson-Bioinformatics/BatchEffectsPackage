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

###############################################################################
###############################################################################

EB_internal<-function(theBeaData,
                      theBatchIdsNotToCorrect,
                      theDoCheckPlotsFlag,
                      theBatchType,
                      theEbType,
                      thePriorFlag,
                      theThreads=1,
                      thePath=NULL,
                      theWriteToFile=FALSE)
{
  logInfo("EB_internal - starting")
  results <- NULL
  if ((1==length(theBatchIdsNotToCorrect))&&(""==theBatchIdsNotToCorrect[[1]]))
  {
    theBatchIdsNotToCorrect<-NULL
  }
  myPath <- thePath
  priorFile <- "notused.png"
  myErrorFile <- NULL
  #setLogging(new("Logging", theFile=cleanFilePath(myPath, "EB.log")))
  if (TRUE==theWriteToFile)
  {
    checkDirForCreation(myPath)
    priorFile <- cleanFilePath(myPath, paste("ANY_Corrections-", theEbType, "_PriorPlots.png", sep=""))
    myErrorFile <- cleanFilePath(myPath, "error.PNG")
  }
  myTime <- system.time({
    results<-BeaEB(theBeaData@mData, theBeaData@mBatches, par.prior=thePriorFlag,
                   by=theBatchType, covariates=theBatchIdsNotToCorrect,
                   prior.plots=theDoCheckPlotsFlag, priorFile, myErrorFile, theThreads)
  })
  logTiming(theEbType, myPath, myTime)
  correctFile <- NULL
  if ((!is.null(results))&&(TRUE==theWriteToFile))
  {
    correctFile <- cleanFilePath(myPath, "adjusted_matrix.tsv")
    writeDataToFile(results, correctFile)
    results <- correctFile
  }
  else if (TRUE==theWriteToFile)
  {
    # do not delete directory, needed for completion flag
    logInfo("EB_internal - no corrections")
    results <- NULL
  }
  logInfo("EB_internal - completed")
  return(results)
}

EB_withNonParametricPriors<-function(theBeaData,
                                  theBatchIdsNotToCorrect,
                                  theDoCheckPlotsFlag,
                                  theBatchType,
                                  theThreads=1,
                                  thePath=NULL,
                                  theDataVersion=NULL, theTestVersion=NULL,
                                  theWriteToFile=FALSE)
{
  newOutputDir <- addVersionsIfNeeded(thePath, theDataVersion, theTestVersion)
  if (!is.null(newOutputDir))
  {
    checkDirForCreation(newOutputDir)
  }
  EB_internal(theBeaData, theBatchIdsNotToCorrect, theDoCheckPlotsFlag,
              theBatchType, "EBwithNonParametricPriors",
              FALSE, theThreads, newOutputDir, theWriteToFile)
}

EB_withParametricPriors<-function(theBeaData,
                                  theBatchIdsNotToCorrect,
                                  theDoCheckPlotsFlag,
                                  theBatchType,
                                  theThreads=1,
                                  thePath=NULL,
                                  theDataVersion=NULL, theTestVersion=NULL,
                                  theWriteToFile=FALSE)
{
  newOutputDir <- addVersionsIfNeeded(thePath, theDataVersion, theTestVersion)
  if (!is.null(newOutputDir))
  {
    checkDirForCreation(newOutputDir)
  }
  EB_internal(theBeaData, theBatchIdsNotToCorrect, theDoCheckPlotsFlag,
              theBatchType, "EBwithParametricPriors",
              TRUE, theThreads, newOutputDir, theWriteToFile)
}

###############################################################################
###############################################################################

MP_Internal<- function(theBeaData,
                       theBatchType,
                       theMpType,
                       theOverallFlag,
                       thePath=NULL,
                       theDataVersion=NULL, theTestVersion=NULL,
                       theWriteToFile=FALSE)
{
  logInfo("MP_Internal - starting")
  results <- NULL
  myPath <- addVersionsIfNeeded(thePath, theDataVersion, theTestVersion)
  myErrorFile <- NULL
  #setLogging(new("Logging", theFile=cleanFilePath(myPath, "MB.log")))
  if (TRUE==theWriteToFile)
  {
    checkDirForCreation(myPath)
    myErrorFile <- cleanFilePath(myPath, "error.PNG")
  }
  myTime <- system.time({
    results<-BeaMP(theBeaData@mData, theBeaData@mBatches, by=theBatchType, overall=theOverallFlag, theIssuesFile=myErrorFile)
  })
  logTiming(theMpType, myPath, myTime)
  correctFile <- NULL
  if ((!is.null(results))&&(TRUE==theWriteToFile))
  {
    correctFile <- cleanFilePath(myPath, "adjusted_matrix.tsv")
    writeDataToFile(results, correctFile)
    results <- correctFile
  }
  else if (TRUE==theWriteToFile)
  {
    # do not delete directory, needed for completion flag
    logInfo("MP_Internal - no corrections")
    results <- NULL
  }
  logInfo("MP_Internal - completed")
  return(results)
}

MP_Overall<- function(theBeaData,
                      thePath=NULL,
                      theDataVersion=NULL, theTestVersion=NULL,
                      theWriteToFile=FALSE)
{
  MP_Internal(theBeaData, "", "MPOverall", TRUE,
              thePath, theDataVersion, theTestVersion, theWriteToFile)
}

MP_ByBatch<- function(theBeaData, theBatchType,
                      thePath=NULL,
                      theDataVersion=NULL, theTestVersion=NULL,
                      theWriteToFile=FALSE)
{
  MP_Internal(theBeaData, theBatchType, "MPByBatch", FALSE,
              thePath, theDataVersion, theTestVersion, theWriteToFile)
}

###############################################################################
###############################################################################

AN_Internal<- function(theBeaData,
                       theBatchType,
                       theAnType,
                       theDoAdjustFlag,
                       thePath=NULL,
                       theDataVersion=NULL, theTestVersion=NULL,
                       theWriteToFile=FALSE)
{
  logInfo("AN_Internal - starting")
  results <- NULL
  myPath <- addVersionsIfNeeded(thePath, theDataVersion, theTestVersion)
  myErrorFile <- NULL
  #setLogging(new("Logging", theFile=cleanFilePath(myPath, "AN.log")))
  if (TRUE==theWriteToFile)
  {
    checkDirForCreation(myPath)
    myErrorFile <- cleanFilePath(myPath, "error.PNG")
  }
  myTime <- system.time({
    results<-BeaAN(theBeaData@mData, theBeaData@mBatches, by=theBatchType, var.adj=theDoAdjustFlag, theIssuesFile=myErrorFile)
  })
  logTiming(theAnType, myPath, myTime)
  correctFile <- NULL
  if ((!is.null(results))&&(TRUE==theWriteToFile))
  {
    correctFile <- cleanFilePath(myPath, "adjusted_matrix.tsv")
    writeDataToFile(results, correctFile)
    results <- correctFile
  }
  else if (TRUE==theWriteToFile)
  {
    # do not delete directory, needed for completion flag
    logInfo("AN_Internal - no corrections")
    results <- NULL
  }
  logInfo("AN_Internal - completed")
  return(results)
}

AN_Adjusted<-function(theBeaData, theBatchType,
                      thePath=NULL,
                      theDataVersion=NULL, theTestVersion=NULL,
                      theWriteToFile=FALSE)
{
  AN_Internal(theBeaData,
              theBatchType,
              "ANAdjusted",
              TRUE,
              thePath=thePath,
              theDataVersion=theDataVersion,
              theTestVersion=theTestVersion,
              theWriteToFile=theWriteToFile)
}

AN_Unadjusted<-function(theBeaData, theBatchType,
                      thePath=NULL,
                      theDataVersion=NULL, theTestVersion=NULL,
                      theWriteToFile=FALSE)
{
  AN_Internal(theBeaData,
              theBatchType,
              "ANUnadjusted",
              FALSE,
              thePath=thePath,
              theDataVersion=theDataVersion,
              theTestVersion=theTestVersion,
              theWriteToFile=theWriteToFile)
}

###############################################################################
###############################################################################
