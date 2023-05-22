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

checkIfTestError <- function()
{
  if (isTRUE(getGlobalMBatchErrorTest()))
  {
    stop("setGlobalMBatchErrorTest is set to TRUE, so stop with this message, to test error writing", call.=FALSE)
  }
}


writeBatchDataTsvForBoxplot <- function(theOutputDir, theBatches)
{
  message("writeBatchDataTsvForBoxplot")
  batchDataFile <- cleanFilePath(theOutputDir, "BatchData.tsv")
  message("writeAsGenericDataframe writeBatchDataTsvForBoxplot")
  message(batchDataFile)
  if (!file.exists(batchDataFile))
  {
    message("Writing BatchData.tsv now")
    writeAsGenericDataframe(batchDataFile, theBatches)
  }
}

#############################################################################
### UMAP
#############################################################################

UMAP_Structures <- function(theData, theTitle, theOutputDir, theDataVersion, theTestVersion,
                            theBatchTypeAndValuePairsToRemove=list(), theBatchTypeAndValuePairsToKeep=list(),
                            theDoDSCFlag=FALSE, theDSCPermutations=0,
                            theDSCThreads=0, theDoDscPermsFileFlag=FALSE, theSeed=314)
{
  theData <- as.numericWithIssues(theData)
  theData <- mbatchFilterData(theData, theBatchTypeAndValuePairsToRemove=theBatchTypeAndValuePairsToRemove,
                              theBatchTypeAndValuePairsToKeep=theBatchTypeAndValuePairsToKeep)
  theTitle <- breakIntoTitle(theTitle)
  # theMatrix, theDataframeBatchData, theTitle, theDataVersion, theTestVersion, theUmapOutputDir, theUmapFileBase
  createBatchEffectsOutput_umap(theData@mData, theData@mBatches, theTitle,
                                theDataVersion, theTestVersion,
                                theOutputDir,
                                theDoDSCFlag, theDSCPermutations,
                                theDSCThreads, theDoDscPermsFileFlag, theSeed)
}

#############################################################################
### Hierarchical Clustering
#############################################################################

HierarchicalClustering_Structures <- function(theData, theTitle, theOutputDir, theDataVersion, theTestVersion,
                                              theBatchTypeAndValuePairsToRemove=list(), theBatchTypeAndValuePairsToKeep=list())
{
  theData <- as.numericWithIssues(theData)
  theData <- mbatchFilterData(theData, theBatchTypeAndValuePairsToRemove=theBatchTypeAndValuePairsToRemove,
                              theBatchTypeAndValuePairsToKeep=theBatchTypeAndValuePairsToKeep)
  theTitle <- breakIntoTitle(theTitle)
  dirVector <- createBatchEffectsOutput_hierclust(theData@mData, theData@mBatches, theTitle, theDataVersion, theTestVersion, "", theHierClustOutputDir=theOutputDir)
  logInfo("HierarchicalClustering_Structures dirVector=", dirVector)
  newOutDir <- dirVector[1]
  rdataFileSamples <- dirVector[2]
  rdataFileFeatures <- dirVector[3]
  logInfo("HierarchicalClustering_Structures rdataFileSamples=", rdataFileSamples)
  logInfo("HierarchicalClustering_Structures rdataFileFeatures=", rdataFileFeatures)
  writeBatchDataTsvForBoxplot(newOutDir, theData@mBatches)
  c(rdataFileSamples, rdataFileFeatures)
}

#############################################################################
### Supervised Clustering
#############################################################################

SupervisedClustering_Batches_Structures <- function(theData, theTitle, theOutputDir, theDataVersion, theTestVersion,
                                                    theBatchTypeAndValuePairsToRemove=list(), theBatchTypeAndValuePairsToKeep=list())
{
	theData <- as.numericWithIssues(theData)
	theData <- mbatchFilterData(theData, theBatchTypeAndValuePairsToRemove=theBatchTypeAndValuePairsToRemove,
	                            theBatchTypeAndValuePairsToKeep=theBatchTypeAndValuePairsToKeep)
	theTitle <- breakIntoTitle(theTitle)
	createBatchEffectsOutput_SupervisedClustering_batches(theData@mData, theData@mBatches, theDataVersion, theTestVersion,
	                                                      theTitle, theOutputDir)
}

SupervisedClustering_Pairs_Structures <- function(theData, theTitle, theOutputDir, theListOfBatchPairs, theDataVersion, theTestVersion,
                                                  theBatchTypeAndValuePairsToRemove=list(), theBatchTypeAndValuePairsToKeep=list())
{
  theData <- as.numericWithIssues(theData)
  theData <- mbatchFilterData(theData, theBatchTypeAndValuePairsToRemove=theBatchTypeAndValuePairsToRemove,
                              theBatchTypeAndValuePairsToKeep=theBatchTypeAndValuePairsToKeep)
  theTitle <- breakIntoTitle(theTitle)
  createBatchEffectsOutput_SupervisedClustering_pairs(theData@mData, theData@mBatches,
                                                      thePairList=theListOfBatchPairs, theDataVersion, theTestVersion,
                                                      theTitle, theOutputDir)}

#############################################################################
### PCA
#############################################################################

PCA_Regular_Structures <- function(theData, theTitle, theOutputDir, theBatchTypeAndValuePairsToRemove, theBatchTypeAndValuePairsToKeep,
                                   theDataVersion, theTestVersion,
                                   theIsPcaTrendFunction=function(...) {FALSE}, theDoCentroidsMtoMFlag=TRUE, theDoPlainMtoMFlag=FALSE,
                                   theDoDSCFlag=TRUE, theDoDscPermsFileFlag=FALSE, theDoSampleLocatorFlag=TRUE,
                                   theListOfComponentsToPlot=c(1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4),
                                   theDSCPermutations=2000, theDSCThreads=1, theMinBatchSize=2,
                                   theSeed=NULL, theMaxGeneCount=20000)
{
  theData <- as.numericWithIssues(theData)
  theData <- mbatchFilterData(theData, theBatchTypeAndValuePairsToRemove=theBatchTypeAndValuePairsToRemove,
                              theBatchTypeAndValuePairsToKeep=theBatchTypeAndValuePairsToKeep)
  theTitle <- breakIntoTitle(theTitle)
  newOutDirVector <- createBatchEffectsOutput_pca(theData@mData, theData@mBatches, theTitle,
                               theDoPlainMtoMFlag,
                               theDoCentroidsMtoMFlag,
                               theDoDSCFlag,
                               theDoDscPermsFileFlag,
                               theDoSampleLocatorFlag,
                               theDataVersion=theDataVersion,
                               theTestVersion=theTestVersion,
                               theIsPcaTrendFunction=theIsPcaTrendFunction,
                               theListOfComponentsToPlot=theListOfComponentsToPlot,
                               theDSCPermutations=theDSCPermutations,
                               theDSCThreads=theDSCThreads,
                               theMinBatchSize=theMinBatchSize,
                               theOutputDir=theOutputDir,
                               theSeed=theSeed,
                               theGeneLimit=theMaxGeneCount)
  for (newOutDir in newOutDirVector)
  {
    writeBatchDataTsvForBoxplot(newOutDir, theData@mBatches)
  }
  #if(TRUE==theDoDSCFlag)
  #{
  #	buildDSCOverviewFile(theOutputDir, "DSCOverview.tsv")
  #}
}

PCA_DualBatch_Structures <- function(theData, theTitle,
                              theBatchTypeAndValuePairsToRemove, theBatchTypeAndValuePairsToKeep,
															theListForDoCentroidDualBatchType,
															theOutputDir, theDataVersion, theTestVersion,
															theIsPcaTrendFunction=function(...) {FALSE},
															theDoDSCFlag=TRUE, theDoDscPermsFileFlag=FALSE, theDoSampleLocatorFlag=TRUE,
															theListOfComponentsToPlot=c(1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4),
															theDSCPermutations=2000, theDSCThreads=1, theMinBatchSize=2,
															theSeed=NULL, theMaxGeneCount=20000)
{
  theData <- as.numericWithIssues(theData)
  theData <- mbatchFilterData(theData, theBatchTypeAndValuePairsToRemove=theBatchTypeAndValuePairsToRemove,
                              theBatchTypeAndValuePairsToKeep=theBatchTypeAndValuePairsToKeep)
  theTitle <- breakIntoTitle(theTitle)
  newOutDirVector <- createBatchEffectsOutput_pca_dualBatch(theData@mData, theData@mBatches,
                                         theListForDoCentroidDualBatchType,
                                         theTitle, theDoDSCFlag, theDoDscPermsFileFlag, theDoSampleLocatorFlag,
                                         theIsPcaTrendFunction, theListOfComponentsToPlot,
                                         theDSCPermutations, theDSCThreads, theMinBatchSize,
                                         theOutputDir, theDataVersion, theTestVersion,
                                         theSeed=theSeed,
                                         theGeneLimit=theMaxGeneCount)
  for (newOutDir in newOutDirVector)
  {
    writeBatchDataTsvForBoxplot(newOutDir, theData@mBatches)
  }
  #if(TRUE==theDoDSCFlag)
  #{
  #  buildDSCOverviewFile(theOutputDir, "DSCOverview.tsv")
  #}
}

#############################################################################
### Boxplot
#############################################################################

Boxplot_AllSamplesData_Structures <- function(theData, theTitle, theOutputDir,
                                              theBatchTypeAndValuePairsToRemove, theBatchTypeAndValuePairsToKeep,
                                              theDataVersion, theTestVersion,
                                              theMaxGeneCount=20000)
{
  theData <- as.numericWithIssues(theData)
  theData <- mbatchFilterData(theData, theBatchTypeAndValuePairsToRemove=theBatchTypeAndValuePairsToRemove,
                              theBatchTypeAndValuePairsToKeep=theBatchTypeAndValuePairsToKeep)
  theData@mData <- mbatchTrimData(theData@mData, theMaxSize=(theMaxGeneCount*ncol(theData@mData)))
  theTitle <- breakIntoTitle(theTitle)
  newOutDir <- createBatchEffectsOutput_BoxPlot_AllSampleData(theData@mData, theData@mBatches,
                                                              theTitle, theOutputDir,
                                                              theDataVersion, theTestVersion,
                                                              thePngFlag=TRUE)
  writeBatchDataTsvForBoxplot(newOutDir, theData@mBatches)
}

Boxplot_AllSamplesRLE_Structures <- function(theData, theTitle, theOutputDir,
                                             theBatchTypeAndValuePairsToRemove, theBatchTypeAndValuePairsToKeep,
                                             theDataVersion, theTestVersion,
                                             theMaxGeneCount=20000)
{
  theData <- as.numericWithIssues(theData)
  theData <- mbatchFilterData(theData, theBatchTypeAndValuePairsToRemove=theBatchTypeAndValuePairsToRemove,
                              theBatchTypeAndValuePairsToKeep=theBatchTypeAndValuePairsToKeep)
  theData@mData <- mbatchTrimData(theData@mData, theMaxSize=(theMaxGeneCount*ncol(theData@mData)))
  theTitle <- breakIntoTitle(theTitle)
  newOutDir <- createBatchEffectsOutput_BoxPlot_AllSampleRLE(theData@mData, theData@mBatches,
                                                             theTitle, theOutputDir,
                                                             theDataVersion, theTestVersion,
                                                             thePngFlag=TRUE)
  writeBatchDataTsvForBoxplot(newOutDir, theData@mBatches)
}

Boxplot_Group_Structures <- function(theData, theTitle, theOutputDir,
                                     theBatchTypeAndValuePairsToRemove, theBatchTypeAndValuePairsToKeep,
                                     theListOfGroupBoxFunction, theListOfGroupBoxLabels,
                                     theDataVersion, theTestVersion,
                                     theMaxGeneCount=20000)
{
  theData <- as.numericWithIssues(theData)
  theData <- mbatchFilterData(theData, theBatchTypeAndValuePairsToRemove=theBatchTypeAndValuePairsToRemove,
                              theBatchTypeAndValuePairsToKeep=theBatchTypeAndValuePairsToKeep)
  theData@mData <- mbatchTrimData(theData@mData, theMaxSize=(theMaxGeneCount*ncol(theData@mData)))
  theTitle <- breakIntoTitle(theTitle)
  newOutDirVector <- createBatchEffectsOutput_BoxPlot_Group(theData@mData, theData@mBatches,
                                                      theTitle, theOutputDir,
                                                      theListOfGroupBoxFunction, theListOfGroupBoxLabels,
                                                      theDataVersion, theTestVersion,
                                                      thePngFlag=TRUE)
  for (newOutDir in newOutDirVector)
  {
    writeBatchDataTsvForBoxplot(newOutDir, theData@mBatches)
  }
}

#############################################################################
#############################################################################

as.numericWithIssues<-function(theData)
{
	warnLevel<-getOption("warn")
	on.exit(options(warn=warnLevel))
	options(warn=3) # warnings are errors
	mode(theData@mData) <- "numeric"  #Sets all elements of the data frame to numeric
	return(theData)
}


