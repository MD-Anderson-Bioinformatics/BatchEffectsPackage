# MBatch Copyright (c) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021 University of Texas MD Anderson Cancer Center
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
# MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>

#############################################################################
### UMAP
#############################################################################

UMAP_Structures <- function(theData, theTitle, theOutputPath, theBatchTypeAndValuePairsToRemove=list(), theBatchTypeAndValuePairsToKeep=list())
{
  theData <- as.numericWithIssues(theData)
  theData <- mbatchFilterData(theData, theBatchTypeAndValuePairsToRemove=theBatchTypeAndValuePairsToRemove,
                              theBatchTypeAndValuePairsToKeep=theBatchTypeAndValuePairsToKeep)
  theTitle <- breakIntoTitle(theTitle)
  createBatchEffectsOutput_umap(theData@mData, theData@mBatches, theTitle, theUmapOutputDir=theOutputPath)
}

#############################################################################
### Hierarchical Clustering
#############################################################################

HierarchicalClustering_Structures <- function(theData, theTitle, theOutputPath, theBatchTypeAndValuePairsToRemove=list(), theBatchTypeAndValuePairsToKeep=list())
{
  theData <- as.numericWithIssues(theData)
  theData <- mbatchFilterData(theData, theBatchTypeAndValuePairsToRemove=theBatchTypeAndValuePairsToRemove,
                              theBatchTypeAndValuePairsToKeep=theBatchTypeAndValuePairsToKeep)
  theTitle <- breakIntoTitle(theTitle)
  createBatchEffectsOutput_hierclust(theData@mData, theData@mBatches, theTitle, "", theHierClustOutputDir=theOutputPath)
}

#############################################################################
### Supervised Clustering
#############################################################################

SupervisedClustering_Batches_Structures <- function(theData, theTitle, theOutputPath,
																										theBatchTypeAndValuePairsToRemove=list(),
																										theBatchTypeAndValuePairsToKeep=list())
{
	theData <- as.numericWithIssues(theData)
	theData <- mbatchFilterData(theData, theBatchTypeAndValuePairsToRemove=theBatchTypeAndValuePairsToRemove,
	                            theBatchTypeAndValuePairsToKeep=theBatchTypeAndValuePairsToKeep)
	theTitle <- breakIntoTitle(theTitle)
	createBatchEffectsOutput_SupervisedClustering_batches(theData@mData, theData@mBatches,
	                                                      theTitle=theTitle,
	                                                      theOutputPath=theOutputPath)
}

SupervisedClustering_Pairs_Structures <- function(theData, theTitle, theOutputPath, theListOfBatchPairs,
																									theBatchTypeAndValuePairsToRemove=list(), theBatchTypeAndValuePairsToKeep=list())
{
  theData <- as.numericWithIssues(theData)
  theData <- mbatchFilterData(theData, theBatchTypeAndValuePairsToRemove=theBatchTypeAndValuePairsToRemove,
                              theBatchTypeAndValuePairsToKeep=theBatchTypeAndValuePairsToKeep)
  theTitle <- breakIntoTitle(theTitle)
  createBatchEffectsOutput_SupervisedClustering_pairs(theData@mData, theData@mBatches,
                                                      thePairList=theListOfBatchPairs,
                                                      theTitle=theTitle, theOutputPath=theOutputPath)}

#############################################################################
### PCA
#############################################################################

PCA_Regular_Structures <- function(theData, theTitle, theOutputPath, theBatchTypeAndValuePairsToRemove, theBatchTypeAndValuePairsToKeep,
															theIsPcaTrendFunction=function(...) {FALSE}, theDoCentroidsMtoMFlag=TRUE, theDoPlainMtoMFlag=FALSE,
															theDoDSCFlag=TRUE, theDoDscPermsFileFlag=FALSE, theDoSampleLocatorFlag=TRUE,
															theListOfComponentsToPlot=c(1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4),
															theDSCPermutations=2000, theDSCThreads=1, theMinBatchSize=2,
															theJavaParameters="-Xms1200m", theSeed=NULL, theMaxGeneCount=20000)
{
  theData <- as.numericWithIssues(theData)
  theData <- mbatchFilterData(theData, theBatchTypeAndValuePairsToRemove=theBatchTypeAndValuePairsToRemove,
                              theBatchTypeAndValuePairsToKeep=theBatchTypeAndValuePairsToKeep)
  theTitle <- breakIntoTitle(theTitle)
  createBatchEffectsOutput_pca(theData@mData, theData@mBatches, theTitle,
                               theDoPlainMtoMFlag,
                               theDoCentroidsMtoMFlag,
                               theDoDSCFlag,
                               theDoDscPermsFileFlag,
                               theDoSampleLocatorFlag,
                               theIsPcaTrendFunction=theIsPcaTrendFunction,
                               theListOfComponentsToPlot=theListOfComponentsToPlot,
                               theDSCPermutations=theDSCPermutations,
                               theDSCThreads=theDSCThreads,
                               theMinBatchSize=theMinBatchSize,
                               theOutputDir=theOutputPath,
                               theSeed=theSeed,
                               theGeneLimit=theMaxGeneCount,
                               theJavaParameters=theJavaParameters)
  #if(TRUE==theDoDSCFlag)
  #{
  #	buildDSCOverviewFile(theOutputPath, "DSCOverview.tsv")
  #}
}

PCA_DualBatch_Structures <- function(theData, theTitle, theOutputPath, theBatchTypeAndValuePairsToRemove, theBatchTypeAndValuePairsToKeep,
															theListForDoCentroidDualBatchType, theIsPcaTrendFunction=function(...) {FALSE},
															theDoDSCFlag=TRUE, theDoDscPermsFileFlag=FALSE, theDoSampleLocatorFlag=TRUE,
															theListOfComponentsToPlot=c(1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4),
															theDSCPermutations=2000, theDSCThreads=1, theMinBatchSize=2,
															theJavaParameters="-Xms1200m", theSeed=NULL, theMaxGeneCount=20000)
{
  theData <- as.numericWithIssues(theData)
  theData <- mbatchFilterData(theData, theBatchTypeAndValuePairsToRemove=theBatchTypeAndValuePairsToRemove,
                              theBatchTypeAndValuePairsToKeep=theBatchTypeAndValuePairsToKeep)
  theTitle <- breakIntoTitle(theTitle)
  createBatchEffectsOutput_pca_dualBatch(theData@mData, theData@mBatches,
                                         theListForDoCentroidDualBatchType,
                                         theTitle, theDoDSCFlag, theDoDscPermsFileFlag, theDoSampleLocatorFlag,
                                         theIsPcaTrendFunction, theListOfComponentsToPlot,
                                         theDSCPermutations, theDSCThreads, theMinBatchSize,
                                         theOutputPath,
                                         theSeed=theSeed,
                                         theGeneLimit=theMaxGeneCount,
                                         theJavaParameters=theJavaParameters)
  #if(TRUE==theDoDSCFlag)
  #{
  #  buildDSCOverviewFile(theOutputPath, "DSCOverview.tsv")
  #}
}

#############################################################################
### Boxplot
#############################################################################

Boxplot_AllSamplesData_Structures <- function(theData, theTitle, theOutputPath, theBatchTypeAndValuePairsToRemove, theBatchTypeAndValuePairsToKeep,
                                              theMaxGeneCount=20000)
{
  theData <- as.numericWithIssues(theData)
  theData <- mbatchFilterData(theData, theBatchTypeAndValuePairsToRemove=theBatchTypeAndValuePairsToRemove,
                              theBatchTypeAndValuePairsToKeep=theBatchTypeAndValuePairsToKeep)
  theData@mData <- mbatchTrimData(theData@mData, theMaxSize=(theMaxGeneCount*ncol(theData@mData)))
  theTitle <- breakIntoTitle(theTitle)
  createBatchEffectsOutput_BoxPlot_AllSampleData(theData@mData, theData@mBatches,
                                                 theTitle,
                                                 theOutputPath,
                                                 thePngFlag=TRUE)
}

Boxplot_AllSamplesRLE_Structures <- function(theData, theTitle, theOutputPath, theBatchTypeAndValuePairsToRemove, theBatchTypeAndValuePairsToKeep,
                                             theMaxGeneCount=20000)
{
  theData <- as.numericWithIssues(theData)
  theData <- mbatchFilterData(theData, theBatchTypeAndValuePairsToRemove=theBatchTypeAndValuePairsToRemove,
                              theBatchTypeAndValuePairsToKeep=theBatchTypeAndValuePairsToKeep)
  theData@mData <- mbatchTrimData(theData@mData, theMaxSize=(theMaxGeneCount*ncol(theData@mData)))
  theTitle <- breakIntoTitle(theTitle)
  createBatchEffectsOutput_BoxPlot_AllSampleRLE(theData@mData, theData@mBatches,
                                                 theTitle,
                                                 theOutputPath,
                                                 thePngFlag=TRUE)
}

Boxplot_Group_Structures <- function(theData, theTitle, theOutputPath,
                                     theBatchTypeAndValuePairsToRemove, theBatchTypeAndValuePairsToKeep,
                                     theListOfGroupBoxFunction, theListOfGroupBoxLabels, theMaxGeneCount=20000)
{
  theData <- as.numericWithIssues(theData)
  theData <- mbatchFilterData(theData, theBatchTypeAndValuePairsToRemove=theBatchTypeAndValuePairsToRemove,
                              theBatchTypeAndValuePairsToKeep=theBatchTypeAndValuePairsToKeep)
  theData@mData <- mbatchTrimData(theData@mData, theMaxSize=(theMaxGeneCount*ncol(theData@mData)))
  theTitle <- breakIntoTitle(theTitle)
  createBatchEffectsOutput_BoxPlot_Group(theData@mData, theData@mBatches,
                                         theTitle, theOutputPath,
                                         theListOfGroupBoxFunction, theListOfGroupBoxLabels,
                                         thePngFlag=TRUE)
}

#############################################################################
### TRINOVA
#############################################################################

TRINOVA_Structures<-function(theData, theTitle, theOutputPath,
                            theBatchTypeAndValuePairsToRemove, theBatchTypeAndValuePairsToKeep,
                            theMaxGeneCount=20000)
{
  theData <- as.numericWithIssues(theData)
  theData <- mbatchFilterData(theData, theBatchTypeAndValuePairsToRemove=theBatchTypeAndValuePairsToRemove,
                              theBatchTypeAndValuePairsToKeep=theBatchTypeAndValuePairsToKeep)
  theData@mData <- mbatchTrimData(theData@mData, theMaxSize=(theMaxGeneCount*ncol(theData@mData)))
  theTitle <- breakIntoTitle(theTitle)
  createBatchEffectsOutput_TRINOVA(theData@mData, theData@mBatches, theTitle, theOutputPath)
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


