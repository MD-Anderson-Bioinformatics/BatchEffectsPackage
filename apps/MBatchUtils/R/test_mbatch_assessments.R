# MBatchUtils Copyright (c) 2011-2024 University of Texas MD Anderson Cancer Center
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
# MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>

loadDataManually <- function(theGeneFile, theBatchFile)
{
	# this reads files into a matrix (for the gene data) and a data frame (for the batch information)
	# as an example of how to get the MBatch object from mbatchLoadStructures, which is used as an argument to the other MBatch functions
	geneMatrix <- readAsGenericMatrix(theGeneFile)
	batchDataframe <- readAsGenericDataframe(theBatchFile)
	mbatchLoadStructures(geneMatrix, batchDataframe)
}

################################################################################
################################################################################
################################################################################
################################################################################

callMBatch_UMAP_Structures <- function(theOutputDir, theDataVersion, theTestVersion,
                                       theDataObject, theTitle)
{
  # output directory
  outdir <- cleanFilePath(theOutputDir, "UMAP")
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  # here, we call UMAP, passing a title and an output path,
  # and to use the data without removing any type/values
  # and keeping any other defaults
  UMAP_Structures(theDataObject, theTitle,
                  outdir, theDataVersion, theTestVersion)
}

################################################################################
################################################################################
################################################################################
################################################################################

callMBatch_HierarchicalClustering <- function(theOutputDir, theDataVersion, theTestVersion,
                                              theGeneFile, theBatchFile, theTitle, theBoxplotMaxGenes)
{
  # load data
  myData <- loadDataManually(theGeneFile, theBatchFile)
  callMBatch_HierarchicalClustering_Structures(theOutputDir, theDataVersion, theTestVersion, myData, theTitle, theBoxplotMaxGenes)
}

callMBatch_HierarchicalClustering_Structures <- function(theOutputDir, theDataVersion, theTestVersion,
                                                         theDataObject, theTitle, theBoxplotMaxGenes)
{
  # output directory
	outdir <- cleanFilePath(theOutputDir, "HierarchicalClustering")
	dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
	theDataObject <- as.numericWithIssues(theDataObject)
	theDataObject@mData <- mbatchTrimData(theDataObject@mData, theMaxSize = (theBoxplotMaxGenes *
	                                                               ncol(theDataObject@mData)))
	# here, we take all the defaults to hierarchical clustering, passing a title and an output path
	message("outdir=",outdir)
	message("theDataVersion=",theDataVersion)
	message("theTestVersion=",theTestVersion)
	HierarchicalClustering_Structures(theData=theDataObject,
									theTitle=theTitle,
									theOutputDir=outdir,
									theDataVersion=theDataVersion,
									theTestVersion=theTestVersion)
}

################################################################################
################################################################################
################################################################################
################################################################################

callMBatch_SupervisedClustering <- function(theOutputDir, theDataVersion, theTestVersion,
                                            theGeneFile, theBatchFile, theTitle,
                                            theShaidyMapGen, theNgchmWidgetJs, theShaidyMapGenJava,
                                            theNGCHMShaidyMem, theSampleType, theNgchmFeatureMapFile)
{
  myData <- loadDataManually(theGeneFile, theBatchFile)
  callMBatch_SupervisedClustering_Structures(theOutputDir, theDataVersion, theTestVersion, myData, theTitle, theSampleType)
}

callMBatch_SupervisedClustering_Structures <- function(theOutputDir, theDataVersion, theTestVersion,
                                                       theDataObject, theTitle,
                                                       theShaidyMapGen, theNgchmWidgetJs, theShaidyMapGenJava, theNGCHMShaidyMem,
                                                       theNgchmRowType,
                                                       theNgchmColumnType,
                                                       theNgchmFeatureMapFile)
{
  # output directory
	outdir <- cleanFilePath(theOutputDir, "SupervisedClustering")
	dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
	message("callMBatch_SupervisedClustering_Structures")
	# here, we call supervised clustering, passing a title and an output path,
	# telling it to generate a heatmap and to use the data without removing any type/values
	# and keeping any other defaults
	message("trim to same size as hierarchical 2")
	SupervisedClustering_Batches_Structures(theData=theDataObject,
									theTitle=theTitle,
									theOutputDir=outdir,
									theDataVersion=theDataVersion,
									theTestVersion=theTestVersion,
									theBatchTypeAndValuePairsToRemove=NULL,
									theBatchTypeAndValuePairsToKeep=NULL)
	for(myBatchType in colnames(theDataObject@mBatches)[2:length(colnames(theDataObject@mBatches))])
  {
	  warnLevel<-getOption("warn")
	  on.exit(options(warn=warnLevel))
	  options(warn=-1)
    message("title ", theTitle)
    message("myBatchType ", myBatchType)
    message("theShaidyMapGen ", theShaidyMapGen)
    message("theNgchmWidgetJs ", theNgchmWidgetJs)
    message("theShaidyMapGenJava ", theShaidyMapGenJava)
    message("dim(myMBatchData@mBatches) ", dim(theDataObject@mBatches))
    message("trim to same size as hierarchical 3")
    dataDirVersioned <- cleanFilePath(outdir, myBatchType)
    dataDirVersioned <- addVersionsIfNeeded(dataDirVersioned, theDataVersion, theTestVersion)
    matrixFile <- cleanFilePath(dataDirVersioned, paste(myBatchType, "SCMatrix.tsv", sep="_"))
    message("matrixFile ", matrixFile)
    myMatrix <- theDataObject@mData
    colDend <- NULL
    rowDend <- NULL
    if (file.exists(matrixFile))
    {
      # matrix from file is transposed
      myMatrix <- readAsGenericMatrix(matrixFile)
      colDend <- cleanFilePath(dataDirVersioned, paste(myBatchType, "uDend.RData", sep="_"))
      rowDend <- cleanFilePath(dataDirVersioned, paste(myBatchType, "uDend_feature.RData", sep="_"))
    }
    message("dim(myMatrix) ", dim(myMatrix))
    message("colDend ", colDend)
    message("rowDend ", rowDend)
    #c("pearson", "ward.D2"),
    # do not use feature map, since supervised clustering is sample ids in both directions
    buildBatchHeatMapFromHC_Structures(theMatrixData=myMatrix,
                                       theBatchData=theDataObject@mBatches,
                                       theTitle=paste(theTitle, "/", myBatchType, "Supervised Clustering", sep=" "),
                                       theOutputDir=cleanFilePath(outdir, myBatchType),
                                       theDataVersion=theDataVersion,
                                       theTestVersion=theTestVersion,
                                       theOutputFilename=paste(myBatchType, "ngchm.ngchm", sep="_"),
                                       theColDendRDataFile=colDend,
                                       theRowDendRDataFile=NULL,
                                       theNgchmFeatureMapFile=theNgchmFeatureMapFile,
                                       theRowType=theNgchmRowType, theColType=theNgchmColumnType,
                                       theRowCluster=c("pearson", "ward.D2"),
                                       theShaidyMapGen=theShaidyMapGen,
                                       theNgchmWidgetJs=theNgchmWidgetJs,
                                       theShaidyMapGenJava=theShaidyMapGenJava,
                                       theShaidyMapGenArgs=c(paste(c("-Xms", "-Xmx"), theNGCHMShaidyMem, sep=""), "-Djava.awt.headless=true"))
	}
}

################################################################################
################################################################################
################################################################################
################################################################################

callMBatch_SupervisedClusteringPairs <- function(theOutputDir, theGeneFile, theBatchFile, theTitle)
{
  # load data
  myData <- loadDataManually(theGeneFile, theBatchFile)
  callMBatch_SupervisedClusteringPairs_Structures(theOutputDir, myData, theTitle)
}

callMBatch_SupervisedClusteringPairs_Structures <- function(theOutputDir, theDataObject, theTitle)
{
  # output directory
	outdir <- cleanFilePath(theOutputDir, "SupervisedClusteringPairs")
	dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
	# here, we call supervised clustering pairs, passing a title and an output path,
	# telling it to generate a heatmap and tell it we want the BatchId and PlateId data paired and the TSS and ShipDate data paired
	# and to use the data without removing any type/values
	# and keeping any other defaults
	SupervisedClustering_Pairs_Structures(theData=theDataObject,
									theTitle=theTitle,
									theOutputDir=outdir,
									theListOfBatchPairs=c("BatchId", "PlateId", "TSS", "ShipDate"),
									theBatchTypeAndValuePairsToRemove=NULL,
									theBatchTypeAndValuePairsToKeep=NULL)
}

################################################################################
################################################################################
################################################################################
################################################################################

isTrendBatch<-function(theBatchTypeName, theListOfBatchIds)
{
  return(is.element(theBatchTypeName, c("ShipDate")))
}

callMBatch_PCA <- function(theOutputDir, theGeneFile, theBatchFile, theTitle, aTest=TRUE)
{
  # load data
  myData <- loadDataManually(theGeneFile, theBatchFile)
  callMBatch_PCA_Structures(theOutputDir, myData, theTitle, isTrendBatch)
  if(TRUE==aTest)
  {
    buildDSCOverviewFile(cleanFilePath(theOutputDir, "PCA"),
                         cleanFilePath(theOutputDir, "DSC"),
                         "DSCOverview.tsv")
    clearDSCOverviewFiles(cleanFilePath(theOutputDir, "PCA"))
  }
}

callMBatch_PCA_Structures <- function(theOutputDir, theDataVersion, theTestVersion,
                                      theDataObject, theTitle,
                                      theIsPcaTrendFunction,
                                      theDSCPermutations=1000,
                                      theDSCThreads=1,
                                      theMinBatchSize=2,
                                      theSeed=314,
                                      theMaxGeneCount=10000)
{
  # output directory
	outdir <- cleanFilePath(theOutputDir, "PCA")
	dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
	# here, we call PCA, passing a title and an output path,
	# and to use the data without removing any type/values
	# and keeping any other defaults
	PCA_Regular_Structures(theData=theDataObject,
							theTitle=theTitle,
							theOutputDir=outdir,
							theDataVersion=theDataVersion,
							theTestVersion=theTestVersion,
							theBatchTypeAndValuePairsToRemove=NULL,
							theBatchTypeAndValuePairsToKeep=NULL,
							theDoDscPermsFileFlag = TRUE,
							theIsPcaTrendFunction=theIsPcaTrendFunction,
							theDSCPermutations=theDSCPermutations,
							theDSCThreads=theDSCThreads,
							theMinBatchSize=theMinBatchSize,
							theSeed=theSeed,
							theMaxGeneCount=theMaxGeneCount)
}

################################################################################
################################################################################
################################################################################
################################################################################

callMBatch_PCADualBatch <- function(theOutputDir, theGeneFile, theBatchFile, theTitle, aTest=TRUE)
{
  myData <- loadDataManually(theGeneFile, theBatchFile)
  callMBatch_PCADualBatch_Structures(theOutputDir, myData, theTitle, isTrendBatch)
  if(TRUE==aTest)
  {
    buildDSCOverviewFile(cleanFilePath(theOutputDir, "PCADualBatch"),
                         cleanFilePath(theOutputDir, "PCADualBatch"),
                         "DSCOverview.tsv")
    clearDSCOverviewFiles(cleanFilePath(theOutputDir, "PCADualBatch"))
  }
}

callMBatch_PCADualBatch_Structures <- function(theOutputDir, theDataObject, theTitle,
                                               theIsPcaTrendFunction,
                                               theDSCPermutations=1000,
                                               theDSCThreads=1,
                                               theMinBatchSize=2,
                                               theSeed=314,
                                               theMaxGeneCount=10000)
{
  # output directory
	outdir <- cleanFilePath(theOutputDir, "PCADualBatch")
	dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
	# here, we call PCA dual batch, passing a title and an output path,
	# and to use the data without removing any type/values
	# and tell it we want the BatchId and PlateId data paired and the TSS and ShipDate data paired
	# and keeping any other defaults
	PCA_DualBatch_Structures(theData=theDataObject,
							theTitle=theTitle,
							theOutputDir=outdir,
							theBatchTypeAndValuePairsToRemove=NULL,
							theBatchTypeAndValuePairsToKeep=NULL,
							theListForDoCentroidDualBatchType=c("BatchId", "PlateId", "TSS", "ShipDate"),
							theIsPcaTrendFunction=theIsPcaTrendFunction,
							theDSCPermutations=theDSCPermutations,
							theDSCThreads=theDSCThreads,
							theMinBatchSize=theMinBatchSize,
							theSeed=theSeed,
							theMaxGeneCount=theMaxGeneCount)
}

################################################################################
################################################################################
################################################################################
################################################################################

callMBatch_BoxplotAllSamplesData <- function(theOutputDir, theGeneFile, theBatchFile, theTitle, theMaxGeneCount = 10000)
{
  # load data
  myData <- loadDataManually(theGeneFile, theBatchFile)
  callMBatch_BoxplotAllSamplesData_Structures(theOutputDir, myData, theTitle, theMaxGeneCount)
}

callMBatch_BoxplotAllSamplesData_Structures <- function(theOutputDir, theDataVersion, theTestVersion, theDataObject, theTitle, theMaxGeneCount = 10000)
{
  # output directory
	outdir <- cleanFilePath(theOutputDir, "BoxPlot")
	dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
	# here, we call boxplot all samples data, passing a title and an output path,
	# and to use the data without removing any type/values
	# and tell it to use 8G of memory
	# and keeping any other defaults
	Boxplot_AllSamplesData_Structures(theData=theDataObject,
                                    theTitle=theTitle,
                    							  theOutputDir=outdir,
                    							  theDataVersion=theDataVersion,
                    							  theTestVersion=theTestVersion,
                    							  theBatchTypeAndValuePairsToRemove=NULL,
                    							  theBatchTypeAndValuePairsToKeep=NULL,
                    							  theMaxGeneCount = theMaxGeneCount)
}

################################################################################
################################################################################
################################################################################
################################################################################

callMBatch_BoxplotAllSamplesRLE <- function(theOutputDir, theGeneFile, theBatchFile, theTitle, theMaxGeneCount = 10000)
{
  # load data
  myData <- loadDataManually(theGeneFile, theBatchFile)
  callMBatch_BoxplotAllSamplesRLE_Structures(theOutputDir, myData, theTitle)
}

callMBatch_BoxplotAllSamplesRLE_Structures <- function(theOutputDir, theDataVersion, theTestVersion, theDataObject, theTitle, theMaxGeneCount = 10000)
{
	# output directory
  outdir <- cleanFilePath(theOutputDir, "BoxPlot")
	dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
	# here, we call boxplot all samples RLE, passing a title and an output path,
	# and to use the data without removing any type/values
	# and tell it to use 8G of memory
	# and keeping any other defaults
	Boxplot_AllSamplesRLE_Structures(theData=theDataObject,
                                   theTitle=theTitle,
                                   theOutputDir=outdir,
                      						 theDataVersion=theDataVersion,
                      						 theTestVersion=theTestVersion,
                      						 theBatchTypeAndValuePairsToRemove=NULL,
                      						 theBatchTypeAndValuePairsToKeep=NULL,
                                	 theMaxGeneCount=theMaxGeneCount)
}

################################################################################
################################################################################
################################################################################
################################################################################

callMBatch_BoxplotGroup <- function(theOutputDir, theGeneFile, theBatchFile, theTitle, theMaxGeneCount = 10000,
                                    theFunction=c(mean), theFunctionName=c("Mean"))
{
  # load data
  myData <- loadDataManually(theGeneFile, theBatchFile)
  callMBatch_BoxplotGroup_Structures(theOutputDir, myData, theTitle,
                                     theMaxGeneCount,
                                     theFunction, theFunctionName)
}

callMBatch_BoxplotGroup_Structures <- function(theOutputDir, theDataVersion, theTestVersion, theDataObject, theTitle, theMaxGeneCount = 10000,
                                               theFunction=list(mean), theFunctionName=list("Mean"))
{
	# output directory
  outdir <- cleanFilePath(theOutputDir, "BoxPlot")
	dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
	# here, we call boxplot group, passing a title and an output path,
	# and to use the data without removing any type/values
	# and tu use the mean function and the label "mean"
	# and tell it to use 8G of memory
	# and keeping any other defaults
	Boxplot_Group_Structures(theData=theDataObject,
	                                      theTitle=theTitle,
                          							theOutputDir=outdir,
                          							theDataVersion=theDataVersion,
                          							theTestVersion=theTestVersion,
                          							theBatchTypeAndValuePairsToRemove=NULL,
                          							theBatchTypeAndValuePairsToKeep=NULL,
                          							theListOfGroupBoxFunction=theFunction,
                          							theListOfGroupBoxLabels=theFunctionName,
                          							theMaxGeneCount=theMaxGeneCount)
}

################################################################################
################################################################################
################################################################################
################################################################################

callMBatch_VolcanoPlot_Structures <- function(theOutputDir, theDataVersion, theTestVersion,
                                              theDataObject, theTitle, theLogFrameFlag, theMaxGeneCount = 10000)
{
  # output directory
  outdir <- cleanFilePath(theOutputDir, "VolcanoPlot")
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  # here, we call volcano plot passing a title and an output path,
  Volcano_Structures(theData=theDataObject,
                     theTitle=theTitle,
                     theOutputDir=outdir,
                     theLogFrameFlag=theLogFrameFlag,
                     theBatchTypeAndValuePairsToRemove=NULL,
                     theBatchTypeAndValuePairsToKeep=NULL,
                     theDataVersion=theDataVersion,
                     theTestVersion=theTestVersion,
                     theMaxFeatureCount=theMaxGeneCount)
}

################################################################################
################################################################################
################################################################################
################################################################################
