## ---- echo=TRUE----------------------------------------------------------
{
  library(MBatch)
  
  # set the paths
  theGeneFile="/bea_testing/MATRIX_DATA/matrix_data-Tumor.tsv"
  theBatchFile="/bea_testing/MATRIX_DATA/batches-Tumor.tsv"
  theOutputDir="/bea_testing/output/PCA_DualBatch_Structures"
  theRandomSeed=314
  
  # make sure the output dir exists and is empty
  unlink(theOutputDir, recursive=TRUE)
  dir.create(theOutputDir, showWarnings=FALSE, recursive=TRUE)
  # load data
  myData <- mbatchLoadFiles(theGeneFile, theBatchFile)
  myData@mData <- mbatchTrimData(myData@mData, 100000)
  # trend function to handle time data
	isTrendBatch<-function(theBatchTypeName, theListOfBatchIds)
	{
		return(is.element(theBatchTypeName, c("ShipDate")))
	}
	###
  # do two plots PlateId versus TSS and BatchId versus TSS
	PCA_DualBatch_Structures(myData,
		"Example Title for Test Data", theOutputDir, list(), list(),
		theIsPcaTrendFunction=isTrendBatch,
		theListForDoCentroidDualBatchType=c("PlateId", "TSS", "BatchId", "TSS"),
		theDoDSCFlag=TRUE, theDoDscPermsFileFlag=TRUE, theDoSampleLocatorFlag=TRUE,
		theListOfComponentsToPlot=c(1, 2, 1, 3), theDSCPermutations=100,
		theDSCThreads=5, theMinBatchSize=0,
		theJavaParameters="-Xms1000m", theSeed=0, theMaxGeneCount=0)
	print(dir(theOutputDir, recursive=TRUE))
}

