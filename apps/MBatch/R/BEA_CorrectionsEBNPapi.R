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


doEBNPlus_internal <- function(theBeaData1Matrix,
                      thePath,
                      theIssue,
                      theType,
                      theBeaData2Matrix,
                      theEBNP_TrimGenesFunction,
                      theEBNP_TrimBarcodesFunction,
                      theEBNP_RemoveRowDuplicatesFunction,
                      theEBNP_RemoveColDuplicatesFunction,
                      theEBNP_Data1BatchId,
                      theEBNP_Data2BatchId,
                      theEBNP_BatchWithZero,
                      theEBNP_FixDataSet,
                      theEBNP_CorrectForZero,
                      theSeed,
                      theEBNP_ValidationRatio,
                      thePriorPlotPath,
                      theDataVersion,
                      theTestVersion,
                      thePriorPlotFile,
                      theEBNP_MinSampleNum=3,
                      theEBNP_AddData1Rows=FALSE,
                      theEBNP_AddData2Rows=FALSE)
{
  logInfo("doEBNPlus - starting")
  corrections<-NULL
  myTime <- system.time({
    corrections<-BeaEBNPlus(theBeaData1=theBeaData1Matrix,
                            theIssue=theIssue,
                            theBeaData2=theBeaData2Matrix,
                            theEBNP_TrimGenesFunction=theEBNP_TrimGenesFunction,
                            theEBNP_TrimBarcodesFunction=theEBNP_TrimBarcodesFunction,
                            theEBNP_RemoveRowDuplicatesFunction=theEBNP_RemoveRowDuplicatesFunction,
                            theEBNP_RemoveColDuplicatesFunction=theEBNP_RemoveColDuplicatesFunction,
                            theEBNP_Data1BatchId=theEBNP_Data1BatchId,
                            theEBNP_Data2BatchId=theEBNP_Data2BatchId,
                            theEBNP_BatchWithZero=theEBNP_BatchWithZero,
                            theEBNP_FixDataSet=theEBNP_FixDataSet,
                            theEBNP_CorrectForZero=theEBNP_CorrectForZero,
                            theSeed=theSeed,
                            theEBNP_ValidationRatio=theEBNP_ValidationRatio,
                            theEBNP_TestRatio=0,
                            thePriorPlotPath=thePriorPlotPath,
                            theDataVersion=theDataVersion,
                            theTestVersion=theTestVersion,
                            thePriorPlotFile=thePriorPlotFile,
                            theEBNP_MinSampleNum=theEBNP_MinSampleNum,
                            theEBNP_AddData1Rows=theEBNP_AddData1Rows,
                            theEBNP_AddData2Rows=theEBNP_AddData2Rows)[["CorrectedResults"]]
  })
  logTiming(theType, thePath, myTime)
  correctFile <- NULL
  if (!is.null(corrections))
  {
    checkDirForCreation(thePath)
    correctFile <- cleanFilePath(thePath, "adjusted_matrix.tsv")
    writeDataToFile(corrections, correctFile)
  }
  else
  {
    # do not delete directory, needed for completion flag
    logInfo("doEBNPlus - no corrections")
    correctFile <- NULL
  }
  logInfo("doEBNPlus - completed")
  return(correctFile)
}

EBNPlus_Correction_Files <- function(theDataFile1, theDataFile2, theOutputDir,
                                     theDataVersion, theTestVersion,
                                     theBatchId1, theBatchId2, theSeed=NULL,
  																	 theEBNP_PriorPlotsFlag=FALSE,
  																	 theEBNP_MinSampleNum=3,
  																	 theEBNP_AddData1Rows=FALSE,
  																	 theEBNP_AddData2Rows=FALSE)
{
  # do not update theOutputDir, as update is done via internal function call
	#############################################################################
	#############################################################################
	#### SETTINGS: working, log directory
	setwd(theOutputDir)
	ebnp_PriorPlotsFile <- "ENBPlus_PriorPlots.PNG"
	logInfo("ebnp_PriorPlotsFile=", ebnp_PriorPlotsFile)
	logInfo(getMBatchVersion())
	#############################################################################
	# mbatchLoadFiles
	#############################################################################
	#### SETTINGS: input file names (data, batch and covariates) if needed
	logInfo("BeaEBNplusFiles -- mbatchLoadFiles data1")
	myData1 <- readAsGenericMatrix(theDataFile1)
	logInfo("BeaEBNplusFiles -- mbatchLoadFiles data2")
	myData2 <- readAsGenericMatrix(theDataFile2)
	correctedFilenames <- list()
	currentNames <- names(correctedFilenames)
	correctedFilenames[[length(correctedFilenames)+1]] <- doEBNPlus_internal(theBeaData1Matrix=myData1,
																																	thePath=theOutputDir,
																																	theIssue=cleanFilePath(theOutputDir, "error.txt"),
																																	theType="EBNPlus",
																																	theBeaData2Matrix=myData2,
																																	theEBNP_TrimGenesFunction=trimGenes,
																																	theEBNP_TrimBarcodesFunction=keepBarcodes,
																																	theEBNP_RemoveRowDuplicatesFunction=removeDuplicatesFromRows,
																																	theEBNP_RemoveColDuplicatesFunction=removeDuplicatesFromColumns,
																																	theEBNP_Data1BatchId=theBatchId1,
																																	theEBNP_Data2BatchId=theBatchId2,
																																	theEBNP_BatchWithZero="1",
																																	theEBNP_FixDataSet=as.numeric(NA),
																																	theEBNP_CorrectForZero=TRUE,
																																	theSeed=theSeed,
																																	theEBNP_ValidationRatio=0,
																																	thePriorPlotPath=theOutputDir,
																																	theDataVersion=theDataVersion,
																																	theTestVersion=theTestVersion,
																																	thePriorPlotFile=ebnp_PriorPlotsFile,
																																	theEBNP_MinSampleNum=theEBNP_MinSampleNum,
																																	theEBNP_AddData1Rows=theEBNP_AddData1Rows,
																																	theEBNP_AddData2Rows=theEBNP_AddData2Rows)
	names(correctedFilenames) <- c(currentNames, "EBNPlus")
	logInfo(correctedFilenames)
	correctedFilenames
}

EBNPlus_CombineBatches <- function(theBeaBatches1,
															theBeaBatches2,
															theEBNP_Data1BatchId,
															theEBNP_Data2BatchId,
															theBarcodeTrimFunction=NULL,
															theSep=".")
{
	logDebug("starting BeaEBNPlusBatches")
	barcodes1 <- as.vector(unlist(theBeaBatches1["Sample"]))
	barcodes2 <- as.vector(unlist(theBeaBatches2["Sample"]))
	if (FALSE==is.null(theBarcodeTrimFunction))
	{
		barcodes1 <- theBarcodeTrimFunction(barcodes1)
		barcodes2 <- theBarcodeTrimFunction(barcodes2)
	}
	barcodes1 <- paste(barcodes1, theEBNP_Data1BatchId, sep=theSep)
	barcodes2 <- paste(barcodes2, theEBNP_Data2BatchId, sep=theSep)

	theBeaBatches1["Sample"] <- barcodes1
	theBeaBatches2["Sample"] <- barcodes2
	newBatches <- rbind(theBeaBatches1, theBeaBatches2)
	newBatches["EBNPlus"] <- c(rep(theEBNP_Data1BatchId, length(barcodes1)), rep(theEBNP_Data2BatchId, length(barcodes2)))
	newBatches
}


#############################################################################
### ENBPlus Correction
#############################################################################

EBNPlus_Correction_Structures <- function(theDataMatrix1, theDataMatrix2, theBatchId1, theBatchId2,
																					theEBNP_BatchWithZero, theEBNP_FixDataSet, theEBNP_CorrectForZero,
																					theEBNP_ParametricPriorsFlag,
																					theSeed=NULL,
																					theOutputDir=NULL,
																					theDataVersion=NULL,
																					theTestVersion=NULL,
																					thePriorFile=NULL,
																					theEBNP_MinSampleNum=3,
																					theEBNP_AddData1Rows=FALSE,
																					theEBNP_AddData2Rows=FALSE)
{
  logDebug("EBNPlus_Correction_Structures theEBNP_AddData1Rows=", theEBNP_AddData1Rows)
  logDebug("EBNPlus_Correction_Structures theEBNP_AddData2Rows=", theEBNP_AddData2Rows)
  BeaEBNPlus(theBeaData1=theDataMatrix1,
						 theIssue=NULL,
						 theBeaData2=theDataMatrix2,
						 theEBNP_TrimGenesFunction=NULL,
						 theEBNP_TrimBarcodesFunction=NULL,
						 theEBNP_RemoveRowDuplicatesFunction=removeDuplicatesFromRows,
						 theEBNP_RemoveColDuplicatesFunction=removeDuplicatesFromColumns,
						 theEBNP_Data1BatchId=theBatchId1,
						 theEBNP_Data2BatchId=theBatchId2,
						 theEBNP_BatchWithZero=theEBNP_BatchWithZero,
						 theEBNP_FixDataSet=theEBNP_FixDataSet,
						 theEBNP_CorrectForZero=theEBNP_CorrectForZero,
						 theSeed=theSeed,
						 theEBNP_ValidationRatio=0,
						 theEBNP_TestRatio=0,
						 theEBNP_ParametricPriorsFlag=theEBNP_ParametricPriorsFlag,
						 thePriorPlotPath=theOutputDir,
						 theDataVersion=theDataVersion,
						 theTestVersion=theTestVersion,
						 thePriorPlotFile=thePriorFile,
						 theEBNP_MinSampleNum=theEBNP_MinSampleNum,
						 theEBNP_AddData1Rows=theEBNP_AddData1Rows,
						 theEBNP_AddData2Rows=theEBNP_AddData2Rows)[["CorrectedResults"]]

}

EBNPlus_TrainAndValidateReplicates_Structures <- function(theDataMatrix1, theDataMatrix2, theBatchId1, theBatchId2,
																													theEBNP_BatchWithZero, theEBNP_FixDataSet, theEBNP_CorrectForZero,
																													theEBNP_ParametricPriorsFlag,
																													theEBNP_ValidationRatio,
																													theEBNP_TestRatio=0,
																													theSeed=NULL,
																													thePriorPlotPath=NULL,
																													theDataVersion=NULL,
																													theTestVersion=NULL,
																													thePriorPlotFile=NULL,
																													theEBNP_MinSampleNum=3,
																													theEBNP_AddData1Rows=FALSE,
																													theEBNP_AddData2Rows=FALSE,
																													theTestSeed=NULL)
{
	# [["TrainingSet1"]] [["TrainingSet2"]] [["TrainingResults"]]
	# [["ValidationSet1"]] [["ValidationSet2"]] [["ValidationResults"]]
	# [["CorrectedResults"]]
	foo <- BeaEBNPlus(theBeaData1=theDataMatrix1,
						 theIssue=NULL,
						 theBeaData2=theDataMatrix2,
						 theEBNP_TrimGenesFunction=NULL,
						 theEBNP_TrimBarcodesFunction=NULL,
						 theEBNP_RemoveRowDuplicatesFunction=removeDuplicatesFromRows,
						 theEBNP_RemoveColDuplicatesFunction=removeDuplicatesFromColumns,
						 theEBNP_Data1BatchId=theBatchId1,
						 theEBNP_Data2BatchId=theBatchId2,
						 theEBNP_BatchWithZero=theEBNP_BatchWithZero,
						 theEBNP_FixDataSet=theEBNP_FixDataSet,
						 theEBNP_CorrectForZero=theEBNP_CorrectForZero,
						 theSeed=theSeed,
						 theEBNP_ValidationRatio=theEBNP_ValidationRatio,
						 theEBNP_TestRatio=theEBNP_TestRatio,
						 theEBNP_ParametricPriorsFlag=theEBNP_ParametricPriorsFlag,
						 thePriorPlotPath=thePriorPlotPath,
						 theDataVersion=theDataVersion,
						 theTestVersion=theTestVersion,
						 thePriorPlotFile=thePriorPlotFile,
						 theEBNP_MinSampleNum=theEBNP_MinSampleNum,
						 theEBNP_AddData1Rows=theEBNP_AddData1Rows,
						 theEBNP_AddData2Rows=theEBNP_AddData2Rows,
						 theTestSeed=theTestSeed)
	foo
}

EBNPlus_TrainAndValidateFromVector_Structures <- function(theDataMatrix1, theDataMatrix2,
                                                          theBatchId1, theBatchId2,
																													theEBNP_PsuedoReplicates1Train, theEBNP_PsuedoReplicates2Train,
																													theEBNP_PsuedoReplicates1Validation, theEBNP_PsuedoReplicates2Validation,
																													theEBNP_BatchWithZero,
																													theEBNP_FixDataSet,
																													theEBNP_CorrectForZero,
																													theEBNP_ParametricPriorsFlag,
																													theEBNP_TestRatio=0,
																													theSeed=NULL,
																													theTestSeed=NULL,
																													thePriorPlotPath=NULL,
																													theDataVersion=NULL,
																													theTestVersion=NULL,
																													thePriorPlotFile=NULL,
																													theEBNP_MinSampleNum=3)
{
  logDebug("EBNPlus_TrainAndValidateFromVector_Structures - start")

  myPriorPlotsFile <- NULL
  if (!is.null(thePriorPlotFile))
  {
    priorPlotPath <- addVersionsIfNeeded(thePriorPlotPath, theDataVersion, theTestVersion)
    checkDirForCreation(priorPlotPath)
    myPriorPlotsFile <- file.path(priorPlotPath, thePriorPlotFile)
  }
  # changes here may need to be reflected in both BEA_CorrectionsEBNP.R (EBNPlus) and BEA_EBNPapi.R (EBNPlus_TrainAndValidateFromVector_Structures)
  # if training is null, validation should also be null (training=NULL means use all for training)
  stopifnotWithLogging("if training is null, both validation replicate vectors should also be null (training=NULL means use all for training)",
                       !(
                         (
                          (is.null(theEBNP_PsuedoReplicates1Train))||(is.null(theEBNP_PsuedoReplicates2Train))
                         )
                         &&
                         (
                          (!is.null(theEBNP_PsuedoReplicates1Validation))||(!is.null(theEBNP_PsuedoReplicates2Validation))
                         )
                        )
                       )
	theDataMatrix1 <- removeDuplicatesFromRows(theDataMatrix1)
	theDataMatrix2 <- removeDuplicatesFromRows(theDataMatrix2)
	theDataMatrix1 <- removeDuplicatesFromColumns(theDataMatrix1)
	theDataMatrix2 <- removeDuplicatesFromColumns(theDataMatrix2)
	# if training vectors are null, then use current colnames
	if (is.null(theEBNP_PsuedoReplicates1Train))
	{
	  logWarn("theEBNP_PsuedoReplicates1Train is null, so using all samples for training set")
	  theEBNP_PsuedoReplicates1Train <- colnames(theDataMatrix1)
	}
	if (is.null(theEBNP_PsuedoReplicates2Train))
	{
	  logWarn("theEBNP_PsuedoReplicates2Train is null, so using all samples for training set")
	  theEBNP_PsuedoReplicates2Train <- colnames(theDataMatrix2)
	}
	logDebug("EBNPlus myPriorPlotsFile=", myPriorPlotsFile)
	resultList <- list()
	logDebug(paste("dim(theDataMatrix1)=", paste(dim(theDataMatrix1), sep = " ", collapse = " "), sep = ""))
	printMatrix(theDataMatrix1)
	printElements(rownames(theDataMatrix1)[1:100])
	logDebug(paste("dim(theDataMatrix2)=", paste(dim(theDataMatrix2), sep = " ", collapse = " "), sep = ""))
	printMatrix(theDataMatrix2)
	printElements(rownames(theDataMatrix2)[1:100])

	logDebug("remove unknown genes, that start with ?")
	theDataMatrix1 <- theDataMatrix1[!grepl("^[?]+", rownames(theDataMatrix1)), ]
	theDataMatrix2 <- theDataMatrix2[!grepl("^[?]+", rownames(theDataMatrix2)), ]

	logDebug(paste("dim(theDataMatrix1)=", paste(dim(theDataMatrix1), sep = " ", collapse = " "), sep = ""))
	printMatrix(theDataMatrix1)
	printElements(rownames(theDataMatrix1)[1:100])
	logDebug(paste("dim(theDataMatrix2)=", paste(dim(theDataMatrix2), sep = " ", collapse = " "), sep = ""))
	printMatrix(theDataMatrix2)
	printElements(rownames(theDataMatrix2)[1:100])

	logDebug("make EBNplus")
	ebObj  <- new("EBNplus",
								theDataMatrix1,
								theDataMatrix2,
								NULL,
								NULL,
								NULL,
								NULL)
	logDebug("after EBNplus")
	logDebug(paste("dim(ebObj@mData1)=", paste(dim(ebObj@mData1), sep = " ", collapse = " "), sep = ""))
	printMatrix(ebObj@mData1)
	printElements(rownames(ebObj@mData1)[1:100])
	logDebug(paste("dim(ebObj@mData2)=", paste(dim(ebObj@mData2), sep = " ", collapse = " "), sep = ""))
	printMatrix(ebObj@mData2)
	printElements(rownames(ebObj@mData2)[1:100])
	EBNPlus_CheckData_Structures(ebObj@mData1, ebObj@mData2, theEBNP_PsuedoReplicates1Train, theEBNP_PsuedoReplicates2Train)

	ebObj@DF1batch <- theBatchId1
	ebObj@DF2batch  <- theBatchId2
	ebObj@batchWith0 <- theEBNP_BatchWithZero
	# here fix the HiSeq data. I didn’t fix it because in theory it may get better results in small replicates number. The sample code is for our paper.
	#ebObj@fixSet <- 1
	ebObj@fixSet <- theEBNP_FixDataSet
	ebObj@correct0 <- theEBNP_CorrectForZero
	logDebug("getBiComOrder")
	ebObj  <- getBiComOrder(ebObj)

	logDebug(paste("dim(ebObj@mat1Com)=",paste(dim(ebObj@mat1Com), sep = " ", collapse = " "),sep = ""))
	printMatrix(ebObj@mat1Com)
	printElements(rownames(ebObj@mat1Com)[1:100])
	logDebug(paste("dim(ebObj@mat2Com)=",paste(dim(ebObj@mat2Com), sep = " ", collapse = " "),sep = ""))
	printMatrix(ebObj@mat2Com)
	printElements(rownames(ebObj@mat2Com)[1:100])

	logDebug("get Validation and Training sets from vectors")
	ebObj@mat1Train <- theDataMatrix1[,colnames(theDataMatrix1) %in% theEBNP_PsuedoReplicates1Train ]
	ebObj@mat2Train <- theDataMatrix2[,colnames(theDataMatrix2) %in% theEBNP_PsuedoReplicates2Train ]

	if (length(theEBNP_PsuedoReplicates1Validation)>0)
	{
		ebObj@mat1Validation <- theDataMatrix1[,colnames(theDataMatrix1) %in% theEBNP_PsuedoReplicates1Validation ]
		ebObj@mat2Validation <- theDataMatrix2[,colnames(theDataMatrix2) %in% theEBNP_PsuedoReplicates2Validation ]
	}

	logDebug(paste("dim(ebObj@mat1Validation)=",paste(dim(ebObj@mat1Validation), sep = " ", collapse = " "),sep = ""))
	if (sum(dim(ebObj@mat1Validation)) > 0)
	{
		printMatrix(ebObj@mat1Validation)
	}
	logDebug(paste("dim(ebObj@mat2Validation)=",paste(dim(ebObj@mat2Validation), sep = " ", collapse = " "),sep = ""))
	if (sum(dim(ebObj@mat2Validation)) > 0)
	{
		printMatrix(ebObj@mat2Validation)
	}
	logDebug(paste("dim(ebObj@mat1Train)=",paste(dim(ebObj@mat1Train), sep = " ", collapse = " "),sep = ""))
	if (sum(dim(ebObj@mat1Train)) > 0)
	{
		printMatrix(ebObj@mat1Train)
	}
	logDebug(paste("dim(ebObj@mat2Train)=",paste(dim(ebObj@mat2Train), sep = " ", collapse = " "),sep = ""))
	if (sum(dim(ebObj@mat2Train)) > 0)
	{
		printMatrix(ebObj@mat2Train)
	}
	if (theEBNP_TestRatio!=0)
	{
	  logDebug("before Test Set")
	  # a second seed is used that applies only to the test set creation. (Make sure that TEST is the last thing done)
	  # This allows the user to keep the same test set each time or change it each time, depending on their needs
	  # Machine learning, for example, needs different test sets each time, but you want the validation set to be constant
	  #
	  # take mat1TRain and mat2Train and split that into a new train1 and train2 and a test1 and test2
	  # this calls getTestSet
	  ebObj <- getTestSet(ebObj, testRatio=theEBNP_TestRatio, theTestSeed=theTestSeed )
	  logDebug("after Test Set")
	}

	# set results
	resultList[["TrainingSet1"]] <- sortMatrix(ebObj@mat1Train)
	resultList[["TrainingSet2"]] <- sortMatrix(ebObj@mat2Train)
	resultList[["ValidationSet1"]] <- sortMatrix(ebObj@mat1Validation)
	resultList[["ValidationSet2"]] <- sortMatrix(ebObj@mat2Validation)
	# test set output
	resultList[["TestSet1"]] <- sortMatrix(ebObj@mat1Test)
	resultList[["TestSet2"]] <- sortMatrix(ebObj@mat2Test)

	logDebug("train")
	# does this need seed to be reset?
	testFlag <- FALSE
	if (theEBNP_TestRatio!=0)
	{
	  testFlag <- TRUE
	}
	objafterEB <- train(ebObj, par.prior = theEBNP_ParametricPriorsFlag, test = testFlag,
	                    theEBNP_PriorPlotsFile=myPriorPlotsFile, minSampleNum=theEBNP_MinSampleNum)
	logDebug("after train")

	logDebug("EBadj")
	wholeEB <- NULL
	trainEB <- NULL
	testEB <- NULL
	validationEB <- NULL
	if (0==length(theEBNP_PsuedoReplicates1Validation))
	{
		logDebug("EBadj whole 2")
		objafterEB  <- EBadj(objafterEB , "whole")
		logDebug("objafterEB@wholeEB")
		wholeEB <- objafterEB@wholeEB
	}
	else
	{
		logDebug("EBadj train")
		objafterEB  <- EBadj(objafterEB , "train")
		logDebug("objafterEB@trainEB")
		trainEB <- objafterEB@trainEB
		logDebug("EBadj validation")
		objafterEB  <- EBadj(objafterEB , "validation")
		logDebug("objafterEB@validationEB")
		validationEB <- objafterEB@validationEB
		#####################################################
		if (theEBNP_TestRatio!=0)
		{
		  logDebug("EBadj train for test results")
		  testResults  <- EBadj(objafterEB , "test")
		  logDebug("testResults@testEB")
		  testEB <- testResults@testEB
		}
		#####################################################
	}

	resultList[["CorrectedResults"]] <- sortMatrix(wholeEB)
	resultList[["TrainingResults"]] <- sortMatrix(trainEB)
	resultList[["ValidationResults"]] <- sortMatrix(validationEB)
	resultList[["TestResults"]] <- sortMatrix(testEB)

	# [["TrainingSet1"]] [["TrainingSet2"]] [["TrainingResults"]]
	# [["ValidationSet1"]] [["ValidationSet2"]] [["ValidationResults"]]
	# [["CorrectedResults"]]

	return(resultList)
}

EBNPlus_CheckData_Structures <- function(theDataMatrix1, theDataMatrix2, theDataReplicates1=NULL, theDataReplicates2=NULL)
{
	#matrix
	stopifnotWithLogging("dataset 1 should be a matrix", is.matrix(theDataMatrix1))
	stopifnotWithLogging("dataset 2 should be a matrix", is.matrix(theDataMatrix2))
	#duplicates (sample ids) columns
	mycolnames1 <- colnames(theDataMatrix1)
	stopifnotWithLogging("dataset 1 should have colnames", length(mycolnames1)>1)
	mycolnames2 <- colnames(theDataMatrix2)
	stopifnotWithLogging("dataset 2 should have colnames", length(mycolnames2)>1)
	#gene names matching (x%?) rows
	myrownames1 <- rownames(theDataMatrix1)
	stopifnotWithLogging("dataset 1 should have rownames", length(myrownames1)>1)
	myrownames2 <- rownames(theDataMatrix2)
	stopifnotWithLogging("dataset 2 should have rownames", length(myrownames2)>1)
	stopifnotWithLogging("dataset rownames should have overlap", intersect(myrownames1, myrownames2)>1)
	#replicates or valid lists of "replicates"
	if ((is.null(theDataReplicates1)) && (is.null(theDataReplicates2)))
	{
		# handle replicates by column names
		stopifnotWithLogging("dataset colnames should have more than 3 replicates", length(intersect(mycolnames1, mycolnames2))>3)
		warnifnotWithLogging("dataset colnames should have more than 16 replicates", length(intersect(mycolnames1, mycolnames2))>16)
	}
	else
	{
		# handle replicates from provded lists
		stopifnotWithLogging("replicate replacements 1 should match colnames", length(intersect(theDataReplicates1, mycolnames1))>0)
		stopifnotWithLogging("replicate replacements 2 should match colnames", length(intersect(theDataReplicates2, mycolnames2))>0)
		stopifnotWithLogging("replicate replacements should have more than 3 replicates in set 1", length(intersect(theDataReplicates1, mycolnames1))>3)
		stopifnotWithLogging("replicate replacements should have more than 3 replicates in set 2", length(intersect(theDataReplicates2, mycolnames2))>3)
		warnifnotWithLogging("replicate replacements should have more than 16 replicates in set 1", length(intersect(theDataReplicates1, mycolnames1))>16)
		warnifnotWithLogging("replicate replacements should have more than 16 replicates in set 2", length(intersect(theDataReplicates2, mycolnames2))>16)
	}
	#all data numeric
	stopifnotWithLogging("dataset 1 should be numeric", is.numeric(theDataMatrix1))
	stopifnotWithLogging("dataset 2 should be numeric", is.numeric(theDataMatrix2))
}
