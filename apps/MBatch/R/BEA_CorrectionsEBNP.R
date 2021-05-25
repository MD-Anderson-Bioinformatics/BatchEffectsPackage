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

BeaEBNPlus <- function(theBeaData1,
											 theIssue,
											 theBeaData2,
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
											 theEBNP_TestRatio,
											 theEBNP_ParametricPriorsFlag=TRUE,
											 theEBNP_PriorPlotsFile=NULL,
											 theEBNP_MinSampleNum=3,
											 theEBNP_AddData1Rows=FALSE,
											 theEBNP_AddData2Rows=FALSE,
											 theTestSeed=NULL)
{
	logDebug("starting BeaEBNPlus")
	logDebug(getMBatchVersion())
	logDebug("BeaEBNPlus theEBNP_AddData1Rows=", theEBNP_AddData1Rows)
	logDebug("BeaEBNPlus theEBNP_AddData2Rows=", theEBNP_AddData2Rows)
	foo <- NULL
	#tryCatch(
		foo <- EBNPlus(
			theData1=theBeaData1,
			theData2=theBeaData2,
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
			theEBNP_ParametricPriorsFlag=theEBNP_ParametricPriorsFlag,
			theEBNP_ValidationRatio=theEBNP_ValidationRatio,
			theEBNP_TestRatio=theEBNP_TestRatio,
			theEBNP_PriorPlotsFile=theEBNP_PriorPlotsFile,
			theEBNP_MinSampleNum=theEBNP_MinSampleNum,
			theEBNP_AddData1Rows=theEBNP_AddData1Rows,
			theEBNP_AddData2Rows=theEBNP_AddData2Rows,
			theTestSeed=theTestSeed
	)
	#,
	#error=function(e)
	#{
	#	handleIssuesFunction(e, theIssue)
	#})
	logDebug("finishing BeaEBNPlus")
	return(foo)
}

printElements <- function(theElements)
{
	# if (is.list(theElements))
	# {
	# 	listSummary(theElements)
	# }
	# else
	# {
	# 	print(theElements[1:min(5,length(theElements))])
	# }
}

printMatrix <- function(theMatrix)
{
	# print(is.matrix(theMatrix))
	# print(dim(theMatrix))
	# rowMax <- dim(theMatrix)[1]
	# colMax <- dim(theMatrix)[2]
	# rowMax <- min(rowMax, 4)
	# colMax <- min(colMax, 4)
	# print(theMatrix[1:rowMax, 1:colMax])
}

listSummary <- function(theList)
{
	# for(myname in names(theList))
	# {
	# 	message(myname)
	# 	if (is.matrix(theList[[myname]]))
	# 	{
	# 		message(myname)
	# 		if (sum(dim(theList[[myname]]))>0)
	# 		{
	# 			print(theList[[myname]][1:min(4,dim(theList[[myname]])[1]), 1:min(4,dim(theList[[myname]])[2])])
	# 		}
	# 	}
	# 	else if (is.data.frame(theList[[myname]]))
	# 	{
	# 		message(myname)
	# 		if (sum(dim(theList[[myname]]))>0)
	# 		{
	# 			print(theList[[myname]][1:min(5,dim(theList[[myname]])[1]),1:min(4,dim(theList[[myname]])[2])])
	# 		}
	# 	}
	# 	else if (is.list(theList[[myname]]))
	# 	{
	# 		message(myname)
	# 		listSummary(theList[[myname]])
	# 	}
	# 	else if (is.vector(theList[[myname]]))
	# 	{
	# 		message(myname)
	# 		print(theList[[myname]][1:min(5,length(theList[[myname]]))])
	# 	}
	# }
}

objectSummary <- function(theObj)
{
	# [rows, cols]
	# for(myname in slotNames(theObj))
	# {
	# 	message(myname)
	# 	if (is.matrix(slot(theObj, myname)))
	# 	{
	# 		message(myname)
	# 		if (sum(dim(slot(theObj, myname)))>0)
	# 		{
	# 			print(slot(theObj, myname)[1:min(4,dim(slot(theObj, myname))[1]), 1:min(4,dim(slot(theObj, myname))[2])])
	# 		}
	# 	}
	# 	else if (is.data.frame(slot(theObj, myname)))
	# 	{
	# 		message(myname)
	# 		if (sum(dim(slot(theObj, myname)))>0)
	# 		{
	# 			print(slot(theObj, myname)[1:min(5,dim(slot(theObj, myname))[1]),1:min(4,dim(slot(theObj, myname))[2])])
	# 		}
	# 	}
	# 	else if (is.list(slot(theObj, myname)))
	# 	{
	# 		message(myname)
	# 		listSummary(slot(theObj, myname))
	# 	}
	# 	else if (is.vector(slot(theObj, myname)))
	# 	{
	# 		message(myname)
	# 		print(slot(theObj, myname)[1:min(5,length(slot(theObj, myname)))])
	# 	}
	# }
}

addMissingRows <- function(theMatAdd, theMatStatic)
{
  myrows <- rownames(theMatAdd)
  cmprow <- rownames(theMatStatic)
  addrows <- setdiff(cmprow, myrows)
  if(length(addrows)>0)
  {
    narow <- rep(x=NA, times=length(colnames(theMatAdd)))
    #print(narow)
    for(newrow in addrows)
    {
      #print(newrow)
      theMatAdd <- rbind(theMatAdd, narow)
    }
    rownames(theMatAdd) <- c(myrows, addrows)
  }
  theMatAdd
}

EBNPlus <- function(theData1,
										theData2,
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
										theEBNP_ParametricPriorsFlag,
										theEBNP_ValidationRatio,
										theEBNP_TestRatio,
										theEBNP_PriorPlotsFile,
										theEBNP_MinSampleNum=3,
										theEBNP_AddData1Rows=FALSE,
										theEBNP_AddData2Rows=FALSE,
										theTestSeed=NULL)
{
  logDebug("EBNPlus theEBNP_AddData1Rows=", theEBNP_AddData1Rows)
  logDebug("EBNPlus theEBNP_AddData2Rows=", theEBNP_AddData2Rows)
  # changes here may need to be reflected in both BEA_CorrectionsEBNP.R (EBNPlus) and BEA_EBNPapi.R (EBNPlus_TrainAndValidateFromVector_Structures)
	logDebug("EBNPlus theEBNP_PriorPlotsFile=", theEBNP_PriorPlotsFile)
	resultList <- list()
	logDebug(paste("dim(theData1)=", paste(dim(theData1), sep = " ", collapse = " "), sep = ""))
	printMatrix(theData1)
	printElements(rownames(theData1)[1:100])
	logDebug(paste("dim(theData2)=", paste(dim(theData2), sep = " ", collapse = " "), sep = ""))
	printMatrix(theData2)
	printElements(rownames(theData2)[1:100])

	logDebug("remove unknown genes, that start with ?")
	theData1 <- theData1[!grepl("^[?]+", rownames(theData1)), ]
	theData2 <- theData2[!grepl("^[?]+", rownames(theData2)), ]

	logDebug("check on adding missing rows")
	if(TRUE==theEBNP_AddData1Rows)
	{
	  logDebug("add missing rows set 1")
	  theData1 <- addMissingRows(theData1, theData2)
	  logDebug("rownames 1=", rownames(theData1))
	}
	if(TRUE==theEBNP_AddData2Rows)
	{
	  logDebug("add missing rows set 2")
	  theData2 <- addMissingRows(theData2, theData1)
	  logDebug("rownames 2=", rownames(theData2))
	}


	logDebug(paste("dim(theData1)=", paste(dim(theData1), sep = " ", collapse = " "), sep = ""))
	printMatrix(theData1)
	printElements(rownames(theData1)[1:100])
	logDebug(paste("dim(theData2)=", paste(dim(theData2), sep = " ", collapse = " "), sep = ""))
	printMatrix(theData2)
	printElements(rownames(theData2)[1:100])

	logDebug("make EBNplus")
	ebObj  <- new("EBNplus",
								theData1,
								theData2,
								theEBNP_TrimBarcodesFunction,
								theEBNP_RemoveColDuplicatesFunction,
								theEBNP_TrimGenesFunction,
								theEBNP_RemoveRowDuplicatesFunction)
	logDebug("after EBNplus")
	logDebug(paste("dim(ebObj@mData1)=", paste(dim(ebObj@mData1), sep = " ", collapse = " "), sep = ""))
	printMatrix(ebObj@mData1)
	printElements(rownames(ebObj@mData1)[1:100])
	logDebug(paste("dim(ebObj@mData2)=", paste(dim(ebObj@mData2), sep = " ", collapse = " "), sep = ""))
	printMatrix(ebObj@mData2)
	printElements(rownames(ebObj@mData2)[1:100])
	EBNPlus_CheckData_Structures(ebObj@mData1, ebObj@mData2)

	ebObj@DF1batch <- theEBNP_Data1BatchId
	ebObj@DF2batch  <- theEBNP_Data2BatchId
	ebObj@batchWith0 <- theEBNP_BatchWithZero
	# here fix the HiSeq data. I didn’t fix it because in theory it may get better results in small replicates number. The sample code is for our paper.
	#ebObj@fixSet <- 1
	ebObj@fixSet <- theEBNP_FixDataSet
	ebObj@correct0 <- theEBNP_CorrectForZero
	logDebug("getBiComOrder")
	ebObj <- getBiComOrder(ebObj)

	logDebug(paste("dim(ebObj@mat1Com)=",paste(dim(ebObj@mat1Com), sep = " ", collapse = " "),sep = ""))
	printMatrix(ebObj@mat1Com)
	printElements(rownames(ebObj@mat1Com)[1:100])
	logDebug(paste("dim(ebObj@mat2Com)=",paste(dim(ebObj@mat2Com), sep = " ", collapse = " "),sep = ""))
	printMatrix(ebObj@mat2Com)
	printElements(rownames(ebObj@mat2Com)[1:100])

	logDebug("getValidationSet")
	#ebObj <- getValidationSet(ebObj, seed = theSeed, validationRatio = 0.33)
	ebObj <- getValidationSet(ebObj, seed = theSeed, validationRatio = theEBNP_ValidationRatio)
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

	#if you want to get the whole adjusted set, then the code should be:
	# use all the replicates as training set to get the maximum power.
	# instead of this code, set theEBNP_ValidationRatio to 0
	#logDebug("copy Com matrixes")
	#ebObj@mat1Train = ebObj@mat1Com
	#ebObj@mat2Train = ebObj@mat2Com
	# set results
	#### sort sample ids and features (genes)
	resultList[["TrainingSet1"]] <- sortMatrix(ebObj@mat1Train)
	resultList[["TrainingSet2"]] <- sortMatrix(ebObj@mat2Train)
	resultList[["ValidationSet1"]] <- sortMatrix(ebObj@mat1Validation)
	resultList[["ValidationSet2"]] <- sortMatrix(ebObj@mat2Validation)
	# add test set output!!! Yay!
	resultList[["TestSet1"]] <- sortMatrix(ebObj@mat1Test)
	resultList[["TestSet2"]] <- sortMatrix(ebObj@mat2Test)

	#logDebug("EBNPlus argument to train")
	#objectSummary(ebObj)
	logDebug("train")
	# does this need seed to be reset?
	testFlag <- FALSE
	if (theEBNP_TestRatio!=0)
	{
	  testFlag <- TRUE
	}
	objafterEB <- train(ebObj, par.prior = theEBNP_ParametricPriorsFlag, test = testFlag,
	                    theEBNP_PriorPlotsFile=theEBNP_PriorPlotsFile, minSampleNum=theEBNP_MinSampleNum)
	logDebug("after train")
	#logDebug("EBNPlus after train")
	#objectSummary(objafterEB)

	logDebug("EBadj")
	wholeEB <- NULL
	trainEB <- NULL
	testEB <- NULL
	validationEB <- NULL
	if (0==theEBNP_ValidationRatio)
	{
		logDebug("EBadj whole 1")
		objafterEB  <- EBadj(objafterEB , "whole")
		logDebug("objafterEB@wholeEB")
		wholeEB <- objafterEB@wholeEB
		#print(wholeEB)
	}
	else
	{
		#####################################################
		logDebug("EBadj validation for validation results")
		validEB  <- EBadj(objafterEB , "validation")
		logDebug("objafterEB@validationEB")
		validationEB <- validEB@validationEB
		#####################################################
		logDebug("EBadj train for training results")
		trainingResults  <- EBadj(objafterEB , "train")
		logDebug("trainingResults@trainEB")
		trainEB <- trainingResults@trainEB
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

	return(resultList)
}

# trim barcodes to 15 -- this gives through sample type,
# so that normals match normals, tumor matches tumor, etc.
# https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode
trimBarcodes <- function(theBarcodes, theBarcodeLength=15)
{
	logDebug("trimBarcodes")
	substr(theBarcodes, 1, theBarcodeLength)
}

keepBarcodes <- function(theBarcodes)
{
	theBarcodes
}

trimGenes <- function(theGenes)
{
	logDebug("trimGenes")
	foo <- as.vector(unlist(
		sapply(theGenes, function(theGene)
		{
			# keep the same if it starts with ?
			if (TRUE==grepl("^[?]+", theGene))
			{
				return(theGene)
			}
			else
			{
				# split on the | and take the first argument
				# this makes no change if no pipe
				return(strsplit(theGene, "|", fixed=TRUE)[[1]][1])
			}
		})
	))
	foo
}

# remove duplicates from within same matrix
# TODO if this really removes a lot, maybe test trimming
# barcode to 16 (vial) or even 19 (portion)
removeDuplicatesFromColumns <- function(theMatrix)
{
	logDebug("removeDuplicatesFromColumns")
	indexOfDuplicates <- which(duplicated(colnames(theMatrix)))
	if (length(indexOfDuplicates) > 0)
	{
		# minus sign uses inverse of indexes
		theMatrix <- theMatrix[ ,-indexOfDuplicates]
	}
	return(theMatrix)
}

removeDuplicatesFromRows <- function(theMatrix)
{
	logDebug("removeDuplicatesFromRows")
	indexOfDuplicates <- which(duplicated(rownames(theMatrix)))
	if (length(indexOfDuplicates) > 0)
	{
		# minus sign uses inverse of indexes
		theMatrix <- theMatrix[-indexOfDuplicates, ]
	}
	return(theMatrix)
}

######################################################################

