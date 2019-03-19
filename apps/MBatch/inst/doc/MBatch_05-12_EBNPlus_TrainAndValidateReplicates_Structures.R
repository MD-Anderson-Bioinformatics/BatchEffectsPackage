## ---- echo=TRUE----------------------------------------------------------
{
  library(MBatch)
  
  theDataFile1="/bea_testing/MATRIX_DATA/brca_rnaseq2_matrix_data.tsv"
  theDataFile2="/bea_testing/MATRIX_DATA/brca_agi4502_matrix_data.tsv"
  theOutputDir="/bea_testing/output/EBNPlus_TrainAndValidateReplicates_Structures"
  theBatchId1="RNASeqV2"
  theBatchId2="Agilent4502"
  theRandomSeed=314
  
  # trim genes to get just gene symbols from standardized data
  trimGenes <- function(theGenes)
  {
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
  
  # remove duplicates from columns (samples)
  removeDuplicatesFromColumns <- function(theMatrix)
  {
    indexOfDuplicates <- which(duplicated(colnames(theMatrix)))
    if (length(indexOfDuplicates) > 0)
    {
      # minus sign uses inverse of indexes
      theMatrix <- theMatrix[ ,-indexOfDuplicates]
    }
    return(theMatrix)
  }
  
  # remove duplicates from rows (genes/probes)
  removeDuplicatesFromRows <- function(theMatrix)
  {
    indexOfDuplicates <- which(duplicated(rownames(theMatrix)))
    if (length(indexOfDuplicates) > 0)
    {
      # minus sign uses inverse of indexes
      theMatrix <- theMatrix[-indexOfDuplicates, ]
    }
    return(theMatrix)
  }
  
  
  printMatrix <- function(theMatrix)
  {
    print(is.matrix(theMatrix))
    print(dim(theMatrix))
    rowMax <- dim(theMatrix)[1]
    colMax <- dim(theMatrix)[2]
    rowMax <- min(rowMax, 4)
    colMax <- min(colMax, 4)
    print(theMatrix[1:rowMax, 1:colMax])
  }
  
  if ((!dir.exists(theDataFile1))&&(!dir.exists(theDataFile2)))
  {
    warnLevel<-getOption("warn")
    on.exit(options(warn=warnLevel))
    # warnings are errors
    options(warn=3)
    # if there is a warning, show the calls leading up to it
    options(showWarnCalls=TRUE)
    # if there is an error, show the calls leading up to it
    options(showErrorCalls=TRUE)
    #
    unlink(theOutputDir, recursive=TRUE)
    dir.create(theOutputDir, showWarnings=FALSE, recursive=TRUE)
    # read the files in. This can be done however you want
    print("read the files")
    theDataMatrix1 <- readAsGenericMatrix(theDataFile1)
    theDataMatrix2 <- readAsGenericMatrix(theDataFile2)
    # this is the reduce genes to just gene symbols, handling those from standardized data
    print("reduce to gene symbols")
    rownames(theDataMatrix1) <- trimGenes(rownames(theDataMatrix1))
    rownames(theDataMatrix2) <- trimGenes(rownames(theDataMatrix2))
    # remove any duplicates (this is a requirement for EBNplus)
    print("remove duplicates")
    theDataMatrix1 <- removeDuplicatesFromColumns(removeDuplicatesFromRows(theDataMatrix1))
    theDataMatrix2 <- removeDuplicatesFromColumns(removeDuplicatesFromRows(theDataMatrix2))
    print("EBNPlus_TrainAndValidateReplicates_Structures")
    resultsList <- EBNPlus_TrainAndValidateReplicates_Structures(
      theDataMatrix1, theDataMatrix2, theBatchId1, theBatchId2,
      theEBNP_BatchWithZero="1",
      theEBNP_FixDataSet=as.numeric(NA),
      theEBNP_CorrectForZero=TRUE,
      theEBNP_ParametricPriorsFlag=TRUE,
      theEBNP_ValidationRatio=0.3,
      theEBNP_TestRatio=0.3,
      theSeed=theRandomSeed,
      theTestSeed=theRandomSeed,
      theEBNP_PriorPlotsFile=file.path(theOutputDir, "priorplots.PNG"))
    print("TestSet1")
    printMatrix(resultsList$TestSet1)
    print("TestSet2")
    printMatrix(resultsList$TestSet2)
    print("TrainingSet1")
    printMatrix(resultsList$TrainingSet1)
    print("TrainingSet2")
    printMatrix(resultsList$TrainingSet2)
    print("TrainingResults")
    printMatrix(resultsList$TrainingResults)
    print("ValidationSet1")
    printMatrix(resultsList$ValidationSet1)
    print("ValidationSet2")
    printMatrix(resultsList$ValidationSet2)
    print("ValidationResults")
    printMatrix(resultsList$ValidationResults)
    print("CorrectedResults")
    printMatrix(resultsList$CorrectedResults)
  }
}

