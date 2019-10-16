#MBatchUtils Copyright ? 2018 University of Texas MD Anderson Cancer Center
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(MBatch)

#############################################################################
#### utility functions

# write matrix files out
writeMatrixToFile <- function(theMatrix, theFile)
{
  write.table(theMatrix, file=theFile, quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)
}

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

convertNulls <- function(theString)
{
  if("null"==theString)
  {
    theString <- NULL
  }
  else if("NULL"==theString)
  {
    theString <- NULL
  }
  else if(""==theString)
  {
    theString <- NULL
  }
  theString
}

convertStringArray <- function(theString)
{
  if (!is.null(theString))
  {
    theString <- strsplit(theString, ",")[[1]]
  }
  theString
}


sortMatrix <- function(theMatrix)
{
  if (!is.null(theMatrix))
  {
    if (dim(theMatrix)[1]>0)
    {
      if (dim(theMatrix)[2]>0)
      {
        # order rows, columns
        theMatrix <- theMatrix[order(rownames(theMatrix)), order(colnames(theMatrix))]
      }
    }
  }
  theMatrix
}

sortDataframe <- function(theDataframe)
{
  if (!is.null(theDataframe))
  {
    if (dim(theDataframe)[1]>0)
    {
      if (dim(theDataframe)[2]>0)
      {
        theDataframe <- theDataframe[order(theDataframe$Sample),]
      }
    }
  }
  theDataframe
}

convertNA <- function(theData)
{
  if (is.na(theData))
  {
    theData <- ""
  }
  theData
}

mbatchRunFromConfig <- function(theConfigFile, theOutputDir, theNaStrings, theShaidyMapGen, theShaidyMapGenJava,
                                theNGCHMShaidyMem="16G", thePCAMem="4800m", theBoxplotMem="16G")
{
  message(mbatchUtilVersion())
  collateOrigValue<-Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
  Sys.setlocale("LC_COLLATE","C")
  message("Changing LC_COLLATE to C for duration of run")
  ####################################################################
  myConfig <- read.csv(theConfigFile, header=FALSE, sep="\t", as.is=TRUE, row.names=1 )
  rbnOnly <- as.logical(convertNulls(myConfig["RBN_Only",]))
  if (isTRUE(rbnOnly))
  {
    runRBNfromConfig(theConfigFile, theOutputDir)
  }
  else
  {
    title <- myConfig["title",]
    sampleidBatchType <- myConfig["sampleidBatchType",]
    batchTypesForMBatch <- myConfig["batchTypesForMBatchArray",]
    filterMaxValue <- myConfig["filterMaxValue",]
    filterLogTransformFlag <- myConfig["filterLogTransformFlag",]
    filterLogTransformFlag2 <- myConfig["filterLogTransformFlag2",]
    filteringBatchType <- myConfig["filteringBatchType",]
    filteringBatches <- convertStringArray(convertNulls(myConfig["filteringBatchesArray",]))
    CDP_Flag <- as.logical(convertNulls(myConfig["CDP_Flag",]))
    message("title ", title)
    message("sampleidBatchType ", sampleidBatchType)
    message("batchTypesForMBatch ", batchTypesForMBatch)
    message("filterMaxValue ", filterMaxValue)
    message("filterLogTransformFlag ", filterLogTransformFlag)
    message("filterLogTransformFlag2 ", filterLogTransformFlag2)
    message("filteringBatchType ", filteringBatchType)
    message("filteringBatches ", filteringBatches)
    message("CDP_Flag ", CDP_Flag)
    selectedBatchToCorrect <- myConfig["selectedBatchToCorrect",]
    selectedCorrection <- myConfig["selectedCorrection",]
    selectedCorrectionMinBatchSize <- myConfig["selectedCorrectionMinBatchSize",]
    message("selectedBatchToCorrect ", selectedBatchToCorrect)
    message("selectedCorrection ", selectedCorrection)
    message("selectedCorrectionMinBatchSize ", selectedCorrectionMinBatchSize)
    selectedDSCPermutations <- myConfig["selectedDSCPermutations",]
    selectedDSCThreads <- myConfig["selectedDSCThreads",]
    selectedDSCMinBatchSize <- myConfig["selectedDSCMinBatchSize",]
    selectedDSCSeed <- myConfig["selectedDSCSeed",]
    selectedDSCMaxGeneCount <- myConfig["selectedDSCMaxGeneCount",]
    message("selectedDSCPermutations ", selectedDSCPermutations)
    message("selectedDSCThreads ", selectedDSCThreads)
    message("selectedDSCMinBatchSize ", selectedDSCMinBatchSize)
    message("selectedDSCSeed ", selectedDSCSeed)
    message("selectedDSCMaxGeneCount ", selectedDSCMaxGeneCount)
    selectedBoxplotMaxGeneCount <- myConfig["selectedBoxplotMaxGeneCount",]
    selectedNgchmFlag <- as.logical(convertNulls(convertNA(myConfig["selectedNgchmFlag",])))
    selectedMutBatchFlag <- as.logical(convertNulls(convertNA(myConfig["selectedMutBatchFlag",])))
    message("selectedBoxplotMaxGeneCount ", selectedBoxplotMaxGeneCount)
    message("selectedNgchmFlag ", selectedNgchmFlag)
    #message("selectedMutBatchFlag ", selectedMutBatchFlag)
    RBN_UseFirstAsInvariant <- as.logical(convertNulls(convertNA(myConfig["RBN_UseFirstAsInvariantFlag",])))
    RBN_InvariantId <- convertNulls(myConfig["RBN_InvariantId",])
    RBN_VariantId <- convertNulls(myConfig["RBN_VariantId",])
    RBN_Matched <- as.logical(convertNulls(myConfig["RBN_Matched",]))
    RBN_InvariantRepsType <- convertNulls(myConfig["RBN_InvariantRepsType",])
    RBN_VariantRepsType <- convertNulls(myConfig["RBN_VariantRepsType",])
    RBN_InvariantReps <- convertStringArray(convertNulls(myConfig["RBN_InvariantRepsArray",]))
    RBN_VariantReps <- convertStringArray(convertNulls(myConfig["RBN_VariantRepsArray",]))
    message("RBN_UseFirstAsInvariant ", RBN_UseFirstAsInvariant)
    message("RBN_InvariantId ", RBN_InvariantId)
    message("RBN_VariantId ", RBN_VariantId)
    message("RBN_Matched ", RBN_Matched)
    message("RBN_InvariantRepsType ", RBN_InvariantRepsType)
    message("RBN_VariantRepsType ", RBN_VariantRepsType)
    message("RBN_InvariantReps ", RBN_InvariantReps)
    message("RBN_VariantReps ", RBN_VariantReps)
    EBNPlus_GroupId1 <- convertNulls(myConfig["EBNPlus_GroupId1",])
    EBNPlus_GroupId2 <- convertNulls(myConfig["EBNPlus_GroupId2",])
    EBNPlus_Seed <- as.integer(convertNulls(myConfig["EBNPlus_Seed",]))
    EBNPlus_MinSamples <- as.integer(convertNulls(myConfig["EBNPlus_MinSamples",]))
    message("EBNPlus_GroupId1 ", EBNPlus_GroupId1)
    message("EBNPlus_GroupId2 ", EBNPlus_GroupId2)
    message("EBNPlus_Seed ", EBNPlus_Seed)
    message("EBNPlus_MinSamples ", EBNPlus_MinSamples)
    ####################################################################
    if("null"==filteringBatchType)
    {
      filteringBatchType <- NULL
    }
    else if(""==filteringBatchType)
    {
      filteringBatchType <- NULL
    }
    ####
    if ((!is.null(filteringBatches)) && ("null"==filteringBatches))
    {
      filteringBatches <- NULL
    }
    else if ((!is.null(filteringBatches)) && (""==filteringBatches))
    {
      filteringBatches <- NULL
    }
    ####
    if(("none"==selectedCorrection)||("null"==selectedCorrection)||(""==selectedCorrection))
    {
      selectedCorrection <- NULL
    }
    ####
    if (!is.null(filteringBatches))
    {
      filteringBatches <- strsplit(filteringBatches, ",")[[1]]
    }
    if (!is.null(batchTypesForMBatch))
    {
      batchTypesForMBatch <- strsplit(batchTypesForMBatch, ",")[[1]]
    }
    filterMaxValue <- as.numeric(filterMaxValue)
    filterLogTransformFlag <- as.logical(convertNulls(convertNA(filterLogTransformFlag)))
    filterLogTransformFlag2 <- as.logical(convertNulls(convertNA(filterLogTransformFlag2)))
    selectedCorrectionMinBatchSize <- as.numeric(selectedCorrectionMinBatchSize)
    selectedDSCPermutations <- as.numeric(selectedDSCPermutations)
    selectedDSCThreads <- as.numeric(selectedDSCThreads)
    selectedDSCMinBatchSize <- as.numeric(selectedDSCMinBatchSize)
    selectedDSCSeed <- as.numeric(selectedDSCSeed)
    selectedDSCMaxGeneCount <- as.numeric(selectedDSCMaxGeneCount)
    selectedBoxplotMaxGeneCount <- as.numeric(selectedBoxplotMaxGeneCount)
    ####################################################################
    sourceDir <-dirname(theConfigFile)
    logFile <- file.path(sourceDir, "mbatch.log")
    datFile <- file.path(sourceDir, "matrix_data.tsv")
    batFile <- file.path(sourceDir, "batches.tsv")
    datFile2 <- file.path(sourceDir, "matrix_data2.tsv")
    batFile2 <- file.path(sourceDir, "batches2.tsv")
    #############################################################################
    # check directories
    message("theOutputDir=", theOutputDir)
    if(!dir.exists(theOutputDir))
    {
      message("create ", theOutputDir)
      dir.create(theOutputDir, showWarnings=FALSE, recursive=TRUE)
    }
    #############################################################################
    setLogging(new("Logging", theFile=logFile))
    #############################################################################
    message("rename data")
    origDatFile <- paste(datFile, ".bak", sep="")
    origBatFile <- paste(batFile, ".bak", sep="")
    origDatFile2 <- paste(datFile2, ".bak", sep="")
    origBatFile2 <- paste(batFile2, ".bak", sep="")
    file.rename(datFile, origDatFile)
    file.rename(batFile, origBatFile)
    if(file.exists(datFile2))
    {
      file.rename(datFile2, origDatFile2)
      file.rename(batFile2, origBatFile2)
    }
    message("preprocess data")
    filterList <- list()
    if (!is.null(filteringBatchType))
    {
      filterList <- lapply(filteringBatches, function(x) { c(filteringBatchType, x)})
    }
    preprocessData(origDatFile, origBatFile, datFile, batFile, theSize=0, theTransformFlag=filterLogTransformFlag, theBatchTypeAndValuePairsToRemove=filterList)
    if(file.exists(origDatFile2))
    {
      preprocessData(origDatFile2, origBatFile2, datFile2, batFile2, theSize=0, theTransformFlag=filterLogTransformFlag2)
    }
    #############################################################################
    #### do data load
    #############################################################################
    message("load data as needed")
    # TODO: HERE FIX FILES
    myOriginalData <- NULL
    myOriginalData2 <- NULL
    myMBatchData <- NULL
    myMBatchData2 <- NULL
    file.exists(datFile)
    if ((!is.null(selectedCorrection)) && ( startsWith(selectedCorrection, "RBN_") ))
    {
      message("Read RBN data")
      if (isTRUE(RBN_UseFirstAsInvariant))
      {
        # TODO: use , theNaStrings=theNaStrings after changing pre-set values
        myMBatchData <- mbatchLoadFiles(datFile, batFile)
        myMBatchData2 <- mbatchLoadFiles(datFile2, batFile2)
      }
      else
      {
        # TODO: use , theNaStrings=theNaStrings after changing pre-set values
        myMBatchData <- mbatchLoadFiles(datFile2, batFile2)
        myMBatchData2 <- mbatchLoadFiles(datFile, batFile)
      }
    }
    else if ((!is.null(selectedCorrection)) && ( startsWith(selectedCorrection, "EBN_") ))
    {
      message("Read EBNplus data")
      # EBNplus cannot use the mbatchLoadFiles function, since it filters meaningfule data
      matrix1 <- readAsGenericMatrix(datFile)
      matrix2 <- readAsGenericMatrix(datFile2)
      df1 <- readAsDataFrame(batFile)
      df2 <- readAsDataFrame(batFile2)
      # make genes nice
      rownames(matrix1) <- trimGenes(rownames(matrix1))
      rownames(matrix2) <- trimGenes(rownames(matrix2))
      # remove any duplicates (this is a requirement for EBNplus)
      matrix1 <- removeDuplicatesFromColumns(removeDuplicatesFromRows(matrix1))
      matrix2 <- removeDuplicatesFromColumns(removeDuplicatesFromRows(matrix2))
      # build objects
      myMBatchData<-new("BEA_DATA", sortMatrix(matrix1), sortDataframe(df1), data.frame())
      myMBatchData2<-new("BEA_DATA", sortMatrix(matrix2), sortDataframe(df2), data.frame())
    }
    else
    {
      message("Read regular data")
      message(datFile)
      message(batFile)
      # TODO: use , theNaStrings=theNaStrings after changing pre-set values
      myMBatchData <- mbatchLoadFiles(datFile, batFile)
    }
    # rename sample id column to Sample
    if("Sample"!=sampleidBatchType)
    {
      names(myMBatchData@mBatches)[names(myMBatchData@mBatches)==sampleidBatchType] <- "Sample"
      sampleidBatchType <- "Sample"
    }
    newBatchTypesForMBatch <- c()
    goodBatchTypes <- names(myMBatchData@mBatches)
    for(val in batchTypesForMBatch)
    {
      if (val %in% goodBatchTypes)
      {
        newBatchTypesForMBatch <- c(newBatchTypesForMBatch, val)
      }
    }
    batchTypesForMBatch <- newBatchTypesForMBatch
    #### do corrections
    if (isTRUE(CDP_Flag))
    {
      myOriginalData <- myMBatchData
      if (!is.null(myMBatchData2))
      {
        myOriginalData2 <- myMBatchData2
      }
    }
    if (!is.null(selectedCorrection))
    {
      myMBatchData <- doCorrectionsFromConfig(theOutputDir=theOutputDir,
                                              theDataObject=myMBatchData,
                                              theDataObject2=myMBatchData2,
                                              theSelectedBatchToCorrect=selectedBatchToCorrect,
                                              theSelectedCorrection=selectedCorrection,
                                              thePermutationThreads=selectedDSCThreads,
                                              theNaStrings=theNaStrings,
                                              theRBN_InvariantId=RBN_InvariantId,
                                              theRBN_VariantId=RBN_VariantId,
                                              theRBN_Matched=RBN_Matched,
                                              theRBN_InvariantReps=RBN_InvariantReps,
                                              theRBN_VariantReps=RBN_VariantReps,
                                              theEBNPlus_GroupId1=EBNPlus_GroupId1,
                                              theEBNPlus_GroupId2=EBNPlus_GroupId2,
                                              theEBNPlus_Seed=EBNPlus_Seed,
                                              theEBNPlus_MinSamples=EBNPlus_MinSamples)
    }
    #### do filtering
    batchesToRemove <- setdiff(colnames(myMBatchData@mBatches), c(sampleidBatchType, batchTypesForMBatch))
    myMBatchData <- mbatchFilterData(theBeaData=myMBatchData,
                                     theBatchTypeAndValuePairsToRemove=list(),
                                     theBatchTypeAndValuePairsToKeep=list(),
                                     theBatchTypesToRemove=batchesToRemove,
                                     theMinIqr=0,
                                     theMinSd=0,
                                     theMinMad=0)
    if (isTRUE(CDP_Flag))
    {
      batchesToRemove <- setdiff(colnames(myOriginalData@mBatches), c(sampleidBatchType, batchTypesForMBatch))
      myOriginalData <- mbatchFilterData(theBeaData=myOriginalData,
                                       theBatchTypeAndValuePairsToRemove=list(),
                                       theBatchTypeAndValuePairsToKeep=list(),
                                       theBatchTypesToRemove=batchesToRemove,
                                       theMinIqr=0,
                                       theMinSd=0,
                                       theMinMad=0)
    }
    #### do assessments
    rdataHC <- doAssessmentsFromConfig(theOutputDir=theOutputDir,
                            theDataObject=myMBatchData,
                            theTitle=title,
                            thePermutations=selectedDSCPermutations,
                            thePermutationThreads=selectedDSCThreads,
                            theMaxSize=filterMaxValue,
                            theSeed=selectedDSCSeed,
                            theDSCMaxGenes=selectedDSCMaxGeneCount,
                            theBoxplotMaxGenes=selectedBoxplotMaxGeneCount,
                            theMinBatchSize=selectedDSCMinBatchSize,
                            thePCAMem=thePCAMem,
                            theBoxplotMem=theBoxplotMem)
    #############################################################################
    if (isTRUE(selectedNgchmFlag))
    {
      if (!is.null(rdataHC))
      {
        warnLevel<-getOption("warn")
        on.exit(options(warn=warnLevel))
        options(warn=-1)
        for(myBatchType in batchTypesForMBatch)
        {
          message("title ", title)
          message("myBatchType ", myBatchType)
          message("theShaidyMapGen ", theShaidyMapGen)
          message("theShaidyMapGenJava ", theShaidyMapGenJava)
          message("dim(myOriginalData@mData) ", dim(myOriginalData@mData))
          message("dim(myOriginalData@mBatches) ", dim(myOriginalData@mBatches))
          buildBatchHeatMapFromHC_Structures(theMatrixData=myOriginalData@mData,
                                       theBatchData=myOriginalData@mBatches,
                                       theTitle=paste(title, myBatchType, sep=" "),
                                       theOutputFile=file.path(theOutputDir, "NGCHM", paste(myBatchType, "ngchm.ngchm", sep="_")),
                                       theUDendRDataFile=rdataHC,
                                       theRowType="labels", theColType="bio.tcga.barcode.sample",
                                       theRowCluster=c("pearson", "ward.D2"),
                                       theShaidyMapGen=theShaidyMapGen,
                                       theShaidyMapGenJava=theShaidyMapGenJava,
                                       theShaidyMapGenArgs=c(paste(c("-Xms", "-Xmx"), theNGCHMShaidyMem, sep=""), "-Djava.awt.headless=true"))
        }
      }
      else
      {
        message("no HC data, so skip NGCHM")
      }
    }
    #############################################################################
    if (isTRUE(CDP_Flag))
    {
      graphics.off()
      warnLevel<-getOption("warn")
      on.exit(options(warn=warnLevel))
      options(warn=-1)
      if (!file.exists(file.path(theOutputDir, "CDP")))
      {
        dir.create(file.path(theOutputDir, "CDP"),showWarnings = FALSE, recursive = TRUE)
      }
      CDP_Structures(file.path(theOutputDir, "CDP", "CDP_Plot_Data1_Diagram.PNG"), myOriginalData@mData, myMBatchData@mData, "CDP Output 1")
      if (!is.null(myOriginalData2))
      {
        CDP_Structures(file.path(theOutputDir, "CDP", "CDP_Plot_Data2_Diagram.PNG"), myOriginalData2@mData, myMBatchData@mData, "CDP Output 2")
      }
    }
  }
  #############################################################################
  message("Write success note to '", file.path(theOutputDir, "MBATCH_SUCCESS.txt"), "'")
  file.create(file.path(theOutputDir, "MBATCH_SUCCESS.txt"))
}


doAssessmentsFromConfig <- function(theOutputDir, theDataObject, theTitle,
                                    thePermutations, thePermutationThreads, theMaxSize,
                                    theSeed, theDSCMaxGenes, theBoxplotMaxGenes,
                                    theMinBatchSize, thePCAMem, theBoxplotMem)
{
  # reduce size
  theDataObject@mData <- mbatchTrimData(theDataObject@mData, theMaxSize)
  ####
  #### HierarchicalClustering
  ####
  rdataHC <- callMBatch_HierarchicalClustering_Structures(theOutputDir, theDataObject, theTitle)
  ####
  #### SupervisedClustering
  ####
  callMBatch_SupervisedClustering_Structures(theOutputDir, theDataObject, theTitle)
  ####
  #### PCAPlus
  ####
  callMBatch_PCA_Structures(theOutputDir, theDataObject, theTitle,
                            theIsPcaTrendFunction=isTrendBatch,
                            theDSCPermutations=thePermutations,
                            theDSCThreads=thePermutationThreads,
                            theMinBatchSize=theMinBatchSize,
                            theJavaParameters=c(paste(c("-Xms", "-Xmx"), thePCAMem, sep=""), "-Djava.awt.headless=true"),
                            theSeed=theSeed,
                            theMaxGeneCount=theDSCMaxGenes)
  ####
  #### BoxPlot
  ###
  pipelineMean <- function(x)
  {
    mean(x, na.rm=TRUE)
  }
  if (theBoxplotMaxGenes>5000)
  {
    theBoxplotMaxGenes <- 5000
  }
  callMBatch_BoxplotAllSamplesData_Structures(theOutputDir, theDataObject, theTitle,
                                              theMaxGeneCount=theBoxplotMaxGenes)
  callMBatch_BoxplotAllSamplesRLE_Structures(theOutputDir, theDataObject, theTitle,
                                             theMaxGeneCount=theBoxplotMaxGenes)
  callMBatch_BoxplotGroup_Structures(theOutputDir, theDataObject, theTitle,
                                     theMaxGeneCount=theBoxplotMaxGenes,
                                     theFunction=list(pipelineMean), theFunctionName=list("Mean"))
  if (!file.exists(file.path(theOutputDir, "BatchData.tsv")))
  {
    writeAsDataframe(file.path(theOutputDir, "BatchData.tsv"), theDataObject@mBatches)
  }
  message("Processing success")
  mbatchWriteSuccessfulLog()
  rdataHC
}

doCorrectionsFromConfig <- function(theOutputDir, theDataObject, theDataObject2,
                                    theSelectedBatchToCorrect, theSelectedCorrection, thePermutationThreads, theNaStrings,
                                    theRBN_InvariantId, theRBN_VariantId, theRBN_Matched, theRBN_InvariantReps, theRBN_VariantReps,
                                    theEBNPlus_GroupId1, theEBNPlus_GroupId2, theEBNPlus_Seed, theEBNPlus_MinSamples)
{
  processed <- FALSE
  myCorrectedFile <- NULL
  dataBatches <- NULL
  dataMatrix <- NULL
  if ("EB_withPara"==theSelectedCorrection)
  {
    processed <- TRUE
    myCorrectedFile <- EB_withParametricPriors(theBeaData=theDataObject,
                                               theBatchIdsNotToCorrect=c(""),
                                               theDoCheckPlotsFlag=TRUE,
                                               theBatchType=theSelectedBatchToCorrect,
                                               theThreads=thePermutationThreads,
                                               thePath=theOutputDir,
                                               theWriteToFile=TRUE)
    dataBatches <- theDataObject@mBatches
    if(!is.null(myCorrectedFile))
    {
      dataMatrix <- readAsGenericMatrix(myCorrectedFile)
    }
  }
  else if ("EB_withNonpara"==theSelectedCorrection)
  {
    processed <- TRUE
    myCorrectedFile <- EB_withNonParametricPriors(theBeaData=theDataObject,
                                                  theBatchIdsNotToCorrect=c(""),
                                                  theDoCheckPlotsFlag=TRUE,
                                                  theBatchType=theSelectedBatchToCorrect,
                                                  theThreads=thePermutationThreads,
                                                  thePath=theOutputDir,
                                                  theWriteToFile=TRUE)
    dataBatches <- theDataObject@mBatches
    if(!is.null(myCorrectedFile))
    {
      dataMatrix <- readAsGenericMatrix(myCorrectedFile)
    }
  }
  else if ("MP_overall"==theSelectedCorrection)
  {
    processed <- TRUE
    myCorrectedFile <- MP_Overall(theBeaData=theDataObject,
                                  thePath=theOutputDir,
                                  theWriteToFile=TRUE)
    dataBatches <- theDataObject@mBatches
    if(!is.null(myCorrectedFile))
    {
      dataMatrix <- readAsGenericMatrix(myCorrectedFile)
    }
  }
  else if ("MP_batch"==theSelectedCorrection)
  {
    processed <- TRUE
    myCorrectedFile <- MP_ByBatch(theBeaData=theDataObject,
                                  theBatchType=theSelectedBatchToCorrect,
                                  thePath=theOutputDir,
                                  theWriteToFile=TRUE)
    dataBatches <- theDataObject@mBatches
    if(!is.null(myCorrectedFile))
    {
      dataMatrix <- readAsGenericMatrix(myCorrectedFile)
    }
  }
  else if ("ANOVA_adj"==theSelectedCorrection)
  {
    processed <- TRUE
    myCorrectedFile <- AN_Adjusted(theBeaData=theDataObject,
                                   theBatchType=theSelectedBatchToCorrect,
                                   thePath=theOutputDir,
                                   theWriteToFile=TRUE)
    dataBatches <- theDataObject@mBatches
    if(!is.null(myCorrectedFile))
    {
      dataMatrix <- readAsGenericMatrix(myCorrectedFile)
    }
  }
  else if ("ANOVA_unadj"==theSelectedCorrection)
  {
    processed <- TRUE
    myCorrectedFile <- AN_Unadjusted(theBeaData=theDataObject,
                                     theBatchType=theSelectedBatchToCorrect,
                                     thePath=theOutputDir,
                                     theWriteToFile=TRUE)
    dataBatches <- theDataObject@mBatches
    if(!is.null(myCorrectedFile))
    {
      dataMatrix <- readAsGenericMatrix(myCorrectedFile)
    }
  }
  else if ("RBN_Replicates"==theSelectedCorrection)
  {
    #theDataObject, theDataObject2, theRBN_InvariantId, theRBN_VariantId, theRBN_Matched,
    processed <- TRUE
    myCorrectedFile <- RBN_Replicates(theInvariantMatrix=theDataObject@mData,
                               theVariantMatrix=theDataObject2@mData,
                               theInvariantGroupId=theRBN_InvariantId,
                               theVariantGroupId=theRBN_VariantId,
                               theMatchedReplicatesFlag=theRBN_Matched,
                               theCombineOnlyFlag=FALSE,
                               thePath=theOutputDir,
                               theWriteToFile=TRUE)
    # TODO: confirm this works for RBN
    dataBatches <- EBNPlus_CombineBatches(theBeaBatches1=theDataObject@mBatches,
                                          theBeaBatches2=theDataObject2@mBatches,
                                          theEBNP_Data1BatchId=theRBN_InvariantId,
                                          theEBNP_Data2BatchId=theRBN_VariantId,
                                          theBarcodeTrimFunction=trimGenes, theSep="-")
    writeAsDataframe(theFile=file.path(theOutputDir, "corrected_batches.tsv"), theDataframe=dataBatches)
    if(!is.null(myCorrectedFile))
    {
      dataMatrix <- readAsGenericMatrix(myCorrectedFile)
    }
  }
  else if ("RBN_Pseudoreps"==theSelectedCorrection)
  {
    #theDataObject, theDataObject2, theRBN_InvariantId, theRBN_VariantId, theRBN_Matched, theRBN_InvariantReps, theRBN_VariantReps,
    processed <- TRUE
    myCorrectedFile <- RBN_Pseudoreplicates(theInvariantMatrix=theDataObject@mData,
                                            theVariantMatrix=theDataObject2@mData,
                                            theInvariantReplicates = theRBN_InvariantReps,
                                            theVariantReplicates = theRBN_VariantReps,
                                            theInvariantGroupId=theRBN_InvariantId,
                                            theVariantGroupId=theRBN_VariantId,
                                            theMatchedReplicatesFlag=theRBN_Matched,
                                            theCombineOnlyFlag=FALSE,
                                            thePath=theOutputDir,
                                            theWriteToFile=TRUE)
    # TODO: confirm this works for RBN
    dataBatches <- EBNPlus_CombineBatches(theBeaBatches1=theDataObject@mBatches,
                                          theBeaBatches2=theDataObject2@mBatches,
                                          theEBNP_Data1BatchId=theRBN_InvariantId,
                                          theEBNP_Data2BatchId=theRBN_VariantId,
                                          theBarcodeTrimFunction=trimGenes, theSep="-")
    writeAsDataframe(theFile=file.path(theOutputDir, "corrected_batches.tsv"), theDataframe=dataBatches)
    if(!is.null(myCorrectedFile))
    {
      dataMatrix <- readAsGenericMatrix(myCorrectedFile)
    }
  }
  else if ("EBN_Plus"==theSelectedCorrection)
  {
    #theDataObject, theDataObject2, theEBNPlus_Batch1, theEBNPlus_Batch2, theEBNPlus_GroupId1, theEBNPlus_GroupId2, theEBNPlus_Seed,
    processed <- TRUE
    correctedFile <- file.path(theOutputDir, "ANY_Corrections-EBNPlus.tsv")
    print(dim(theDataObject@mData))
    print(dim(theDataObject2@mData))
    dataMatrix <- EBNPlus_Correction_Structures(theDataMatrix1=theDataObject@mData,
                                                theDataMatrix2=theDataObject2@mData,
                                                theBatchId1=theEBNPlus_GroupId1,
                                                theBatchId2=theEBNPlus_GroupId2,
                                                theEBNP_BatchWithZero="1",
                                                theEBNP_FixDataSet=as.numeric(NA),
                                                theEBNP_CorrectForZero=TRUE,
                                                theEBNP_ParametricPriorsFlag=TRUE,
                                                theSeed=theEBNPlus_Seed,
                                                theEBNP_MinSampleNum=theEBNPlus_MinSamples,
                                                theEBNP_PriorPlotsFile=file.path(theOutputDir, "priorplots_Diagram.PNG"))
    writeAsMatrix(correctedFile, dataMatrix)
    dataBatches <- EBNPlus_CombineBatches(theBeaBatches1=theDataObject@mBatches,
                                      theBeaBatches2=theDataObject2@mBatches,
                                      theEBNP_Data1BatchId=theEBNPlus_GroupId1,
                                      theEBNP_Data2BatchId=theEBNPlus_GroupId2,
                                      theBarcodeTrimFunction=trimGenes)
    writeAsDataframe(theFile=file.path(theOutputDir, "corrected_batches.tsv"), theDataframe=dataBatches)
  }
  else
  {
    stop("Unrecognized correction")
  }

  if ((!is.null(dataMatrix))&&(TRUE==processed))
  {
    # if have data info and processing was done, return new data object
    theDataObject <- mbatchLoadStructures(dataMatrix, dataBatches)
  }
  dataMatrix <- NULL
  dataBatches <- NULL
  theDataObject
}
