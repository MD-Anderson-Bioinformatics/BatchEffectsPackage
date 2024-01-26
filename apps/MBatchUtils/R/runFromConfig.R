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

#############################################################################
#### utility functions

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
  if (is.na(theString))
  {
    theString <- NULL
  }
  else if("null"==theString)
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

copyDataFiles <- function(theDatFile, theBatFile, theDatFile2, theBatFile2)
{
  message("rename data")
  origDatFile <- file.path(dirname(theDatFile), paste("original-", basename(theDatFile), sep=""))
  origBatFile <- file.path(dirname(theBatFile), paste("original-", basename(theBatFile), sep=""))
  # use copy instead of rename, since on mounted filesystems,
  # rename seems to not work properly.
  print(paste("copy ", theDatFile, " to ", origDatFile, sep=""))
  file.copy(theDatFile, origDatFile, overwrite=TRUE)
  print(paste("copy ", theBatFile, " to ", origBatFile, sep=""))
  file.copy(theBatFile, origBatFile, overwrite=TRUE)
  if(file.exists(theDatFile2))
  {
    origDatFile2 <- file.path(dirname(theDatFile2), paste("original-", basename(theDatFile2), sep=""))
    origBatFile2 <- file.path(dirname(theBatFile2), paste("original-", basename(theBatFile2), sep=""))
    print(paste("copy ", theDatFile2, " to ", origDatFile2, sep=""))
    file.copy(theDatFile2, origDatFile2, overwrite=TRUE)
    print(paste("copy ", theBatFile2, " to ", origBatFile2, sep=""))
    file.copy(theBatFile2, origBatFile2, overwrite=TRUE)
  }
}

populateNeededFiles <- function(theInfoDir, theJobId)
{
  message("populateNeededFiles - with default values")
  # if needed for insertion
  timestamp <- format(Sys.time(), "%Y_%m_%d_%H%M")
  dataStampFile <- file.path(theInfoDir, "data_stamp.txt")
  if (!file.exists(dataStampFile))
  {
    value <- paste("DATA_", timestamp, sep="")
    message("Set default test_stamp.txt to ", value)
    write(value, file=dataStampFile)
  }
  jobIdFile <- file.path(theInfoDir, "job_id.txt")
  if (!file.exists(jobIdFile))
  {
    value <- paste("JOB_", theJobId, sep="")
    message("Set default job_id.txt to ", value)
    write(value, file=jobIdFile)
  }
  originalDataFile <- file.path(theInfoDir, "original_data.json")
  if (!file.exists(originalDataFile))
  {
    value <- paste("{\n",
                   "\"source\": \"USER\",\n",
                   "\"program\": \"USER\",\n",
                   "\"project\": \"USER\",\n",
                   "\"category\": \"USER\",\n",
                   "\"platform\": \"USER\",\n",
                   "\"data\": \"USER\",\n",
                   "\"details\": \"USER\",\n",
                   "\"version\": \"",timestamp ,"\"\n",
                   "}", sep="")
    message("Set default original_data.json to ", value)
    write(value, file=originalDataFile)
  }
  sourceIdFile <- file.path(theInfoDir, "source_id.txt")
  if (!file.exists(sourceIdFile))
  {
    value <- paste("SRC_", timestamp, sep="")
    message("Set default source_id.txt to ", value)
    write(value, file=sourceIdFile)
  }
  testStampFile <- file.path(theInfoDir, "test_stamp.txt")
  if (!file.exists(testStampFile))
  {
    value <- paste("TEST_", timestamp, sep="")
    message("Set default test_stamp.txt to ", value)
    write(value, file=testStampFile)
  }
  versionTypeFile <- file.path(theInfoDir, "version_type.txt")
  if (!file.exists(versionTypeFile))
  {
    message("Set default version_type.txt to Original-Analyzed")
    write("Original-Analyzed", file=versionTypeFile)
  }
}

notEmpty <- function(theValue)
{
  ((theValue!="") && (!is.null(theValue)) && (!is.na(theValue)))
}


mbatchRunFromConfig <- function(theConfigFile, theMatrixFile, theBatchesFile,
                                theZipDataDir, theZipResultsDir, theNaStrings,
                                theShaidyMapGen, theNgchmWidgetJs, theShaidyMapGenJava,
                                theNGCHMShaidyMem="16G", theRunPostFlag=FALSE,
                                theMatrixFile2=NULL, theBatchesFile2=NULL)
{
  # used to flag if other than full success written
  otherNote <- FALSE
  message(mbatchUtilVersion())
  collateOrigValue<-Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
  Sys.setlocale("LC_COLLATE","C")
  message("Changing LC_COLLATE to C for duration of run")
  ####################################################################
  stopifnot(file.exists(theShaidyMapGen))
  stopifnot(file.exists(theNgchmWidgetJs))
  stopifnot(file.exists(theShaidyMapGenJava))
  ####################################################################
  myConfig <- read.csv(theConfigFile, header=FALSE, sep="\t", as.is=TRUE, row.names=1, quote="" )
  myDataVersion <- convertNulls(myConfig["DataVersion",])
  # empty means don't use stopifnot(notEmpty(myDataVersion))
  myTestVersion <- convertNulls(myConfig["TestVersion",])
  # empty means don't use stopifnot(notEmpty(myTestVersion))
  # get flag for replacing NAs with zero
  replaceNAs <- as.logical(convertNulls(convertNA(myConfig["replaceNAs",])))
  ########################
  # setup data directories
  ########################
  zipDataFiles <- setupZipDataDir(theMatrixFile, theBatchesFile, theZipDataDir, myDataVersion,
                                  theMatrixFile2, theBatchesFile2)
  # matrix data file
  currentMatrixFile <- convertNulls(zipDataFiles[1])
  # batches data file
  currentBatchesFile <- convertNulls(zipDataFiles[2])
  # file to write filtered pipeline matrix
  pipelineMatrixFile <- convertNulls(zipDataFiles[3])
  pipelineLogTransformFile <- file.path(dirname(pipelineMatrixFile), "logTransform.txt")
  # file to write filtered pipeline batches
  pipelineBatchesFile <- convertNulls(zipDataFiles[4])
  # matrix data file
  currentMatrixFile2 <- convertNulls(zipDataFiles[5])
  # batches data file
  currentBatchesFile2 <- convertNulls(zipDataFiles[6])
  # file to write filtered pipeline matrix
  pipelineMatrixFile2 <- convertNulls(zipDataFiles[7])
  pipelineLogTransformFile2 <- file.path(dirname(pipelineMatrixFile2), "logTransform2.txt")
  # file to write filtered pipeline batches
  pipelineBatchesFile2 <- convertNulls(zipDataFiles[8])
  # current original and pipeline directories
  zipDataCurrentOrigDir <- convertNulls(zipDataFiles[9])
  zipDataCurrentPipeDir <- convertNulls(zipDataFiles[10])
  # directory to be created and populated with original and pipeline data files
  zipDataVersionDir <- convertNulls(zipDataFiles[11])
  ########################
  # TODO: add correction and batch type to setupZipResultsDir call
  zipResultsFiles <- setupZipResultsDir(dirname(theConfigFile), theZipResultsDir, myDataVersion,
                                        myTestVersion, myConfig["selectedBatchToCorrect",],
                                        myConfig["selectedCorrection",],
                                        myConfig["selectedCorrectionReason",])
  # directory for logging and other run specific files
  loggingDir <- zipResultsFiles[1]
  # directory for correction files (path includes versions) not created yet
  correctionDir <- zipResultsFiles[2]
  # directory for analysis output
  analysisOutputDir <- zipResultsFiles[3]
  # config file
  versionedConfigFile <- zipResultsFiles[4]
  # info/<version>/title.txt (to be written)
  versionedTitleFile <- zipResultsFiles[5]
  # reasons dataframe
  reasonsDF <- data.frame (
    key = c("BatchType", "Correction", "Reason"),
    value = c(myConfig["selectedBatchToCorrect",],
              myConfig["selectedCorrection",],
              myConfig["selectedCorrectionReason",]))
  ########################
  # setup results directories
  myJobId <- convertNulls(myConfig["jobId",])
  message("Populate needed files with default values")
  neededFileDir <- dirname(versionedConfigFile)
  populateNeededFiles(neededFileDir, myJobId)
  ####################################################################
  rbnOnly <- as.logical(convertNulls(myConfig["RBN_Only",]))
  mutationsMutbatchFlag <- as.logical(convertNulls(myConfig["mutationsMutbatchFlag",]))
  if (isTRUE(rbnOnly))
  {
    #
    # RBN only
    #
    runRBNfromConfig(versionedConfigFile, correctionDir, loggingDir,
                     currentMatrixFile, currentMatrixFile2,
                     myDataVersion, myTestVersion)
  }
  else if (isTRUE(mutationsMutbatchFlag))
  {
    #
    # MutBatch Mutations MAF (multi-dataset) Flag
    #
    batchTypesForMBatch <- myConfig["batchTypesForMBatchArray",]
    if (!is.null(batchTypesForMBatch))
    {
      batchTypesForMBatch <- strsplit(batchTypesForMBatch, ",")[[1]]
    }
    mutBatchMem <- myConfig["mutBatchMem",]
    mutBatchThreads <- as.integer(myConfig["mutBatchThreads",])
    mutBatchPvalueCutoff <- as.numeric(myConfig["mutBatchPvalueCutoff",])
    mutBatchZscoreCutoff <- as.numeric(myConfig["mutBatchZscoreCutoff",])
    mutBatchDataBaseDir <- myConfig["mutBatchDataBaseDir",]
    mutBatchExtractDir <- myConfig["mutBatchExtractDir",]
    mutBatchOutputDir <- myConfig["mutBatchOutputDir",]
    message("mutBatchMem ", mutBatchMem)
    message("batchTypesForMBatch ", batchTypesForMBatch)
    message("mutBatchThreads ", mutBatchThreads)
    message("mutBatchPvalueCutoff ", mutBatchPvalueCutoff)
    message("mutBatchZscoreCutoff ", mutBatchZscoreCutoff)
    message("mutBatchDataBaseDir ", mutBatchDataBaseDir)
    message("mutBatchExtractDir ", mutBatchExtractDir)
    message("mutBatchOutputDir ", mutBatchOutputDir)
    mutationBatchExtract(theMafDir=mutBatchDataBaseDir, theTypeCountDir=mutBatchExtractDir, theGDCflag=TRUE)
    mutationBatchAssess(theTypeCountDir=mutBatchExtractDir,
                        theOutputDir=mutBatchOutputDir,
                        theJavaArgs=mutBatchMem,
                        theThreads=mutBatchThreads,
                        thePvalueCutoff=mutBatchPvalueCutoff,
                        theZScoreCutoff=mutBatchZscoreCutoff,
                        thePCAflag=FALSE,
                        theBatchTypes=batchTypesForMBatch)
  }
  else
  {
    #
    # General MBatch Run
    #
    title <- myConfig["Title",]
    if ((is.na(title)) || ("" == title))
    {
      message("Title entry empty - check for original.json file")
      jsonFile <- file.path(dirname(theConfigFile), "original_data.json")
      if (file.exists(jsonFile))
      {
        message("found json file ", jsonFile)
        json <- read_json(jsonFile)
        dataVersion <- myConfig["DataVersion",]
        testVersion <- myConfig["TestVersion",]
        title <- paste(json$source, " / ",
                       json$program, " / ",
                       json$project, " / ",
                       json$category, " / ",
                       json$platform, " / ",
                       json$data,
                       if (""!=json$details) {paste(" / ", json$details, sep="")} else {""},
                       " / ", dataVersion,
                       " / ", testVersion, sep="")
      }
      else
      {
        message("no json ", jsonFile)
        title <- ""
      }
    }
    titleFile <- file.path(dirname(theConfigFile), "title.txt")
    message("write versionedTitleFile ", versionedTitleFile)
    writeLines(title, versionedTitleFile)
    sampleidBatchType <- myConfig["sampleidBatchType",]
    ngchmRowType <- myConfig["ngchmRowType",]
    ngchmColumnType <- myConfig["ngchmColumnType",]
    if (is.null(ngchmRowType))
    {
      ngchmRowType <- "scholar"
    }
    if (is.null(ngchmColumnType))
    {
      ngchmColumnType <- "bio.tcga.barcode.sample"
    }
    batchTypesForMBatch <- myConfig["batchTypesForMBatchArray",]
    if (!is.null(batchTypesForMBatch))
    {
      batchTypesForMBatch <- strsplit(batchTypesForMBatch, ",")[[1]]
    }
    filterMaxValue <- myConfig["filterMaxValue",]
    filterLogTransformFlag <- myConfig["filterLogTransformFlag",]
    filterLogTransformFlag2 <- myConfig["filterLogTransformFlag2",]
    filteringBatchType <- myConfig["filteringBatchType",]
    filteringBatches <- convertStringArray(convertNulls(myConfig["filteringBatchesArray",]))
    CDP_Flag <- as.logical(convertNulls(myConfig["CDP_Flag",]))
    message("title ", title)
    message("sampleidBatchType ", sampleidBatchType)
    message("batchTypesForMBatch ", paste(batchTypesForMBatch, sep = " = "))
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
    RBN_InvariantReps <- convertStringArray(convertNulls(myConfig["RBN_InvariantRepsArray",]))
    RBN_VariantReps <- convertStringArray(convertNulls(myConfig["RBN_VariantRepsArray",]))
    message("RBN_UseFirstAsInvariant ", RBN_UseFirstAsInvariant)
    message("RBN_InvariantId ", RBN_InvariantId)
    message("RBN_VariantId ", RBN_VariantId)
    message("RBN_Matched ", RBN_Matched)
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
    # MutBatch run settings
    # MutBatch single-dataset is now part of normal run
    mutBatchMem <- myConfig["mutBatchMem",]
    mutBatchThreads <- as.integer(myConfig["mutBatchThreads",])
    mutBatchPvalueCutoff <- as.numeric(myConfig["mutBatchPvalueCutoff",])
    mutBatchZscoreCutoff <- as.numeric(myConfig["mutBatchZscoreCutoff",])
    message("mutBatchMem ", mutBatchMem)
    message("mutBatchThreads ", mutBatchThreads)
    message("mutBatchPvalueCutoff ", mutBatchPvalueCutoff)
    message("mutBatchZscoreCutoff ", mutBatchZscoreCutoff)
    runUpdateOnly <- as.logical(convertNulls(convertNA(myConfig["runUpdateOnly",])))
    message("runUpdateOnly ", runUpdateOnly)
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
    filterMaxValue <- as.numeric(filterMaxValue)
    filterLogTransformFlag <- as.logical(convertNulls(convertNA(filterLogTransformFlag)))
    theLogFrameFlag <- FALSE
    if (filterLogTransformFlag)
    {
      theLogFrameFlag <- TRUE
    }
    filterLogTransformFlag2 <- as.logical(convertNulls(convertNA(filterLogTransformFlag2)))
    selectedCorrectionMinBatchSize <- as.numeric(selectedCorrectionMinBatchSize)
    selectedDSCPermutations <- as.numeric(selectedDSCPermutations)
    selectedDSCThreads <- as.numeric(selectedDSCThreads)
    selectedDSCMinBatchSize <- as.numeric(selectedDSCMinBatchSize)
    selectedDSCSeed <- as.numeric(selectedDSCSeed)
    selectedDSCMaxGeneCount <- as.numeric(selectedDSCMaxGeneCount)
    selectedBoxplotMaxGeneCount <- as.numeric(selectedBoxplotMaxGeneCount)
    ####################################################################
    ngchmFeatureMapFile <- NULL
    nmapFile <- cleanFilePath(loggingDir, "ngchm_link_map.tsv")
    if (file.exists(nmapFile))
    {
      ngchmFeatureMapFile <- nmapFile
    }
    #############################################################################
    # check directories
    message("analysisOutputDir=", analysisOutputDir)
    if(!dir.exists(analysisOutputDir))
    {
      message("create ", analysisOutputDir)
      dir.create(analysisOutputDir, showWarnings=FALSE, recursive=TRUE)
    }
    #############################################################################
    logFile <- cleanFilePath(loggingDir, "mbatch.log")
    setLogging(new("Logging", theFile=logFile))
    #############################################################################
    message("update data to use sample")
    batchesFileToUse <- currentBatchesFile
    batchesFileToUse2 <- currentBatchesFile2
    {
      # rename sample id column to Sample
      if("Sample"!=sampleidBatchType)
      {
        batchs <- readAsGenericDataframe(currentBatchesFile)
        names(batchs)[names(batchs)==sampleidBatchType] <- "Sample"
        batchesWithSampleColumn <- cleanFilePath(loggingDir, "BatchesWithSampleCol.tsv")
        message("copy ", currentBatchesFile, " to ", batchesWithSampleColumn)
        file.copy(currentBatchesFile, batchesWithSampleColumn, overwrite=TRUE)
        writeAsGenericDataframe(batchesWithSampleColumn, batchs)
        batchesFileToUse <- batchesWithSampleColumn
        if(!is.null(currentBatchesFile2))
        {
          batchs <- readAsGenericDataframe(currentBatchesFile2)
          names(batchs)[names(batchs)==sampleidBatchType] <- "Sample"
          batchesWithSampleColumn2 <- paste(loggingDir, "BatchesNewSampleCol2.tsv", sep="")
          message("copy ", currentBatchesFile2, " to ", batchesWithSampleColumn2)
          file.copy(currentBatchesFile2, batchesWithSampleColumn2, overwrite=TRUE)
          writeAsGenericDataframe(batchesWithSampleColumn2, batchs)
          batchesFileToUse2 <- batchesWithSampleColumn2
        }
      }
    }
    message("preprocess data")
    filterList <- list()
    if (!is.null(filteringBatchType))
    {
      filterList <- lapply(filteringBatches, function(x) { c(filteringBatchType, x)})
    }
    # also writes new matrix.tsv and batches.tsv files
    preprocessData(currentMatrixFile, batchesFileToUse,
                   pipelineMatrixFile, pipelineBatchesFile,
                   theSize=0, theTransformFlag=filterLogTransformFlag,
                   theBatchTypeAndValuePairsToRemove=filterList,
                   theLogTransformFile=pipelineLogTransformFile,
                   theReplaceNAs=replaceNAs)
    if(!is.null(currentMatrixFile2))
    {
      preprocessData(currentMatrixFile2, batchesFileToUse2,
                     pipelineMatrixFile2, pipelineBatchesFile2,
                     theSize=0, theTransformFlag=filterLogTransformFlag2,
                     theLogTransformFile=pipelineLogTransformFile2,
                     theReplaceNAs=replaceNAs)
    }
    #############################################################################
    #### do data load
    #############################################################################
    message("load data as needed")
    myOriginalData <- NULL
    myOriginalData2 <- NULL
    myMBatchData <- NULL
    myMBatchData2 <- NULL
    if ((!is.null(selectedCorrection)) && ( startsWith(selectedCorrection, "RBN_") ))
    {
      message("Read RBN data")
      if (isTRUE(RBN_UseFirstAsInvariant))
      {
        # TODO: use , theNaStrings=theNaStrings after changing pre-set values
        myMBatchData <- mbatchLoadFiles(pipelineMatrixFile, pipelineBatchesFile)
        myMBatchData2 <- mbatchLoadFiles(pipelineMatrixFile2, pipelineBatchesFile2)
      }
      else
      {
        # TODO: use , theNaStrings=theNaStrings after changing pre-set values
        myMBatchData <- mbatchLoadFiles(pipelineMatrixFile, pipelineBatchesFile)
        myMBatchData2 <- mbatchLoadFiles(pipelineMatrixFile2, pipelineBatchesFile2)
      }
    }
    else if ((!is.null(selectedCorrection)) && ( startsWith(selectedCorrection, "EBN_") ))
    {
      message("Read EBNplus data")
      # EBNplus cannot use the mbatchLoadFiles function, since it filters meaningfule data
      matrix1 <- readAsGenericMatrix(pipelineMatrixFile)
      matrix2 <- readAsGenericMatrix(pipelineMatrixFile2)
      df1 <- readAsGenericDataframe(pipelineBatchesFile)
      df2 <- readAsGenericDataframe(pipelineBatchesFile2)
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
      message(pipelineMatrixFile)
      message(pipelineBatchesFile)
      # TODO: use , theNaStrings=theNaStrings after changing pre-set values
      myMBatchData <- mbatchLoadFiles(pipelineMatrixFile, pipelineBatchesFile)
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
    # reduce size as desired
    # myMBatchData@mData <- mbatchTrimData(myMBatchData@mData, filterMaxValue)
    # if (!is.null(myMBatchData2))
    # {
    #   myMBatchData2@mData <- mbatchTrimData(myMBatchData2@mData, filterMaxValue)
    # }
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
      myMBatchData <- doCorrectionsFromConfig(theOutputDir=correctionDir,
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
                                              theEBNPlus_MinSamples=EBNPlus_MinSamples,
                                              theDataVersion=myDataVersion,
                                              theTestVersion=myTestVersion,
                                              theReasons=reasonsDF)
    }
    if ((!is.null(selectedCorrection))&&(is.null(myMBatchData)))
    {
      # correction not performed, one batch
      # set flag for other than full success
      otherNote <- TRUE
      #############################################################################
      message("Write completed note to '", cleanFilePath(analysisOutputDir, "MBATCH_COMPLETED.txt"), "'")
      file.create(cleanFilePath(analysisOutputDir, "MBATCH_COMPLETED.txt"))
    }
    else
    {
      #### do filtering
      # do not use sampleidBatchType, since that column has been changed to "Sample"
      #batchesToRemove <- setdiff(colnames(myMBatchData@mBatches), c(sampleidBatchType, batchTypesForMBatch))
      batchesToRemove <- setdiff(colnames(myMBatchData@mBatches), c("Sample", batchTypesForMBatch))
      myMBatchData <- mbatchFilterData(theBeaData=myMBatchData,
                                       theBatchTypeAndValuePairsToRemove=list(),
                                       theBatchTypeAndValuePairsToKeep=list(),
                                       theBatchTypesToRemove=batchesToRemove,
                                       theMinIqr=0,
                                       theMinSd=0,
                                       theMinMad=0)
      if (isTRUE(CDP_Flag))
      {
        # do not use sampleidBatchType, since that column has been changed to "Sample"
        #batchesToRemove <- setdiff(colnames(myOriginalData@mBatches), c(sampleidBatchType, batchTypesForMBatch))
        batchesToRemove <- setdiff(colnames(myOriginalData@mBatches), c("Sample", batchTypesForMBatch))
        myOriginalData <- mbatchFilterData(theBeaData=myOriginalData,
                                         theBatchTypeAndValuePairsToRemove=list(),
                                         theBatchTypeAndValuePairsToKeep=list(),
                                         theBatchTypesToRemove=batchesToRemove,
                                         theMinIqr=0,
                                         theMinSd=0,
                                         theMinMad=0)
      }
      if (isTRUE(runUpdateOnly))
      {
        # DO UPDATE ONLY
        updateFromConfig(theOutputDir=analysisOutputDir,
                         theDataVersion=myDataVersion,
                         theTestVersion = myTestVersion,
                         theDataObject = myMBatchData,
                         theTitle = title,
                         thePermutations = selectedDSCPermutations,
                         thePermutationThreads = selectedDSCThreads,
                         theMaxSize = filterMaxValue,
                         theSeed = selectedDSCSeed,
                         theDSCMaxGenes = selectedDSCMaxGeneCount,
                         theBoxplotMaxGenes = selectedBoxplotMaxGeneCount,
                         theMinBatchSize = selectedDSCMinBatchSize,
                         theNgchmRowType = ngchmRowType,
                         theNgchmColumnType = ngchmColumnType,
                         theShaidyMapGen = theShaidyMapGen,
                         theNgchmWidgetJs = theNgchmWidgetJs,
                         theShaidyMapGenJava = theShaidyMapGenJava,
                         theNGCHMShaidyMem = theNGCHMShaidyMem,
                         theNgchmFeatureMapFile = ngchmFeatureMapFile,
                         theLogFrameFlag = theLogFrameFlag)
      }
      else
      {
        #
        # MutBatch single-dataset is now part of normal run
        #
        mutBatchSingle(myMBatchData, title, analysisOutputDir, myDataVersion, myTestVersion,
                       theThreads=mutBatchThreads,
                       thePvalueCutoff=mutBatchPvalueCutoff,
                       theZScoreCutoff=mutBatchZscoreCutoff,
                       theBatchTypes=batchTypesForMBatch)
        #### do assessments
        # do not use sampleidBatchType, since that column has been changed to "Sample"
        rdataHCvector <- doAssessmentsFromConfig(
          theOutputDir = analysisOutputDir,
          theDataVersion = myDataVersion,
          theTestVersion = myTestVersion,
          theDataObject = myMBatchData,
          theTitle = title,
          thePermutations = selectedDSCPermutations,
          thePermutationThreads = selectedDSCThreads,
          theMaxSize = filterMaxValue,
          theSeed = selectedDSCSeed,
          theDSCMaxGenes = selectedDSCMaxGeneCount,
          theBoxplotMaxGenes = selectedBoxplotMaxGeneCount,
          theMinBatchSize = selectedDSCMinBatchSize,
          theNgchmRowType = ngchmRowType,
          theNgchmColumnType = ngchmColumnType,
          theShaidyMapGen = theShaidyMapGen,
          theNgchmWidgetJs = theNgchmWidgetJs,
          theShaidyMapGenJava = theShaidyMapGenJava,
          theNGCHMShaidyMem = theNGCHMShaidyMem,
          theNgchmFeatureMapFile = ngchmFeatureMapFile,
          theLogFrameFlag = theLogFrameFlag)
        message("rdataHCvector=", paste(rdataHCvector, sep=", ", collase=", "))
        rdataHCsample <- rdataHCvector[[1]]
        rdataHCfeature <- rdataHCvector[[2]]
        message("rdataHCsample=", rdataHCsample)
        message("rdataHCfeature=", rdataHCfeature)
        if(isTRUE(theRunPostFlag))
        {
          pcaOutputDir <- cleanFilePath(analysisOutputDir, "PCA")
          dscOutputDir <- cleanFilePath(analysisOutputDir, "DSC")
          dscOutputDir <- addVersionsIfNeeded(dscOutputDir, myDataVersion, myTestVersion)
          message("buildDSCOverviewFile analysisOutputDir = ", analysisOutputDir)
          message("buildDSCOverviewFile dscOutputDir = ", dscOutputDir)
          buildDSCOverviewFile(pcaOutputDir, dscOutputDir, analysisOutputDir, "DSCOverview.tsv")
        }
        #############################################################################
        if (isTRUE(selectedNgchmFlag))
        {
          warnLevel<-getOption("warn")
          on.exit(options(warn=warnLevel))
          options(warn=-1)
          myBatchType <- "All"
          message("title ", title)
          message("myBatchType ", myBatchType)
          message("theShaidyMapGen ", theShaidyMapGen)
          message("theShaidyMapGenJava ", theShaidyMapGenJava)
          message("dim(myMBatchData@mData) ", dim(myMBatchData@mData))
          message("dim(myMBatchData@mBatches) ", dim(myMBatchData@mBatches))
          message("trim to same size as hierarchical 1")
          myMBatchData <- as.numericWithIssues(myMBatchData)
          ngchmData <- mbatchTrimData(myMBatchData@mData, theMaxSize = (selectedBoxplotMaxGeneCount * ncol(myMBatchData@mData)))
          buildBatchHeatMapFromHC_Structures(theMatrixData=ngchmData,
                                       theBatchData=myMBatchData@mBatches,
                                       theTitle=paste(title, "/", myBatchType, "NGCHM", sep=" "),
                                       theOutputDir=cleanFilePath(analysisOutputDir, "NGCHM"),
                                       theDataVersion=myDataVersion,
                                       theTestVersion=myTestVersion,
                                       theOutputFilename=paste(myBatchType, "ngchm.ngchm", sep="_"),
                                       theColDendRDataFile=rdataHCsample,
                                       theRowDendRDataFile=NULL,
                                       theNgchmFeatureMapFile=ngchmFeatureMapFile,
                                       theRowType=ngchmRowType, theColType=ngchmColumnType,
                                       theRowCluster=c("pearson", "ward.D2"),
                                       theShaidyMapGen=theShaidyMapGen,
                                       theNgchmWidgetJs=theNgchmWidgetJs,
                                       theShaidyMapGenJava=theShaidyMapGenJava,
                                       theShaidyMapGenArgs=c(paste(c("-Xms", "-Xmx"), theNGCHMShaidyMem, sep=""), "-Djava.awt.headless=true"))
        }
        #############################################################################
        if (isTRUE(CDP_Flag))
        {
          graphics.off()
          warnLevel<-getOption("warn")
          on.exit(options(warn=warnLevel))
          options(warn=-1)
          outdir <- addVersionsIfNeeded(cleanFilePath(analysisOutputDir, "CDP"), myDataVersion, myTestVersion)
          if (!file.exists(outdir))
          {
            dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
          }
          CDP_Structures(cleanFilePath(analysisOutputDir, "CDP"), myDataVersion,
                         myTestVersion, "CDP_Plot_Data1_Diagram.PNG", myOriginalData@mData,
                         myMBatchData@mData, "Primary-Dataset", theTitle=title)
          if (!is.null(myOriginalData2))
          {
            CDP_Structures(cleanFilePath(analysisOutputDir, "CDP"), myDataVersion,
                           myTestVersion, "CDP_Plot_Data2_Diagram.PNG", myOriginalData2@mData,
                           myMBatchData@mData, "Secondary-Dataset", theTitle=title)
          }
        }
      }
    }
  }
  if (isFALSE(otherNote))
  {
    #############################################################################
    message("Write success note to '", cleanFilePath(analysisOutputDir, "MBATCH_SUCCESS.txt"), "'")
    file.create(cleanFilePath(analysisOutputDir, "MBATCH_SUCCESS.txt"))
  }
  # copy current ZIP-DATA to Versioned ZIP-DATA
  # current original and pipeline directories
  # zipDataCurrentOrigDir <- convertNulls(zipDataFiles[9])
  # zipDataCurrentPipeDir <- convertNulls(zipDataFiles[10])
  # directory to be created and populated with original and pipeline data files
  # zipDataVersionDir <- convertNulls(zipDataFiles[10])
  dir.create(zipDataVersionDir, recursive=TRUE, showWarnings=FALSE)
  message("Copy from ", zipDataCurrentOrigDir, " to ", zipDataVersionDir)
  file.copy(zipDataCurrentOrigDir, zipDataVersionDir, recursive=TRUE)
  message("Copy from ", zipDataCurrentPipeDir, " to ", zipDataVersionDir)
  file.copy(zipDataCurrentPipeDir, zipDataVersionDir, recursive=TRUE)
  # copy theZipDataDir/original/mutations.tsv (if it exists) to zipDataVersionDir/original/mutations.tsv
  sorcMutations <- cleanFilePath(cleanFilePath(theZipDataDir, "original"), "mutations.tsv")
  if (file.exists(sorcMutations))
  {
    destMutations <- cleanFilePath(cleanFilePath(zipDataVersionDir, "original"), "mutations.tsv")
    message("Copy mutations.tsv from ", sorcMutations, " to ", destMutations)
    file.copy(sorcMutations, destMutations)
  }
}

updateFromConfig <- function(theOutputDir, theDataVersion, theTestVersion,
                             theDataObject, theTitle,
                             thePermutations, thePermutationThreads, theMaxSize,
                             theSeed, theDSCMaxGenes, theBoxplotMaxGenes,
                             theMinBatchSize,
                             theNgchmRowType,
                             theNgchmColumnType,
                             theShaidyMapGen,
                             theNgchmWidgetJs,
                             theShaidyMapGenJava, theNGCHMShaidyMem,
                             theNgchmFeatureMapFile, theLogFrameFlag)
{
  if (theBoxplotMaxGenes>5000)
  {
    message("Limit boxplot to at most 5000 genes");
    theBoxplotMaxGenes <- 5000
  }
  # reduce size
  theDataObject@mData <- mbatchTrimData(theDataObject@mData, theMaxSize)
  ####
  ####
  #### Volcano Plot
  ###
  message("call callMBatch_VolcanoPlot_Structures")
  print("call callMBatch_VolcanoPlot_Structures")
  callMBatch_VolcanoPlot_Structures(theOutputDir, theDataVersion, theTestVersion, theDataObject,
                                    theTitle, theLogFrameFlag, theBoxplotMaxGenes)
  ####
  #### write success
  ###
  message("Processing success")
  mbatchWriteSuccessfulLog()
  list()
}


doAssessmentsFromConfig <- function(theOutputDir, theDataVersion, theTestVersion,
                                    theDataObject, theTitle,
                                    thePermutations, thePermutationThreads, theMaxSize,
                                    theSeed, theDSCMaxGenes, theBoxplotMaxGenes,
                                    theMinBatchSize,
                                    theNgchmRowType,
                                    theNgchmColumnType,
                                    theShaidyMapGen,
                                    theNgchmWidgetJs,
                                    theShaidyMapGenJava, theNGCHMShaidyMem,
                                    theNgchmFeatureMapFile, theLogFrameFlag)
{
  if (theBoxplotMaxGenes>5000)
  {
    message("Limit boxplot to at most 5000 genes");
    theBoxplotMaxGenes <- 5000
  }
  # reduce size
  theDataObject@mData <- mbatchTrimData(theDataObject@mData, theMaxSize)
  ####
  #### UMAP
  ####
  callMBatch_UMAP_Structures(theOutputDir, theDataVersion, theTestVersion, theDataObject, theTitle)
  ####
  #### PCAPlus
  ####
  callMBatch_PCA_Structures(theOutputDir, theDataVersion, theTestVersion, theDataObject, theTitle,
                            theIsPcaTrendFunction=isTrendBatch,
                            theDSCPermutations=thePermutations,
                            theDSCThreads=thePermutationThreads,
                            theMinBatchSize=theMinBatchSize,
                            theSeed=theSeed,
                            theMaxGeneCount=theDSCMaxGenes)
  ####
  #### HierarchicalClustering
  ####
  rdataHCVector <- callMBatch_HierarchicalClustering_Structures(theOutputDir, theDataVersion, theTestVersion,
                                                          theDataObject, theTitle, theBoxplotMaxGenes)
  rdataHCsample <- rdataHCVector[[1]]
  rdataHCfeature <- rdataHCVector[[2]]
  ####
  #### SupervisedClustering
  ####
  message("call callMBatch_SupervisedClustering_Structures")
  {
    superMatrix <- mbatchTrimData(theDataObject@mData, theMaxSize = (theBoxplotMaxGenes * ncol(theDataObject@mData)))
    mySuperData<-new("BEA_DATA", superMatrix, theDataObject@mBatches, data.frame())
    callMBatch_SupervisedClustering_Structures(theOutputDir, theDataVersion, theTestVersion, mySuperData, theTitle,
                                               theShaidyMapGen, theNgchmWidgetJs, theShaidyMapGenJava, theNGCHMShaidyMem,
                                               theNgchmRowType,
                                               theNgchmColumnType,
                                               theNgchmFeatureMapFile)
  }
  ####
  #### BoxPlot
  ###
  pipelineMean <- function(x)
  {
    mean(x, na.rm=TRUE)
  }
  message("call callMBatch_BoxplotAllSamplesData_Structures")
  print("call callMBatch_BoxplotAllSamplesData_Structures")
  callMBatch_BoxplotAllSamplesData_Structures(theOutputDir, theDataVersion, theTestVersion, theDataObject, theTitle,
                                              theMaxGeneCount=theBoxplotMaxGenes)
  message("call callMBatch_BoxplotAllSamplesRLE_Structures")
  print("call callMBatch_BoxplotAllSamplesRLE_Structures")
  callMBatch_BoxplotAllSamplesRLE_Structures(theOutputDir, theDataVersion, theTestVersion, theDataObject, theTitle,
                                             theMaxGeneCount=theBoxplotMaxGenes)
  message("call callMBatch_BoxplotGroup_Structures")
  print("call callMBatch_BoxplotGroup_Structures")
  callMBatch_BoxplotGroup_Structures(theOutputDir, theDataVersion, theTestVersion, theDataObject, theTitle,
                                     theMaxGeneCount=theBoxplotMaxGenes,
                                     theFunction=list(pipelineMean), theFunctionName=list("Mean"))
  ####
  #### Volcano Plot
  ###
  message("call callMBatch_VolcanoPlot_Structures")
  print("call callMBatch_VolcanoPlot_Structures")
  callMBatch_VolcanoPlot_Structures(theOutputDir, theDataVersion, theTestVersion, theDataObject,
                                    theTitle, theLogFrameFlag, theBoxplotMaxGenes)
  ####
  #### write success
  ###
  message("Processing success")
  mbatchWriteSuccessfulLog()
  list(rdataHCsample, rdataHCfeature)
}

doCorrectionsFromConfig <- function(theOutputDir, theDataObject, theDataObject2,
                                    theSelectedBatchToCorrect, theSelectedCorrection,
                                    thePermutationThreads, theNaStrings,
                                    theRBN_InvariantId, theRBN_VariantId, theRBN_Matched,
                                    theRBN_InvariantReps, theRBN_VariantReps,
                                    theEBNPlus_GroupId1, theEBNPlus_GroupId2,
                                    theEBNPlus_Seed, theEBNPlus_MinSamples,
                                    theDataVersion, theTestVersion, theReasons)
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
    writeAsGenericDataframe(theFile=cleanFilePath(theOutputDir, "corrected_batches.tsv"), theDataframe=dataBatches)
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
    writeAsGenericDataframe(theFile=cleanFilePath(theOutputDir, "corrected_batches.tsv"), theDataframe=dataBatches)
    if(!is.null(myCorrectedFile))
    {
      dataMatrix <- readAsGenericMatrix(myCorrectedFile)
    }
  }
  else if ("EBN_Plus"==theSelectedCorrection)
  {
    #theDataObject, theDataObject2, theEBNPlus_Batch1, theEBNPlus_Batch2, theEBNPlus_GroupId1, theEBNPlus_GroupId2, theEBNPlus_Seed,
    processed <- TRUE
    correctedFile <- cleanFilePath(theOutputDir, "adjusted_matrix.tsv")
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
                                                theOutputDir=theOutputDir,
                                                theDataVersion=theDataVersion,
                                                theTestVersion=theTestVersion,
                                                thePriorFile="priorplots_Diagram.PNG")
    writeAsGenericMatrix(correctedFile, dataMatrix)
    dataBatches <- EBNPlus_CombineBatches(theBeaBatches1=theDataObject@mBatches,
                                      theBeaBatches2=theDataObject2@mBatches,
                                      theEBNP_Data1BatchId=theEBNPlus_GroupId1,
                                      theEBNP_Data2BatchId=theEBNPlus_GroupId2,
                                      theBarcodeTrimFunction=trimGenes)
    writeAsGenericDataframe(theFile=cleanFilePath(theOutputDir, "corrected_batches.tsv"), theDataframe=dataBatches)
  }
  else
  {
    stop("Unrecognized correction")
  }
  if ((!is.null(dataMatrix))&&(TRUE==processed))
  {
    # if log transformed, save an untransformed version
    # ../../../ZIP-DATA/current/pipeline/logTransform.txt
    logTransformFile <- file.path(dirname(dirname(dirname(theOutputDir))), "ZIP-DATA", "current", "pipeline", "logTransform.txt")
    if (file.exists(logTransformFile))
    {
      pquantile <- read_file(logTransformFile)
      pquantile <- as.numeric(pquantile)
      unlogt <- 2^dataMatrix
      unlogt <- unlogt - pquantile
      writeAsGenericDataframe(theFile=cleanFilePath(theOutputDir, "correction_unlogtransformed.tsv"), theDataframe=unlogt)
    }
  }
  if ((!is.null(dataMatrix))&&(TRUE==processed))
  {
    # if have data info and processing was done, return new data object
    theDataObject <- mbatchLoadStructures(dataMatrix, dataBatches)
  }
  if (TRUE==processed)
  {
    writeAsGenericDataframe(theFile=cleanFilePath(theOutputDir, "adjusted_reasons.txt"), theDataframe=theReasons)
  }
  dataMatrix <- NULL
  dataBatches <- NULL
  theDataObject
}
