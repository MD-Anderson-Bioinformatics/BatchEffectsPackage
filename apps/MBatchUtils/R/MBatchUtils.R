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

mbatchUtilVersion <- function()
{
  return("MBatchUtils BEA_VERSION_TIMESTAMP")
}


####################################################################
###
####################################################################

# index file
# current
# 	original	(matrix, batches, clinical, dataset info, etc)
# 	pipeline	(same as above, but as used for processing)
# versions	(broken up by versions)
# 	DATA-YYYY-MM-DD-HHmm
# 		original	(matrix, batches, clinical, dataset info, etc)
# 		pipeline	(same as above, but as used for processing)

setupZipDataDir <- function(theMatrixFile, theBatchFile, theNewDir, theDataVersion,
                            theMatrixFile2=NULL, theBatchFile2=NULL)
{
  message("theMatrixFile=",theMatrixFile)
  message("theBatchFile=",theBatchFile)
  message("theNewDir=",theNewDir)
  message("theDataVersion=",theDataVersion)
  message("theMatrixFile2=",theMatrixFile2)
  message("theBatchFile2=",theBatchFile2)
  # make directories
  currentDir <- cleanFilePath(theNewDir, "current")
  versionsDir <- cleanFilePath(theNewDir, "versions")
  message("dir.create ",currentDir)
  dir.create(currentDir, recursive=TRUE)
  message("dir.create ",versionsDir)
  dir.create(versionsDir, recursive=TRUE)
  originalDir <- cleanFilePath(currentDir, "original")
  pipelineDir <- cleanFilePath(currentDir, "pipeline")
  message("unlink ",originalDir)
  unlink(originalDir, recursive=TRUE)
  message("unlink ",pipelineDir)
  unlink(pipelineDir, recursive=TRUE)
  message("dir.create ",originalDir)
  dir.create(originalDir, recursive=TRUE)
  message("dir.create ",pipelineDir)
  dir.create(pipelineDir, recursive=TRUE)
  dataVersDir <- addVersionsIfNeeded(versionsDir, theDataVersion, NULL)
  message("dir.create ",dataVersDir)
  dir.create(dataVersDir, recursive=TRUE)
  # copy files
  currentMatrixFile <- file.path(originalDir, "matrix.tsv")
  message("file.copy ",theMatrixFile, " to ", currentMatrixFile)
  file.copy(theMatrixFile, currentMatrixFile)
  currentBatchesFile <- ""
  if (!is.null(theBatchFile))
  {
    currentBatchesFile <- file.path(originalDir, "batches.tsv")
    message("file.copy ",theBatchFile, " to ", currentBatchesFile)
    file.copy(theBatchFile, currentBatchesFile)
  }
  currentMatrixFile2 <- ""
  if (!is.null(theMatrixFile2))
  {
    currentMatrixFile2 <- file.path(originalDir, "matrix2.tsv")
    message("file.copy ",theMatrixFile2, " to ", currentMatrixFile2)
    file.copy(theMatrixFile2, currentMatrixFile2)
  }
  currentBatchesFile2 <- ""
  if (!is.null(theBatchFile2))
  {
    currentBatchesFile2 <- file.path(originalDir, "batches2.tsv")
    message("file.copy ",theBatchFile2, " to ", currentBatchesFile2)
    file.copy(theBatchFile2, currentBatchesFile2)
  }
  # pipeline versions
  pipelineMatrixFile <- file.path(pipelineDir, "matrix.tsv")
  pipelineBatchesFile <- file.path(pipelineDir, "batches.tsv")
  pipelineMatrixFile2 <- file.path(pipelineDir, "matrix2.tsv")
  pipelineBatchesFile2 <- file.path(pipelineDir, "batches2.tsv")
  message("pipeline versions pipelineMatrixFile=", pipelineMatrixFile)
  message("pipeline versions pipelineBatchesFile=", pipelineBatchesFile)
  c(currentMatrixFile, currentBatchesFile,
    pipelineMatrixFile, pipelineBatchesFile,
    currentMatrixFile2, currentBatchesFile2,
    pipelineMatrixFile2, pipelineBatchesFile2,
    originalDir, pipelineDir, dataVersDir)
}

# Results Archive Directory Structure
# index file
# info
# 	non-version dataset information file(s)
# 	TEST-YYYY-MM-DD-HHmm
# 		MBatchConfig used for processing by version
# 		version-specific dataset info and log file(s)
# correction
# 	DATA-YYYY-MM-DD-HHmm
# 		TEST-YYYY-MM-DD-HHmm
# 			correction results files (info, matrix, and batches)
# analysis
# 	<analysis-type>
# 		<sub directories>	(depends on analysis type)
# 			DATA-YYYY-MM-DD-HHmm
# 				TEST-YYYY-MM-DD-HHmm	(broken up by versions)
# 					static	(static analysis results – png files)
# 					dynamic	(dynamic analysis results – tsv files for JavaScript/D3)

setupZipResultsDir <- function(theConfigDir, theNewDir, theDataVersion, theTestVersion,
                               theBatchType, theCorrectionType, theCorrectionReason)
{
  # make directories
  infoDir <- cleanFilePath(theNewDir, "info")
  analysisDir <- cleanFilePath(theNewDir, "analysis")
  dir.create(infoDir, recursive=TRUE, showWarnings=FALSE)
  dir.create(analysisDir, recursive=TRUE, showWarnings=FALSE)
  infoTestVerDir <- addVersionsIfNeeded(infoDir, NULL, theTestVersion)
  dir.create(infoTestVerDir, recursive=TRUE, showWarnings=FALSE)
  # copy MBatchConfig.tsv
  versionedConfigFile <- file.path(infoTestVerDir, "MBatchConfig.tsv")
  file.copy(file.path(theConfigDir, "MBatchConfig.tsv"), versionedConfigFile)
  # copy other files
  file.copy(file.path(theConfigDir, "version_type.txt"), file.path(infoTestVerDir, "version_type.txt"))
  file.copy(file.path(theConfigDir, "version_stamp.txt"), file.path(infoTestVerDir, "version_stamp.txt"))
  file.copy(file.path(theConfigDir, "source_id.txt"), file.path(infoTestVerDir, "source_id.txt"))
  file.copy(file.path(theConfigDir, "original_data.json"), file.path(infoTestVerDir, "original_data.json"))
  # optional directory path for corrections
  correctionDir <- cleanFilePath(theNewDir, "correction")
  correctionDir <- addVersionsIfNeeded(correctionDir, NULL, theTestVersion)
  # then add combined correction and batch type to correction directory path
  correctionDir <- cleanFilePath(correctionDir, paste(theCorrectionType, "-",
                                                      cleanAndShortPath(theBatchType), "-",
                                                      theCorrectionReason, sep="", collapse=""))
  c(infoTestVerDir, correctionDir, analysisDir, versionedConfigFile, file.path(infoTestVerDir, "title.txt"))
}

####################################################################
###
####################################################################

getJava <- function()
{
  baseDir <- Sys.getenv("MBATCHUTILS_TEST_JAVA")
  if(!file.exists(baseDir))
  {
    baseDir <- "/software/x86_64/JDK/1.8.0_25/bin/java"
    if(!file.exists(baseDir))
    {
      baseDir <- "/usr/bin/java"
      if(!file.exists(baseDir))
      {
        if(!file.exists(baseDir))
        {
          baseDir <- NULL
        }
      }
    }
  }
  baseDir
}

getTestInputDir <- function()
{
  baseDir <- Sys.getenv("MBATCHUTILS_TEST_INPUT")
  if(!file.exists(baseDir))
  {
    baseDir <- "/BEA/BatchEffectsPackage_data/testing_static/MBatchUtils"
  }
  baseDir
}

getTestInputDirForMBatch <- function()
{
  value <- Sys.getenv("MBATCH_TEST_INPUT")
  if (!isTRUE(file.exists(value)))
  {
    value <- "/BEA/BatchEffectsPackage_data/testing_static/MATRIX_DATA"
  }
  value
}

getTestOutputDir <- function()
{
  baseDir <- Sys.getenv("MBATCHUTILS_TEST_OUTPUT")
  if(!file.exists(baseDir))
  {
    baseDir <- "/BEA/BatchEffectsPackage_data/testing_dynamic/MBatchUtils"
  }
  baseDir
}

####################################################################
###
####################################################################

compareTwoMatrices <- function(theCorrected, theCompare)
{
  correctedRows <- dim(theCorrected)[1]
  correctedCols <- dim(theCorrected)[2]
  compareRows <- dim(theCompare)[1]
  compareCols <- dim(theCompare)[2]
  message("correctedRows=", correctedRows)
  message("correctedCols=", correctedCols)
  message("compareRows=", compareRows)
  message("compareCols=", compareCols)

  stopifnot(correctedRows==compareRows)
  stopifnot(correctedCols==compareCols)

  for(myCol in 1:correctedCols)
  {
    if (!(colnames(theCorrected)[myCol]==colnames(theCompare)[myCol]))
    {
      message("FALSE mismatch 1")
      message("myCol=", myCol)
      message("colnames(theCorrected)[myCol]=", colnames(theCorrected)[myCol])
      message("colnames(theCompare)[myCol]=", colnames(theCompare)[myCol])
      return(FALSE)
    }
    for(myRow in 1:correctedRows)
    {
      if (!(rownames(theCorrected)[myRow]==rownames(theCompare)[myRow]))
      {
        message("FALSE mismatch 2")
        message("myRow=", myRow)
        message("rownames(theCorrected)[myRow]=", rownames(theCorrected)[myRow])
        message("rownames(theCompare)[myRow]=", rownames(theCompare)[myRow])
        return(FALSE)
      }
      # if ((myRow==2762)&&(myCol==1))
      # {
      #   message("myRow=", myRow)
      #   message("myCol=", myCol)
      #   message("theCorrected[myRow, myCol]=", theCorrected[myRow, myCol])
      #   message("theCompare[myRow, myCol]=", theCompare[myRow, myCol])
      #   browser();
      # }
      #message("checking myRow=", myRow, " and myCol=", myCol, " theCorrected[myRow, myCol]=", theCorrected[myRow, myCol], " theCompare[myRow, myCol]=",theCompare[myRow, myCol])
      if (!is.finite(theCorrected[myRow, myCol]) ||
          !is.finite(theCompare[myRow, myCol]) ||
          is.na(theCorrected[myRow, myCol]) ||
          is.na(theCompare[myRow, myCol]))
      {
        # ignore
      }
      else if (!(isTRUE(all.equal(theCorrected[myRow, myCol], theCompare[myRow, myCol]))))
      {
        message("FALSE mismatch 3")
        message("myRow=", myRow)
        message("myCol=", myCol)
        message("theCorrected myRow name =", rownames(theCorrected)[myRow])
        message("theCompare myRow name   =", rownames(theCompare)[myRow])
        message("theCorrected myCol name =", colnames(theCorrected)[myCol])
        message("theCompare myCol name   =", colnames(theCompare)[myCol])
        message("theCorrected[myRow, myCol]=", theCorrected[myRow, myCol])
        message("theCompare[myRow, myCol]  =", theCompare[myRow, myCol])
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

####################################################################
###
####################################################################

removeSamplesWithListedBatches<-function(theDataframeFilteredSamplesToBatches, theBatchTypeAndValuePairsToRemove)
{
  for(myPair in theBatchTypeAndValuePairsToRemove)
  {
    theDataframeFilteredSamplesToBatches <- removeBatchAndValue(theDataframeFilteredSamplesToBatches, myPair[[1]], myPair[[2]])
  }
  return(theDataframeFilteredSamplesToBatches)
}

preprocessData <- function(theInputMatrixFile, theInputBatchFile,
                           theOutputMatrixFile, theOutputBatchFile,
                           theSize, theTransformFlag,
                           theBatchTypeAndValuePairsToRemove=NULL,
                           theLogTransformFile=NULL,
                           theReplaceNAs=FALSE)
{
  print(paste("read ", theInputMatrixFile, sep=""))
  mymatrix <- readAsGenericMatrix(theInputMatrixFile)
  print(paste("read ", theInputBatchFile, sep=""))
  batchs <- readAsGenericDataframe(theInputBatchFile)
  #############################################
  # add missing samples from matrix to batchs
  #############################################
  matrixSamples <- colnames(mymatrix)
  batchsSamples <- batchs$Sample
  missingSamples <- setdiff(matrixSamples, batchsSamples)
  print(paste("length(matrixSamples)=", length(matrixSamples), sep=""))
  print(paste("length(batchsSamples)=", length(batchsSamples), sep=""))
  print(paste("length(missingSamples)=", length(missingSamples), sep=""))
  for (newSample in missingSamples)
  {
    # assign Unknown to combined batches
    batchs <- rbind(batchs, c(newSample, rep("Unknown", length.out=length(batchsSamples)-1)))
  }
  unusedSamples <- setdiff(batchsSamples, matrixSamples)
  for (notSample in unusedSamples)
  {
    batchs <- batchs[batchs$Sample != notSample,]
  }
  #############################################
  # sort, etc the data
  #############################################
  mbatchStr <- mbatchLoadStructures(mymatrix, batchs)
  mymatrix <- mbatchStr@mData
  batchs <- mbatchStr@mBatches
  #####################################################################################################################
  ### remove any samples with given batch types and values
  if (length(theBatchTypeAndValuePairsToRemove)>0)
  {
    print("preprocessData Before removing requested batches, batch data has ", length(batchs$Sample), " samples")
    #print("preprocessData dataframe $Sample = ", paste(batchs$Sample, collapse=", "))
    print("preprocessData Before removing requested batches, gene data has ", length(colnames(mymatrix)), " samples")
    print("preprocessData Remove requested batches (", paste(theBatchTypeAndValuePairsToRemove, sep="", collapse=", "), ")")
    batchs <- removeSamplesWithListedBatches(batchs, theBatchTypeAndValuePairsToRemove)
    print("preprocessData After removing requested batches, batch data has ", length(batchs$Sample), " samples")
    #print("preprocessData dataframe $Sample = ", paste(batchs$Sample, collapse=", "))
    mymatrix <- setGenesToSameAsSamples(mymatrix, batchs)
    print("preprocessData After removing requested batches, gene data has ", length(colnames(mymatrix)), " samples")
  }
  #############################################
  # replace NAs with zeros
  #############################################
  if (isTRUE(theReplaceNAs))
  {
    print("****replace NAs with zeros****")
    mymatrix[is.na(mymatrix)] <- 0.0
  }
  #############################################
  # log transform (normalize)
  #############################################
  if (isTRUE(theTransformFlag))
  {
    print("****Log transform data****")
    # convert to vector
    myVector <- as.vector(mymatrix)
    # so we can remove all zero values
    myVector <- myVector[myVector>0]
    # so when we remove NAs and look for .1 quantile,
    # we get a non-zero answer
    qt <- quantile(myVector, .1, na.rm=TRUE)
    print("****qt****")
    print(qt)
    cat(qt, file=theLogTransformFile)
    # that gives us non-zero (and non-infinite) values
    # within the transformed matrix
    mymatrix <- log2(mymatrix+qt)
  }
  #############################################
  # filter to reduce size
  #############################################
  mymatrix <- mbatchTrimData(mymatrix, theSize)
  #############################################
  # write to files
  #############################################
  print("****theOutputMatrixFile****")
  print(theOutputMatrixFile)
  print(paste("dim ", dim(mymatrix), collapse=",", sep=""))
  writeAsGenericMatrix(theOutputMatrixFile, mymatrix)
  print(theOutputBatchFile)
  print(paste("dim ", dim(batchs), collapse=",", sep=""))
  writeAsGenericDataframe(theOutputBatchFile, batchs)
}

bidirectionalCentering <- function(theMatrix, theType=(function(x){median(x, na.rm=TRUE)}), theSampleFeatureFlag=TRUE)
{
  if (isTRUE(theSampleFeatureFlag))
  {
    sampleVals = apply(theMatrix, 2, theType)
    theMatrix = scale(theMatrix, center=sampleVals, scale=FALSE)
    theMatrix = t(theMatrix)
    featureVals = apply(theMatrix, 2, theType)
    theMatrix = scale(theMatrix, center=featureVals, scale=FALSE)
    theMatrix = t(theMatrix)
  }
  else
  {
    theMatrix = t(theMatrix)
    featureVals = apply(theMatrix, 2, theType)
    theMatrix = scale(theMatrix, center=featureVals, scale=FALSE)
    theMatrix = t(theMatrix)
    sampleVals = apply(theMatrix, 2, theType)
    theMatrix = scale(theMatrix, center=sampleVals, scale=FALSE)
  }
  theMatrix
}

####################################################################
###
####################################################################
