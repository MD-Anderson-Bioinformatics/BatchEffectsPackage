# MBatchUtils Copyright (c) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021 University of Texas MD Anderson Cancer Center
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
    baseDir <- "/BatchEffectsPackage_data/testing_static/MBatchUtils"
  }
  baseDir
}

getTestOutputDir <- function()
{
  baseDir <- Sys.getenv("MBATCHUTILS_TEST_OUTPUT")
  if(!file.exists(baseDir))
  {
    baseDir <- "/BatchEffectsPackage_data/testing_dynamic/MBatchUtils"
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
                           theBatchTypeAndValuePairsToRemove=NULL)
{
  mymatrix <- readAsGenericMatrix(theInputMatrixFile)
  batchs <- readAsGenericDataframe(theInputBatchFile)
  #############################################
  # add missing samples from matrix to batchs
  #############################################
  matrixSamples <- colnames(mymatrix)
  batchsSamples <- batchs$Sample
  missingSamples <- setdiff(matrixSamples, batchsSamples)
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
  #writeAsGenericMatrix(theOutputMatrixFile, mymatrix)
  print("****theOutputMatrixFile****")
  print(theOutputMatrixFile)
  writeAsGenericMatrix(theOutputMatrixFile,mymatrix)
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
