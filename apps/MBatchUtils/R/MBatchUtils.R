#MBatchUtils Copyright ? 2018 University of Texas MD Anderson Cancer Center
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

mbatchUtilVersion <- function()
{
  return("MBatchUtil 2019-02-06-1030")
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
      message("myCol=", myCol)
      message("colnames(theCorrected)[myCol]=", colnames(theCorrected)[myCol])
      message("colnames(theCompare)[myCol]=", colnames(theCompare)[myCol])
      return(FALSE)
    }
    for(myRow in 1:correctedRows)
    {
      if (!(rownames(theCorrected)[myRow]==rownames(theCompare)[myRow]))
      {
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
      if (!is.finite(theCorrected[myRow, myCol])&&!is.finite(theCompare[myRow, myCol]))
      {
        # ignore
      }
      else if (!(all.equal(theCorrected[myRow, myCol], theCompare[myRow, myCol])))
      {
        message("myRow=", myRow)
        message("myCol=", myCol)
        message("theCorrected[myRow, myCol]=", theCorrected[myRow, myCol])
        message("theCompare[myRow, myCol]=", theCompare[myRow, myCol])
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

####################################################################
###
####################################################################

preprocessData <- function(theInputMatrixFile, theInputBatchFile, theOutputMatrixFile, theOutputBatchFile, theSize, theTransformFlag)
{
  mymatrix <- readAsGenericMatrix(theInputMatrixFile)
  batchs <- readAsDataFrame(theInputBatchFile)
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
  #############################################
  # sort, etc the data
  #############################################
  mbatchStr <- mbatchLoadStructures(mymatrix, batchs)
  mymatrix <- mbatchStr@mData
  batchs <- mbatchStr@mBatches
  #############################################
  # filter to reduce size
  #############################################
  matrix <- mbatchTrimData(mymatrix, theSize)
  #############################################
  # log transform (normalize)
  #############################################
  if (isTRUE(theTransformFlag))
  {
    # convert to vector
    myVector <- as.vector(mymatrix)
    # so we can remove all zero values
    myVector <- myVector[myVector>0]
    # so when we remove NAs and look for .1 quantile,
    # we get a non-zero answer
    qt <- quantile(myVector, .1, na.rm=TRUE)
    # that gives us non-zero (and non-infinite) values
    # within the transformed matrix
    mymatrix <- log2(mymatrix+qt)
  }
  #############################################
  # write to files
  #############################################
  #writeAsMatrix(theOutputMatrixFile, mymatrix)
  write.table(x=mymatrix, file=theOutputMatrixFile, sep="\t", na="NA", quote=FALSE, col.names=NA)
  writeAsDataframe(theOutputBatchFile, batchs)
}


####################################################################
###
####################################################################
