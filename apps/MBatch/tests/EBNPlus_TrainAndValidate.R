# MBatch Copyright (c) 2011-2024 University of Texas MD Anderson Cancer Center
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
# MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>

require(MBatch)

inputDir <- getTestInputDir()
outputDir <- getTestOutputDir()
compareDir <- getTestCompareDir()

theDataFile1=cleanFilePath(inputDir, "brca_rnaseq2_matrix_data.tsv")
theDataFile2=cleanFilePath(inputDir, "brca_agi4502_matrix_data.tsv")
theOutputDir=cleanFilePath(outputDir, "EBNPlus_TrainAndValidateReplicates_Structures3")
theCompareFile=cleanFilePath(compareDir, "EBNPlus_TrainAndValidate.tsv")

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
    thePriorPlotPath=theOutputDir,
    theDataVersion="DATA_2022-09-09-1600",
    theTestVersion="TEST_2022-10-10-1300",
    thePriorPlotFile="priorplots.PNG")
  printMatrix(resultsList$ValidationResults)
  testMe <- resultsList$ValidationResults
  compareMe <- readAsGenericMatrix(theCompareFile)
  compared <- compareTwoMatrices(testMe, compareMe)
  print("compared")
  print(compared)
  compared
}

TRUE
