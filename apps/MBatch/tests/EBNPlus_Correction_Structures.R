# MBatch Copyright (c) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020 University of Texas MD Anderson Cancer Center
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
# MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>

library(MBatch)

inputDir <- getTestInputDir()
outputDir <- getTestOutputDir()
compareDir <- getTestCompareDir()

theDataFile1=file.path(inputDir, "brca_rnaseq2_matrix_data.tsv")
theDataFile2=file.path(inputDir, "brca_agi4502_matrix_data.tsv")
theOutputDir=file.path(outputDir, "ebnplus")
theCompareFile=file.path(compareDir, "EBNPlus_Correction_Structures.tsv.zip")
theCompareFilename="EBNPlus_Correction_Structures.tsv"
theBatchId1="RNASeqV2"
theBatchId2="Agilent4502"
theRandomSeed=314
#myRandomSeed <- 314
#myTestSeed <- 42


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

if (!is.null(inputDir))
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
  outdir <- file.path(theOutputDir, "EBNPlus_Correction_Structures")
  unlink(outdir, recursive=TRUE)
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  # read the files in. This can be done however you want
  theDataMatrix1 <- readAsGenericMatrix(theDataFile1)
  theDataMatrix2 <- readAsGenericMatrix(theDataFile2)
  # this is the reduce genes to just gene symbols, handling those from standardized data
  rownames(theDataMatrix1) <- trimGenes(rownames(theDataMatrix1))
  rownames(theDataMatrix2) <- trimGenes(rownames(theDataMatrix2))
  # remove any duplicates (this is a requirement for EBNplus)
  theDataMatrix1 <- removeDuplicatesFromColumns(removeDuplicatesFromRows(theDataMatrix1))
  theDataMatrix2 <- removeDuplicatesFromColumns(removeDuplicatesFromRows(theDataMatrix2))
  correctedMatrix <- EBNPlus_Correction_Structures(theDataMatrix1, theDataMatrix2, theBatchId1, theBatchId2,
                                       theEBNP_BatchWithZero="1",
                                       theEBNP_FixDataSet=as.numeric(NA),
                                       theEBNP_CorrectForZero=TRUE,
                                       theEBNP_ParametricPriorsFlag=TRUE,
                                       theSeed=theRandomSeed,
                                       theEBNP_PriorPlotsFile=file.path(outdir, "priorplots.PNG"))
  outputFile <- file.path(outdir, "corrected.tsv")
  writeAsMatrix(outputFile, correctedMatrix)
  compareMatrix <- as.matrix(read.delim(unz(theCompareFile, theCompareFilename), header=TRUE, sep="\t", as.is=TRUE, check.names=FALSE, stringsAsFactors=FALSE, row.names=1))
  compared <- compareTwoMatrices(correctedMatrix, compareMatrix)
  print(compared)
  compared
} else {
  message("No test data. Skip test.")
  TRUE
}
