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

require(MBatch)

inputDir <- getTestInputDir()
outputDir <- getTestOutputDir()
compareDir <- getTestCompareDir()

invariantFile=file.path(inputDir, "rbn-test6-iset.tsv")
variantFile=file.path(inputDir, "rbn-test6-vset.tsv")
theOutputDir=file.path(outputDir, "RBN_Replicates")
theCompareFile=file.path(compareDir, "rbn-test6-output.tsv")
theRandomSeed=314
#myRandomSeed <- 314
#myTestSeed <- 42

resolveDuplicates <- function(theNames)
{
  # keep first instance of a name
  # number subsequent ones starting with .1
  make.unique(theNames)
}

readRPPAdataAsMatrix_WithTab <- function(theFile)
{
  # read RPPA data as a dataframe
  # column rppaDF[,1] contains row names that may contain duplicates
  rppaDF <- readAsGenericDataframe(theFile)
  # resolve duplicates in row names here
  myRownames <- rppaDF[,1]
  myRownames <- resolveDuplicates(myRownames)
  # convert to numeric, since in R4+ data.matrix converts character(string) to factor, and then to integer
  numDF <- as.data.frame(lapply(rppaDF[,-1], as.numeric))
  # convert to matrix
  ## data.matrix from data.frame converts character to factor to integer in R4+
  myMatrix <- data.matrix(numDF)
  rownames(myMatrix) <- myRownames
  t(myMatrix)
}

readRPPAdataAsMatrix_NoInitialTab <- function(theFile)
{
  # read RPPA data as a dataframe
  # column rppaDF[,1] contains row names that may contain duplicates
  rppaDF <- read.table(theFile, header=TRUE, sep="\t", as.is=TRUE,
                       check.names=FALSE, stringsAsFactors=FALSE,
                       colClasses="character", na.strings="NA",
                       row.names=NULL)
  # resolve duplicates in row names here
  myRownames <- rppaDF[,1]
  myRownames <- resolveDuplicates(myRownames)
  # convert to numeric, since in R4+ data.matrix converts character(string) to factor, and then to integer
  numDF <- as.data.frame(lapply(rppaDF[,-1], as.numeric))
  # convert to matrix
  ## data.matrix from data.frame converts character to factor to integer in R4+
  myMatrix <- data.matrix(numDF)
  rownames(myMatrix) <- myRownames
  t(myMatrix)
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
  unlink(theOutputDir, recursive=TRUE)
  dir.create(theOutputDir, showWarnings=FALSE, recursive=TRUE)

  message("Reading invariant file")
  invMatrix = readRPPAdataAsMatrix_WithTab(invariantFile)
  message("Reading variant file")
  varMatrix = readRPPAdataAsMatrix_WithTab(variantFile)
  filename <- RBN_Replicates(theInvariantMatrix=invMatrix,
                             theVariantMatrix=varMatrix,
                             theInvariantGroupId="Grp1",
                             theVariantGroupId="Grp2",
                             theMatchedReplicatesFlag=TRUE,
                             theCombineOnlyFlag=FALSE,
                             thePath=theOutputDir,
                             theWriteToFile=TRUE)
  correctedMatrix <- readAsGenericMatrix(file.path(theOutputDir, "ANY_Corrections-RBN_Replicates.tsv"))
  compareMatrix <- readAsGenericMatrix(theCompareFile)
  message("correctedMatrix")
  print(dim(correctedMatrix))
  print(correctedMatrix[1:4,1:3])
  message("compareMatrix")
  print(dim(compareMatrix))
  print(compareMatrix[1:4,1:3])
  compared <- compareTwoMatrices(correctedMatrix, compareMatrix)
  print(compared)
  compared
} else {
  message("No test data. Skip test.")
  TRUE
}
