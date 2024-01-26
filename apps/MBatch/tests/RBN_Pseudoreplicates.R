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

invariantFile=cleanFilePath(inputDir, "rbn-pseudo-iset.tsv")
variantFile=cleanFilePath(inputDir, "rbn-pseudo-vset.tsv")
theOutputDir=cleanFilePath(outputDir, "RBN_Pseudoreplicates")
theCompareFile=cleanFilePath(compareDir, "rbn-pseudo-output.tsv")
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
  sampleIds <- rppaDF[,1]
  sampleIds <- resolveDuplicates(sampleIds)
  # save features here, since as.data.frame does "r-conversion" on the strings
  # drop leading empty string
  featureIds <- colnames(rppaDF)[-1]
  featureIds <- resolveDuplicates(featureIds)
  # convert to numeric, since in R4+ data.matrix converts character(string) to factor, and then to integer
  numDF <- as.data.frame(lapply(rppaDF[,-1], as.numeric))
  # convert to matrix
  myMatrix <- data.matrix(numDF)
  rownames(myMatrix) <- sampleIds
  colnames(myMatrix) <- featureIds
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
  message(invariantFile)
  invMatrix = readRPPAdataAsMatrix_WithTab(invariantFile)
  message("invMatrix")
  print(dim(invMatrix))
  print(invMatrix[1:4,1:3])
  message("Reading variant file")
  message(variantFile)
  varMatrix = readRPPAdataAsMatrix_WithTab(variantFile)
  message("varMatrix")
  print(dim(varMatrix))
  print(varMatrix[1:4,1:3])
  invPseudo <- c("BN", "BO", "BP", "BQ", "BR", "BS")
  varPseudo <- c("AN", "AO", "AP", "AQ", "AR", "AS")
  filename <- RBN_Pseudoreplicates(theInvariantMatrix=invMatrix,
                                   theVariantMatrix=varMatrix,
                                   theInvariantReplicates = invPseudo,
                                   theVariantReplicates = varPseudo,
                                   theInvariantGroupId="Grp1",
                                   theVariantGroupId="Grp2",
                                   theMatchedReplicatesFlag=FALSE,
                                   theCombineOnlyFlag=FALSE,
                                   thePath=theOutputDir,
                                   theDataVersion="DATA_2022-09-09-1600",
                                   theTestVersion="TEST_2022-10-10-1300",
                                   theWriteToFile=TRUE)
  correctedMatrix <- readAsGenericMatrix(cleanFilePath(cleanFilePath(cleanFilePath(theOutputDir, "DATA_2022-09-09-1600"), "TEST_2022-10-10-1300"), "adjusted_matrix.tsv"))
  compareMatrix <- readAsGenericMatrix(theCompareFile)
  message("correctedMatrix")
  print(dim(correctedMatrix))
  print(correctedMatrix[1:4,1:3])
  message("compareMatrix")
  print(theCompareFile)
  print(dim(compareMatrix))
  print(compareMatrix[1:4,1:3])
  message("correctedMatrix columns")
  print(colnames(correctedMatrix))
  message("compareMatrix columns")
  print(colnames(compareMatrix))
  compared <- compareTwoMatrices(correctedMatrix, compareMatrix)
  print(compared)
  compared
} else {
  message("No test data. Skip test.")
  TRUE
}
