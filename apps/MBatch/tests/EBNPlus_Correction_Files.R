# MBatch Copyright (c) 2011-2022 University of Texas MD Anderson Cancer Center
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

print(getwd())

inputDir <- getTestInputDir()
outputDir <- getTestOutputDir()
compareDir <- getTestCompareDir()
print(inputDir)
print(outputDir)
print(compareDir)

theDataFile1=cleanFilePath(inputDir, "brca_rnaseq2_matrix_data.tsv")
theDataFile2=cleanFilePath(inputDir, "brca_agi4502_matrix_data.tsv")
theOutputDir=cleanFilePath(outputDir, "ebnplus")
dir.create(theOutputDir, recursive=TRUE, showWarnings=FALSE)
theCompareFile=cleanFilePath(compareDir, "EBNPlus_Correction_Files.tsv.zip")
theCompareFilename="EBNPlus_Correction_Files.tsv"
theBatchId1="RNASeqV2"
theBatchId2="Agilent4502"
theRandomSeed=314
#myRandomSeed <- 314
#myTestSeed <- 42

if (!is.null(inputDir))
{
  message("EBNPlus_Correction_Files")
  warnLevel<-getOption("warn")
  on.exit(options(warn=warnLevel))
  # warnings are errors
  options(warn=3)
  # if there is a warning, show the calls leading up to it
  options(showWarnCalls=TRUE)
  # if there is an error, show the calls leading up to it
  options(showErrorCalls=TRUE)
  #
  outdir <- cleanFilePath(theOutputDir, "EBNPlus_Correction_Files")
  unlink(outdir, recursive=TRUE)
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  #setLogging(new("Logging", theFile=cleanFilePath(outdir, "mbatch.log")))
  # this is an MDA function that starts with and processes standardized data files
  myCorrectedFile <- EBNPlus_Correction_Files(
    theDataFile1=theDataFile1,
    theDataFile2=theDataFile2,
    theOutputDir=outdir,
    theDataVersion="DATA_2022-09-09-1600",
    theTestVersion="TEST_2022-10-10-1300",
    theBatchId1=theBatchId1,
    theBatchId2=theBatchId2,
    theSeed=theRandomSeed,
    theEBNP_PriorPlotsFlag=TRUE)
  message("after correction-load file")
  myCorrectedFile <- myCorrectedFile[[1]]
  myRenamedFile <- cleanFilePath(dirname(myCorrectedFile), "corrected.tsv")
  file.rename(myCorrectedFile, myRenamedFile)
  correctedMatrix <- readAsGenericMatrix(myRenamedFile)
  compareMatrix <- as.matrix(read.delim(unz(theCompareFile, theCompareFilename), header=TRUE, sep="\t", as.is=TRUE, check.names=FALSE, stringsAsFactors=FALSE, row.names=1))
  message("myRenamedFile=",myRenamedFile)
  message("theCompareFile=",theCompareFile)
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
