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

theGeneFile=cleanFilePath(inputDir, "brca_agi4502_matrix_data.tsv")
theBatchFile=cleanFilePath(inputDir, "brca_agi4502_batches.tsv")
theOutputDir=cleanFilePath(outputDir, "TRINOVA_Structures")
theCompareFile=cleanFilePath(compareDir, "TRINOVA_Structures.tsv")
theRandomSeed=314

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
  unlink(theOutputDir, recursive=TRUE, force=TRUE)
  dir.create(theOutputDir, showWarnings=FALSE, recursive=TRUE)
  # load data
  myData <- mbatchLoadFiles(theGeneFile, theBatchFile)
  myData@mData <- mbatchTrimData(myData@mData, 100000)
  myData@mBatches <- myData@mBatches[,c("Sample", "BatchId", "TSS", "ShipDate")]
  # here, we take most defaults
  outFile = TRINOVA_Structures(theData=myData,
                           theTitle="Test",
                           theOutputPath=theOutputDir,
                           theBatchTypeAndValuePairsToRemove=NULL,
                           theBatchTypeAndValuePairsToKeep=NULL,
                           theMaxGeneCount=10000)
  newDF <- readAsGenericDataframe(outFile)
  compareDF <- readAsGenericDataframe(theCompareFile)
  message("newDF dim")
  print(dim(newDF))
  message("compareDF dim")
  print(dim(compareDF))
  message("newDF content")
  print(newDF[1:4,1:2])
  message("compareDF content")
  print(compareDF[1:4,1:2])
  compared <- compareTwoMatrices(newDF, compareDF)
  print(compared)
  compared
} else {
  message("No test data. Skip test.")
  TRUE
}
