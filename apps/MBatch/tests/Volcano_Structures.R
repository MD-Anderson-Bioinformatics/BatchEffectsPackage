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

theGeneFile=cleanFilePath(inputDir, "matrix_data-Tumor.tsv")
theBatchFile=cleanFilePath(inputDir, "batches-Tumor.tsv")
theOutputDir=cleanFilePath(outputDir, "Volcano_Structures")
theCompareFile=cleanFilePath(compareDir, "Volcano-Data-OR_-_University_of_Michigan.json")
theDynamicFile=cleanFilePath(cleanFilePath(cleanFilePath(cleanFilePath(theOutputDir, "TSS"), "DATA_2022-09-09-1600"), "TEST_2022-10-10-1300"), "Volcano-Data-OR_-_University_of_Michigan.json")
theRandomSeed=314
#myRandomSeed <- 314
#myTestSeed <- 42
#theBatchType="TSS"

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
  # here, we take most defaults
  Volcano_Structures(theData=myData,
                     theTitle="Test",
                     theOutputDir=theOutputDir,
                     theLogFrameFlag=FALSE,
                     theBatchTypeAndValuePairsToRemove=NULL,
                     theBatchTypeAndValuePairsToKeep=NULL,
                     theDataVersion="DATA_2022-09-09-1600",
                     theTestVersion="TEST_2022-10-10-1300",
                     theMaxFeatureCount=50000)
  # read Volcano-Data-OR_-_University_of_Michigan.json
  message("read files to compare")
  print(theCompareFile)
  staticLines <- readLines(theCompareFile, warn=FALSE)
  print(theDynamicFile)
  dynamicLines <- readLines(theDynamicFile, warn=FALSE)
  same <- FALSE
  message("compare")
  if (all.equal(staticLines, dynamicLines))
  {
    message("matched")
    same <- TRUE
  }
  same
} else {
  message("No test data. Skip test.")
  TRUE
}
