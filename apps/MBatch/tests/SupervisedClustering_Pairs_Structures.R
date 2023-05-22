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

inputDir <- getTestInputDir()
outputDir <- getTestOutputDir()
compareDir <- getTestCompareDir()

theGeneFile=cleanFilePath(inputDir, "matrix_data-Tumor.tsv")
theBatchFile=cleanFilePath(inputDir, "batches-Tumor.tsv")
theOutputDir=cleanFilePath(outputDir, "SupervisedClustering_Pairs_Structures")
theRandomSeed=314
theBatchType="TSS"


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
  # load data
  myData <- mbatchLoadFiles(theGeneFile, theBatchFile)
  myData@mData <- mbatchTrimData(myData@mData, 100000)
  # here, we take most defaults
  SupervisedClustering_Pairs_Structures(theData=myData,
                                        theTitle="Test Supervised Clustering",
                                        theOutputDir=theOutputDir,
                                        theListOfBatchPairs=c("BatchId", "PlateId", "TSS", "ShipDate"),
                                        theDataVersion="DATA_2022-09-09-1600",
                                        theTestVersion="TEST_2022-10-10-1300",
                                        theBatchTypeAndValuePairsToRemove=NULL,
                                        theBatchTypeAndValuePairsToKeep=NULL)
  message("No comparable output--no error means 'OK' for now.")
  TRUE
} else {
  message("No test data. Skip test.")
  TRUE
}
