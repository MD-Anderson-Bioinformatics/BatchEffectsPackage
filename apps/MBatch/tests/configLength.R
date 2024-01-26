# MBatchUtils Copyright (c) 2011-2024 University of Texas MD Anderson Cancer Center
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

if (!is.null(getTestInputDir()))
{
  #warnLevel<-getOption("warn")
  #on.exit(options(warn=warnLevel))
  # warnings are errors
  #options(warn=3)
  # if there is a warning, show the calls leading up to it
  #options(showWarnCalls=TRUE)
  # if there is an error, show the calls leading up to it
  #options(showErrorCalls=TRUE)
  ########################################################
  ########################################################
  # writes to input directory, so copy files to output
  inDir <- dirname(getTestInputDir())
  inDir <- cleanFilePath(inDir, "MBatchUtils")
  originalData <- cleanFilePath(inDir, "configLength")
  matrixFile <- cleanFilePath(originalData, "matrix_data.tsv")
  batchFile <- cleanFilePath(originalData, "batches.tsv")
  message("matrixFile =",matrixFile)
  message("batchFile  =",batchFile)
  ########################################################
  # DATA and TEST version are in config file
  myMBatchData <- mbatchLoadFiles(matrixFile, batchFile)
  nrowPre <- nrow(myMBatchData@mData)
  myMBatchData <- mbatchFilterData(theBeaData=myMBatchData,
                                   theBatchTypeAndValuePairsToRemove=list(),
                                   theBatchTypeAndValuePairsToKeep=list(),
                                   theBatchTypesToRemove=list(),
                                   theMinIqr=0,
                                   theMinSd=0,
                                   theMinMad=0)
  nrowPost <- nrow(myMBatchData@mData)
  message("nrowPre=", nrowPre)
  message("nrowPost=", nrowPost)
  nrowPost==nrowPre
} else {
  message("No test data. Skip test.")
  TRUE
}
