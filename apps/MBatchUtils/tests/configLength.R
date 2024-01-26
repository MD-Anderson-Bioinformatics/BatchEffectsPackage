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

require(MBatchUtils)

message("This test generates errors by using the setGlobalMBatchErrorTest to generate artificial errors")

if (!is.null(getTestOutputDir()))
{
  setGlobalMBatchErrorTest(TRUE)
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
  theOutputResults=cleanFilePath(cleanFilePath(getTestOutputDir(), "configLength"), "ZIP-RESULTS")
  theDataDir=cleanFilePath(cleanFilePath(getTestOutputDir(), "configLength"), "ZIP-DATA")
  print(theOutputResults)
  unlink(theOutputResults, recursive=TRUE)
  dir.create(theOutputResults, recursive=TRUE, showWarnings=FALSE)
  print(theDataDir)
  unlink(theDataDir, recursive=TRUE)
  originalData <- cleanFilePath(getTestInputDir(), "configLength")
  configFile <- cleanFilePath(originalData, "MBatchConfig.tsv")
  matrixFile <- cleanFilePath(originalData, "matrix_data.tsv")
  batchFile <- cleanFilePath(originalData, "batches.tsv")
  message("configFile =",configFile)
  message("matrixFile =",matrixFile)
  message("batchFile  =",batchFile)
  ########################################################
  baseTestDir=getTestInputDir()
  jarDir=cleanFilePath(baseTestDir, "exe")
  javaExe=getJava()
  jarFile=cleanFilePath(jarDir, "ShaidyMapGen.jar")
  jsFile=cleanFilePath(jarDir, "ngchmWidget-min.js")
  ########################################################
  # DATA and TEST version are in config file
  mbatchRunFromConfig(theConfigFile=configFile,
                      theMatrixFile=matrixFile,
                      theBatchesFile=batchFile,
                      theZipDataDir=theDataDir,
                      theZipResultsDir=theOutputResults,
                      theNaStrings="NA",
                      theShaidyMapGen = jarFile,
                      theNgchmWidgetJs = jsFile,
                      theShaidyMapGenJava = javaExe,
                      theRunPostFlag=TRUE)
  # theRunPostFlag builds DSCOverview.tsv files
  # clear the RData files for DSC
  clearDSCOverviewFiles(theOutputResults)
  setGlobalMBatchErrorTest(FALSE)
  file.exists(cleanFilePath(theOutputResults, "MBATCH_SUCCESS.txt"))
} else {
  message("No test data. Skip test.")
  TRUE
}
