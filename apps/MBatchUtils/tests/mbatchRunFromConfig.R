# MBatchUtils Copyright (c) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021 University of Texas MD Anderson Cancer Center
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

if (!is.null(getTestOutputDir()))
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
  theOutputDirMBatch=file.path(getTestOutputDir(), "configout", "ZIP-RESULTS")
  theOutputDirData=file.path(getTestOutputDir(), "configout", "ZIP-DATA")
  originalData=file.path(theOutputDirData, "original")
  print(theOutputDirMBatch)
  unlink(theOutputDirMBatch, recursive=TRUE)
  dir.create(theOutputDirMBatch, recursive=TRUE, showWarnings=FALSE)
  print(theOutputDirData)
  unlink(theOutputDirData, recursive=TRUE)
  dir.create(originalData, recursive=TRUE, showWarnings=FALSE)
  file.copy(file.path(getTestInputDir(), "config", "MBatchConfig.tsv"), file.path(theOutputDirMBatch, "MBatchConfig.tsv"))
  print(file.exists(file.path(theOutputDirMBatch, "MBatchConfig.tsv")))
  file.copy(file.path(getTestInputDir(), "config", "matrix_data.tsv"), file.path(originalData, "matrix_data.tsv"))
  file.copy(file.path(getTestInputDir(), "config", "batches.tsv"), file.path(originalData, "batches.tsv"))
  configFile=file.path(theOutputDirMBatch, "MBatchConfig.tsv")
  ########################################################
  baseTestDir=getTestInputDir()
  jarDir=file.path(baseTestDir, "exe")
  javaExe=getJava()
  jarFile=file.path(jarDir, "ShaidyMapGen.jar")
  jsFile=file.path(jarDir, "ngchmWidget-min.js")
  ########################################################
  mbatchRunFromConfig(theConfigFile=configFile,
                      theDataDir=originalData,
                      theOutputDir=theOutputDirMBatch,
                      theNaStrings="NA",
                      theShaidyMapGen = jarFile,
                      theNgchmWidgetJs = jsFile,
                      theShaidyMapGenJava = javaExe)
  file.exists(file.path(theOutputDirMBatch, "MBATCH_SUCCESS.txt"))
} else {
  message("No test data. Skip test.")
  TRUE
}
