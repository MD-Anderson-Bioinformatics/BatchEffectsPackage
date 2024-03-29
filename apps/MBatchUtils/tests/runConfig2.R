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

#
# Is also a test for having Title EMPTY in the MBatchConfig.tsv file
#

if (!is.null(getTestOutputDir()))
{
  # set the paths
  baseTestDir=getTestInputDir()
  baseOutputDir=getTestOutputDir()
  theOrigConfigFile=cleanFilePath(cleanFilePath(baseTestDir, "vignettes"), "MBatchConfig.tsv")
  theOrigGeneFile=cleanFilePath(cleanFilePath(baseTestDir, "vignettes"), "matrix_data.tsv")
  theOrigBatchFile=cleanFilePath(cleanFilePath(baseTestDir, "vignettes"), "batches.tsv")
  theOutputDir=cleanFilePath(baseOutputDir, "ZIP_Archive")
  theOutputDirMBatch=cleanFilePath(theOutputDir, "ZIP-RESULTS")
  theOutputDirData=cleanFilePath(theOutputDir, "ZIP-DATA")
  jarDir=cleanFilePath(baseTestDir, "exe")
  javaExe=getJava()
  jarFile=cleanFilePath(jarDir, "ShaidyMapGen.jar")
  jsFile=cleanFilePath(jarDir, "ngchmWidget-min.js")
  theDestConfigFile=cleanFilePath(theOutputDirMBatch, "MBatchConfig.tsv")
  theDestGeneFile=cleanFilePath(theOutputDirData, "matrix_data.tsv")
  theDestBatchFile=cleanFilePath(theOutputDirData, "batches.tsv")
  ####
  # make sure the output dir exists and is empty
  unlink(theOutputDir, recursive=TRUE)
  dir.create(theOutputDir, showWarnings=FALSE, recursive=TRUE)
  dir.create(theOutputDirMBatch, showWarnings=FALSE, recursive=TRUE)
  dir.create(theOutputDirData, showWarnings=FALSE, recursive=TRUE)
  # copy files to ZIP-RESULTS and ZIP-DATA
  file.copy(theOrigConfigFile, theDestConfigFile)
  file.copy(theOrigGeneFile, theDestGeneFile)
  file.copy(theOrigBatchFile, theDestBatchFile)
  # config
  mbatchRunFromConfig(
    theConfigFile = theDestConfigFile,
    theMatrixFile = theOrigGeneFile,
    theBatchesFile = theOrigBatchFile,
    theZipDataDir = theOutputDirData,
    theZipResultsDir = theOutputDirMBatch,
    theNaStrings = c("null", "NA"),
    theShaidyMapGen = jarFile,
    theNgchmWidgetJs = jsFile,
    theShaidyMapGenJava = javaExe,
    theNGCHMShaidyMem = "8G",
    theRunPostFlag = TRUE)
  buildSingleArchive(theOutputDirMBatch, theOutputDirData, theOutputDir)
} else {
  message("No test data. Skip test.")
  TRUE
}
