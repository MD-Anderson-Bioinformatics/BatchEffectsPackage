# MBatchUtils Copyright (c) 2011-2022 University of Texas MD Anderson Cancer Center
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
  # set the paths
  baseTestDir <- getTestInputDir()
  baseOutputDir <- getTestOutputDir()
  myOrigConfigFile <- cleanFilePath(cleanFilePath(baseTestDir, "corrections"), "EBN_Plus_config.tsv")
  dataDir <- cleanFilePath(dirname(baseTestDir), "MATRIX_DATA")

  data1file <- cleanFilePath(dataDir, "brca_agi4502_matrix_data.tsv")
  batch1file <- cleanFilePath(dataDir, "brca_agi4502_batches.tsv")
  data2file <- cleanFilePath(dataDir, "brca_rnaseq2_matrix_data.tsv")
  batch2file <- cleanFilePath(dataDir, "brca_rnaseq2_batches.tsv")

  myOutputDir <- cleanFilePath(baseOutputDir, "correction_EBN_Plus")
  myOutputDirMBatch <- cleanFilePath(myOutputDir, "ZIP-RESULTS")
  myOutputDirData <- cleanFilePath(myOutputDir, "ZIP-DATA")
  jarDir <- cleanFilePath(baseTestDir, "exe")
  javaExe <- getJava()
  jarFile <- cleanFilePath(jarDir, "ShaidyMapGen.jar")
  jsFile <- cleanFilePath(jarDir, "ngchmWidget-min.js")
  myDestConfigFile <- cleanFilePath(myOutputDirMBatch, "MBatchConfig.tsv")
  myDestGeneFile <- cleanFilePath(myOutputDirData, "matrix_data.tsv")
  myDestBatchFile <- cleanFilePath(myOutputDirData, "batches.tsv")
  myDestGeneFile2 <- cleanFilePath(myOutputDirData, "matrix_data2.tsv")
  myDestBatchFile2 <- cleanFilePath(myOutputDirData, "batches2.tsv")
  ####
  # make sure the output dir exists and is empty
  unlink(myOutputDir, recursive=TRUE)
  dir.create(myOutputDir, showWarnings=FALSE, recursive=TRUE)
  dir.create(myOutputDirMBatch, showWarnings=FALSE, recursive=TRUE)
  dir.create(myOutputDirData, showWarnings=FALSE, recursive=TRUE)
  # copy files to ZIP-RESULTS and ZIP-DATA
  file.copy(myOrigConfigFile, myDestConfigFile)
  file.copy(data1file, myDestGeneFile)
  file.copy(batch1file, myDestBatchFile)
  file.copy(data2file, myDestGeneFile2)
  file.copy(batch2file, myDestBatchFile2)
  # config
  mbatchRunFromConfig(
    theConfigFile = myDestConfigFile,
    theMatrixFile = myDestGeneFile,
    theBatchesFile = myDestBatchFile,
    theZipDataDir = myOutputDirData,
    theZipResultsDir = myOutputDirMBatch,
    theNaStrings = c("null", "NA"),
    theShaidyMapGen = jarFile,
    theNgchmWidgetJs = jsFile,
    theShaidyMapGenJava = javaExe,
    theNGCHMShaidyMem = "8G",
    theRunPostFlag = TRUE,
    theMatrixFile2=myDestGeneFile2,
    theBatchesFile2=myDestBatchFile2)
  unlink(myDestConfigFile)
  buildSingleArchive(myOutputDirMBatch, myOutputDirData, myOutputDir)
} else {
  message("No test data. Skip test.")
  TRUE
}
