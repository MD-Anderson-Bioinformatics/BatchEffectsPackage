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
  # set the paths
  baseTestDir=getTestInputDir()
  baseOutputDir=getTestOutputDir()
  theOrigConfigFile=file.path(baseTestDir, "vignettes", "MBatchConfig.tsv")
  theOrigGeneFile=file.path(baseTestDir, "vignettes", "matrix_data.tsv")
  theOrigBatchFile=file.path(baseTestDir, "vignettes", "batches.tsv")
  theOutputDir=file.path(baseOutputDir, "ZIP_Archive")
  theOutputDirMBatch=file.path(theOutputDir, "ZIP-RESULTS")
  theOutputDirData=file.path(theOutputDir, "ZIP-DATA")
  jarDir=file.path(baseTestDir, "exe")
  javaExe=getJava()
  jarFile=file.path(jarDir, "ShaidyMapGen.jar")
  jsFile=file.path(jarDir, "ngchmWidget-min.js")
  theDestConfigFile=file.path(theOutputDirMBatch, "MBatchConfig.tsv")
  theDestGeneFile=file.path(theOutputDirData, "matrix_data.tsv")
  theDestBatchFile=file.path(theOutputDirData, "batches.tsv")
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
    theDataDir = theOutputDirData,
    theOutputDir = theOutputDirMBatch,
    theNaStrings = c("null", "NA"),
    theShaidyMapGen = jarFile,
    theNgchmWidgetJs = jsFile,
    theShaidyMapGenJava = javaExe,
    theNGCHMShaidyMem = "8G",
    thePCAMem = "4800m",
    theBoxplotMem = "8G",
    theRunPostFlag = TRUE
  )
  theMbatchID <- "none"
  theResultDir <- theOutputDirMBatch
  theDataDir <- theOutputDirData
  theZipDir <- theOutputDir
  buildSingleArchive(theMbatchID, theResultDir, theDataDir, theZipDir)
} else {
  message("No test data. Skip test.")
  TRUE
}
