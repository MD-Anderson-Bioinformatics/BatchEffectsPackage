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
  jarDir=file.path(getTestInputDir(), "exe")
  baseTestDir=file.path(getTestInputDir(), "heatmap")
  baseOutputDir=getTestOutputDir()
  javaExe=getJava()
  ####
  jarFile=file.path(jarDir, "ShaidyMapGen.jar")
  jsFile=file.path(jarDir, "ngchmWidget-min.js")
  sourceMatrix=file.path(baseTestDir, "brca_agi4502_matrix_data.tsv")
  sourceBatches=file.path(baseTestDir, "brca_agi4502_batches.tsv")
  dir.create(file.path(baseOutputDir, "NGCHM"), showWarnings=FALSE)
  outputFile=file.path(baseOutputDir, "NGCHM", "buildBatchHeatMap_Structures.ngchm")
  ########################################################
  #############################l###########################
  print(sourceMatrix)
  matData <- readAsGenericMatrix(sourceMatrix)
  print(sourceBatches)
  # need this to keep log in readAsDataFrame from failing
  setLogging(NULL)
  dfData <- readAsDataFrame(sourceBatches)
  print("call buildBatchHeatMap_Structures")
  buildBatchHeatMap_Structures(theMatrixData=matData,
                               theBatchData=dfData,
                               theTitle="Test for BRCA AGI4502 Data",
                               theOutputFile=outputFile,
                               theSortByType="BatchId",
                               theRowType="scholar", theColType="bio.tcga.barcode.sample",
                               theRowCluster=NULL, theColCluster=NULL,
                               theShaidyMapGen=jarFile,
                               theNgchmWidgetJs=jsFile,
                               theShaidyMapGenJava=javaExe,
                               theShaidyMapGenArgs="-Xmx16G")
  TRUE
}
