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

testUrl <- getTestDapiURL()
if (""!=testUrl)
{
  outputDir <- getTestOutputDir()
  if (!is.null(outputDir))
  {
    theOutputDir=cleanFilePath(outputDir, "dapiquery")
    ##############################################################################
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
    ##############################################################################
    myUrl <- testUrl
    pyObj <- getDapiQuery(myUrl)
    pyObj$selected_projects <- append(pyObj$selected_projects, "TCGA-LUSC")
    pyObj$selected_jobtype <- append(pyObj$selected_jobtype, "Original")
    pyObj$selected_data <- append(pyObj$selected_data, "STAR - Counts")
    updateDapiQuery(pyObj)
    # should be 3
    length(pyObj$available_datasets)
    # list is 1 indexed for Python object in R, also remove object from list
    queryEntry <- pyObj$available_datasets[[1]]
    retList <- getDownloadableForDataset(pyObj, queryEntry)
    # should be 1
    versionVector <- retList[[1]]
    length(versionVector)
    # should be 12
    ngchmVector <- retList[[2]]
    length(ngchmVector)
    # download matrix files
    datasetId <- "8ff96a845261f8455b5b2698b5db776f~2022-12-12~2022_12_28_1300"
    versionStr <- "DATA_2022-12-12"
    downloadFile <- file.path(theOutputDir, "from_r_matrix_original.tsv")
    downloadDataMatrix(pyObj, downloadFile, datasetId, versionStr, TRUE)
    downloadFile <- file.path(theOutputDir, "from_r_matrix_pipeline.tsv")
    downloadDataMatrix(pyObj, downloadFile, datasetId, versionStr, FALSE)
    downloadFile <- file.path(theOutputDir, "from_r_batches_original.tsv")
    downloadDataBatches(pyObj, downloadFile, datasetId, versionStr, TRUE)
    downloadFile <- file.path(theOutputDir, "from_r_batches_pipeline.tsv")
    downloadDataBatches(pyObj, downloadFile, datasetId, versionStr, FALSE)
    # download batch files
    datasetId <- "8ff96a845261f8455b5b2698b5db776f~2022-12-12~2022_12_28_1300"
    downloadFile <- file.path(theOutputDir, "batch_id_ngchm.ngchm")
    zipFilePath <- "/analysis/NGCHM/DATA_2022-12-12/TEST_2022_12_28_1300/All_ngchm.ngchm.html"
    downloadNgchmNgchm(pyObj, downloadFile, datasetId, zipFilePath)
    downloadFile <- file.path(theOutputDir, "batch_id_ngchm.html")
    downloadNgchmHtml(pyObj, downloadFile, datasetId, zipFilePath)
    return(TRUE)
  }
}
