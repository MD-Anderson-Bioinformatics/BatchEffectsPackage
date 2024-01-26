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

getDapiQuery <- function(theUrl)
{
  dapiQueryImport <- import("mbatch.dapi.query")
  queryPyObj <- dapiQueryImport$DapiQuery(theUrl)
  queryPyObj$update_from_selected()
  return(queryPyObj)
}

updateDapiQuery <- function(theQueryPyObj)
{
  theQueryPyObj$update_from_selected()
  return(theQueryPyObj)
}

getDownloadableForDataset <- function(theQueryPyObj, theDataEntry)
{
  downloadableList <- theQueryPyObj$get_downloadable(theDataEntry)
  # list[[1]] is vector of data versions available
  # list[[2]] is vector of NGHCMs available
  return(downloadableList)
}

downloadDataMatrix <- function(theQueryPyObj, theDownloadFile, theDatasetId, theVersion, theOriginalFlag)
{
  theQueryPyObj$download_data_batches_to_file(theDownloadFile, theDatasetId, theVersion, theOriginalFlag)
}

downloadDataBatches <- function(theQueryPyObj, theDownloadFile, theDatasetId, theVersion, theOriginalFlag)
{
  theQueryPyObj$download_data_matrix_to_file(theDownloadFile, theDatasetId, theVersion, theOriginalFlag)
}

downloadNgchmNgchm <- function(theQueryPyObj, theDownloadFile, theDatasetId, theNgchmZipPath)
{
  theQueryPyObj$download_ngchm_ngchm_to_file(theDownloadFile, theDatasetId, theNgchmZipPath)
}

downloadNgchmHtml <- function(theQueryPyObj, theDownloadFile, theDatasetId, theNgchmZipPath)
{
  theQueryPyObj$download_ngchm_html_to_file(theDownloadFile, theDatasetId, theNgchmZipPath)
}
