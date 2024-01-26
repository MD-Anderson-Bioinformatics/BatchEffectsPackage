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


convertSeuratRDStoStdDataFiles <- function(theSeuratRdsFilePath, theMatrixFilePath, theBatchFilePath,
                                         theAssayToUse, theDataSlotToUse)
{
  logDebug("convertSeuratRDStoStdDataFiles")
  collateOrigValue<-Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
  Sys.setlocale("LC_COLLATE","C")
  logDebug("Changing LC_COLLATE to C for duration of run")
  checkPackageSettings()
  if(require(Seurat, warn.conflicts=FALSE))
  {
    logDebug("convertSeuratRDStoStdDataFiles converting")
    logDebug("convertSeuratRDStoStdDataFiles theSeuratRdsFilePath ", theSeuratRdsFilePath)
    logDebug("convertSeuratRDStoStdDataFiles theMatrixFilePath ", theMatrixFilePath)
    logDebug("convertSeuratRDStoStdDataFiles theBatchFilePath ", theBatchFilePath)
    logDebug("convertSeuratRDStoStdDataFiles theAssayToUse ", theAssayToUse)
    logDebug("convertSeuratRDStoStdDataFiles theDataSlotToUse ", theDataSlotToUse)
    # READ THE SEURAT FILE
    seuObj <- readRDS(theSeuratRdsFilePath)

    # BUILD THE DATA MATRIX
    assayObj <- seuObj@assays[[theAssayToUse]]
    # select slot with matrix info (usually: counts, data, scale.data)
    assayData <- methods::slot(assayObj, theDataSlotToUse)
    assayMatrix <- as.matrix(assayData)

    # BUILD THE BATCHES DATAFRAME
    metaData <- seuObj@meta.data
    # convert matrix to dataframe
    metaDF <- as.data.frame(metaData)
    # add label for first column
    metaDF$Sample <- rownames(metaDF)
    # move Row column to the front
    metaDF <- metaDF[, c(ncol(metaDF), 1:(ncol(metaDF)-1))]
    # remove rownames
    rownames(metaDF) <- NULL

    # SAVE IN STANDARDIZED DATA FORMAT
    writeAsGenericMatrix(theMatrixFilePath, assayMatrix)
    writeAsGenericDataframe(theBatchFilePath, metaDF)
    return(TRUE)
  }
  else
  {
    logDebug("convertSeuratRDStoStdDataFiles Seurat package not available")
    return(FALSE)
  }
}

convertSeuratObjtoStdDataObj <- function(theSeuratObj, theAssayToUse, theDataSlotToUse)
{
  logDebug("convertSeuratObjtoStdDataObj")
  collateOrigValue<-Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
  Sys.setlocale("LC_COLLATE","C")
  logDebug("Changing LC_COLLATE to C for duration of run")
  checkPackageSettings()
  if(isNamespaceLoaded("Seurat"))
  {
    logDebug("convertSeuratObjtoStdDataObj converting")
    logDebug("convertSeuratObjtoStdDataObj theAssayToUse ", theAssayToUse)
    logDebug("convertSeuratObjtoStdDataObj theDataSlotToUse ", theDataSlotToUse)
    # BUILD THE DATA MATRIX
    assayObj <- theSeuratObj@assays[[theAssayToUse]]
    # select slot with matrix info (usually: counts, data, scale.data)
    assayData <- methods::slot(assayObj, theDataSlotToUse)
    assayMatrix <- as.matrix(assayData)

    # BUILD THE BATCHES DATAFRAME
    metaData <- theSeuratObj@meta.data
    # convert matrix to dataframe
    metaDF <- as.data.frame(metaData)
    # add label for first column
    metaDF$Sample <- rownames(metaDF)
    # move Row column to the front
    metaDF <- metaDF[, c(ncol(metaDF), 1:(ncol(metaDF)-1))]
    # remove rownames
    rownames(metaDF) <- NULL

    # CONVERT TO STANDARDIZED DATA FORMAT
    # R objects get clunky with optional slots, so pass empty dataframe to covariates
    myData <- new("BEA_DATA", sortMatrix(assayMatrix), sortDataframe(metaDF), data.frame())
    return(myData)
  }
  else
  {
    logDebug("convertSeuratObjtoStdDataObj Seurat package not available")
    return(NULL)
  }
}

updateOrBuildSeuratObjectFromStdDataFiles <- function(theSeuratObj, theMatrixFilePath, theBatchFilePath,
                                                  theAssayToUse, theDataSlotToUse)
{
  logDebug("updateOrBuildSeuratObjectFromStdDataFiles")
  collateOrigValue<-Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
  Sys.setlocale("LC_COLLATE","C")
  logDebug("Changing LC_COLLATE to C for duration of run")
  checkPackageSettings()
  if(isNamespaceLoaded("Seurat"))
  {
    logDebug("updateOrBuildSeuratObjectFromStdDataFiles converting")
    logDebug("updateOrBuildSeuratObjectFromStdDataFiles theSeuratObj ", theSeuratObj)
    logDebug("updateOrBuildSeuratObjectFromStdDataFiles theMatrixFilePath ", theMatrixFilePath)
    logDebug("updateOrBuildSeuratObjectFromStdDataFiles theBatchFilePath ", theBatchFilePath)
    logDebug("updateOrBuildSeuratObjectFromStdDataFiles theAssayToUse ", theAssayToUse)
    logDebug("updateOrBuildSeuratObjectFromStdDataFiles theDataSlotToUse ", theDataSlotToUse)
    # use "counts" for unnormalized data, use "data" for normalized data
    # if you want "scale.data", you need an existing Seurat object and assay to which to add it
    dataMatrix <- readAsGenericMatrix(theMatrixFilePath)
    if (is.null(theSeuratObj))
    {
      if ("counts"==theDataSlotToUse)
      {
        theSeuratObj <- Seurat::CreateSeuratObject(counts=dataMatrix)
      }
      else if ("data"==theDataSlotToUse)
      {
        theSeuratObj <- Seurat::CreateSeuratObject(data=dataMatrix)
      }
    }
    else
    {
      theSeuratObj <- Seurat::SetAssayData(object=theSeuratObj, assay=theAssayToUse, slot=theDataSlotToUse, new.data=dataMatrix)
    }
    if (!is.null(theBatchFilePath))
    {
      batchesDF <- readAsGenericDataframe(theBatchFilePath)
      # set rownames for dataframe
      rownames(batchesDF) <- batchesDF$Sample
      # remove the Sample column
      batchesDF$Sample <- NULL
      theSeuratObj <- Seurat::AddMetaData(object=theSeuratObj, metadata=batchesDF)
    }
    return(theSeuratObj)
  }
  else
  {
    logDebug("updateOrBuildSeuratObjectFromStdDataFiles Seurat package not available")
    return(NULL)
  }
}

updateOrBuildSeuratObjectFromStdDataObjs <- function(theSeuratObj, theDataMatrix, theBatchDF,
                                                     theAssayToUse, theDataSlotToUse)
{
  logDebug("updateOrBuildSeuratObjectFromStdDataObjs")
  collateOrigValue<-Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
  Sys.setlocale("LC_COLLATE","C")
  logDebug("Changing LC_COLLATE to C for duration of run")
  checkPackageSettings()
  if(isNamespaceLoaded("Seurat"))
  {
    logDebug("updateOrBuildSeuratObjectFromStdDataObjs converting")
    logDebug("updateOrBuildSeuratObjectFromStdDataObjs theSeuratObj ", theSeuratObj)
    logDebug("updateOrBuildSeuratObjectFromStdDataObjs theDataMatrix ", theDataMatrix)
    logDebug("updateOrBuildSeuratObjectFromStdDataObjs theBatchDF ", theBatchDF)
    logDebug("updateOrBuildSeuratObjectFromStdDataObjs theAssayToUse ", theAssayToUse)
    logDebug("updateOrBuildSeuratObjectFromStdDataObjs theDataSlotToUse ", theDataSlotToUse)
    # use "counts" for unnormalized data, use "data" for normalized data
    # if you want "scale.data", you need an existing Seurat object and assay to which to add it
    if (is.null(theSeuratObj))
    {
      if ("counts"==theDataSlotToUse)
      {
        theSeuratObj <- Seurat::CreateSeuratObject(counts=theDataMatrix)
      }
      else if ("data"==theDataSlotToUse)
      {
        theSeuratObj <- Seurat::CreateSeuratObject(data=theDataMatrix)
      }
    }
    else
    {
      theSeuratObj <- Seurat::SetAssayData(object=theSeuratObj, assay=theAssayToUse, slot=theDataSlotToUse, new.data=theDataMatrix)
    }
    if (!is.null(theBatchDF))
    {
      # set rownames for dataframe
      rownames(theBatchDF) <- theBatchDF$Sample
      # remove the Sample column
      theBatchDF$Sample <- NULL
      theSeuratObj <- Seurat::AddMetaData(object=theSeuratObj, metadata=theBatchDF)
    }
    return(theSeuratObj)
  }
  else
  {
    logDebug("updateOrBuildSeuratObjectFromStdDataObjs Seurat package not available")
    return(NULL)
  }
}
