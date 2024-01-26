\name{updateOrBuildSeuratObjectFromStdDataFiles}
\alias{updateOrBuildSeuratObjectFromStdDataFiles}
\title{updateOrBuildSeuratObjectFromStdDataFiles}
\description{Update or build a Seurat object from Standardized Data files.}
\usage{
updateOrBuildSeuratObjectFromStdDataFiles(theSeuratObj,
theMatrixFilePath, theBatchFilePath, theAssayToUse, theDataSlotToUse)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theSeuratObj}{Seurat object to convert to a Standardized Data object. If NULL,}
  \item{theMatrixFilePath}{String: full path and filename to matrix TSV file to write.}
  \item{theBatchFilePath}{String: full path and filename to batches TSV file to write.}
  \item{theAssayToUse}{String: assay name to use, as per code seuObj@assays[[theAssayToUse]]. Not sure of possible total list other than: RNA, RNA3, RNA5}
  \item{theDataSlotToUse}{String: within assay, name of type of data to use. Also slot name in assay object. Not all types will be present. Known types are: counts, data, scale.data}
}
\details{
Update or build a Seurat object from Standardized Data files.
If selecting scale.data for theDataSlotToUse, you must pass in a good Seurat object.
}
\value{Seurat object. Or NULL.}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (theSeuratObj, theMatrixFilePath, theBatchFilePath,
    theAssayToUse, theDataSlotToUse)
{
    logDebug("updateOrBuildSeuratObjectFromStdDataFiles")
    collateOrigValue <- Sys.getlocale("LC_COLLATE")
    on.exit(Sys.setlocale("LC_COLLATE", collateOrigValue), add = TRUE)
    Sys.setlocale("LC_COLLATE", "C")
    logDebug("Changing LC_COLLATE to C for duration of run")
    checkPackageSettings()
    if (isNamespaceLoaded("Seurat")) {
        logDebug("updateOrBuildSeuratObjectFromStdDataFiles converting")
        logDebug("updateOrBuildSeuratObjectFromStdDataFiles theSeuratObj ",
            theSeuratObj)
        logDebug("updateOrBuildSeuratObjectFromStdDataFiles theMatrixFilePath ",
            theMatrixFilePath)
        logDebug("updateOrBuildSeuratObjectFromStdDataFiles theBatchFilePath ",
            theBatchFilePath)
        logDebug("updateOrBuildSeuratObjectFromStdDataFiles theAssayToUse ",
            theAssayToUse)
        logDebug("updateOrBuildSeuratObjectFromStdDataFiles theDataSlotToUse ",
            theDataSlotToUse)
        dataMatrix <- readAsGenericMatrix(theMatrixFile)
        if (is.null(theSeuratObj)) {
            if ("counts" == theDataSlotToUse) {
                theSeuratObj <- CreateSeuratObject(counts = dataMatrix)
            }
            else if ("data" == theDataSlotToUse) {
                theSeuratObj <- CreateSeuratObject(data = dataMatrix)
            }
        }
        else {
            theSeuratObj <- SetAssayData(object = theSeuratObj,
                assay = theAssayToUse, slot = theDataSlotToUse,
                new.data = dataMatrix)
        }
        if (!is.null(theBatchesFile)) {
            batchesDF <- readAsGenericDataframe(theBatchesFile)
            rownames(batchesDF) <- batchesDF$Sample
            batchesDF$Sample <- NULL
            theSeuratObj <- AddMetaData(object = theSeuratObj,
                metadata = batchesDF)
        }
        return(theSeuratObj)
    }
    else {
        logDebug("updateOrBuildSeuratObjectFromStdDataFiles Seurat package not available")
        return(NULL)
    }
  }
}
