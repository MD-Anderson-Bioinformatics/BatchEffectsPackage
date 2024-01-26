\name{updateOrBuildSeuratObjectFromStdDataObjs}
\alias{updateOrBuildSeuratObjectFromStdDataObjs}
\title{updateOrBuildSeuratObjectFromStdDataObjs}
\description{Update or build a Seurat object from Standardized Data objects.}
\usage{
updateOrBuildSeuratObjectFromStdDataObjs(theSeuratObj, theDataMatrix,
theBatchDF, theAssayToUse, theDataSlotToUse)
}
\arguments{
  \item{theSeuratObj}{Seurat object to convert to a Standardized Data object. If NULL,}
  \item{theDataMatrix}{Matrix object: Standardized Data matrix of data information.}
  \item{theBatchDF}{Dataframe object: Standardized Data dataframe of batch information.}
  \item{theAssayToUse}{String: assay name to use, as per code seuObj@assays[[theAssayToUse]]. Not sure of possible total list other than: RNA, RNA3, RNA5}
  \item{theDataSlotToUse}{String: within assay, name of type of data to use. Also slot name in assay object. If selecting scale.data, you must pass in a good Seurat object. Known types are: counts, data, scale.data}
}
\details{
Update or build a Seurat object from Standardized Data objects.
If selecting scale.data for theDataSlotToUse, you must pass in a good Seurat object.
}
\value{Seurat object. Or NULL.}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (theSeuratObj, theDataMatrix, theBatchDF, theAssayToUse,
    theDataSlotToUse)
{
    logDebug("updateOrBuildSeuratObjectFromStdDataObjs")
    collateOrigValue <- Sys.getlocale("LC_COLLATE")
    on.exit(Sys.setlocale("LC_COLLATE", collateOrigValue), add = TRUE)
    Sys.setlocale("LC_COLLATE", "C")
    logDebug("Changing LC_COLLATE to C for duration of run")
    checkPackageSettings()
    if (isNamespaceLoaded("Seurat")) {
        logDebug("updateOrBuildSeuratObjectFromStdDataObjs converting")
        logDebug("updateOrBuildSeuratObjectFromStdDataObjs theSeuratObj ",
            theSeuratObj)
        logDebug("updateOrBuildSeuratObjectFromStdDataObjs theDataMatrix ",
            theDataMatrix)
        logDebug("updateOrBuildSeuratObjectFromStdDataObjs theBatchDF ",
            theBatchDF)
        logDebug("updateOrBuildSeuratObjectFromStdDataObjs theAssayToUse ",
            theAssayToUse)
        logDebug("updateOrBuildSeuratObjectFromStdDataObjs theDataSlotToUse ",
            theDataSlotToUse)
        if (is.null(theSeuratObj)) {
            if ("counts" == theDataSlotToUse) {
                theSeuratObj <- CreateSeuratObject(counts = theDataMatrix)
            }
            else if ("data" == theDataSlotToUse) {
                theSeuratObj <- CreateSeuratObject(data = theDataMatrix)
            }
        }
        else {
            theSeuratObj <- SetAssayData(object = theSeuratObj,
                assay = theAssayToUse, slot = theDataSlotToUse,
                new.data = theDataMatrix)
        }
        if (!is.null(theBatchDF)) {
            rownames(theBatchDF) <- theBatchDF$Sample
            theBatchDF$Sample <- NULL
            theSeuratObj <- AddMetaData(object = theSeuratObj,
                metadata = theBatchDF)
        }
        return(theSeuratObj)
    }
    else {
        logDebug("updateOrBuildSeuratObjectFromStdDataObjs Seurat package not available")
        return(NULL)
    }
  }
}
