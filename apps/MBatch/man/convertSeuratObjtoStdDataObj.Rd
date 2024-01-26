\name{convertSeuratObjtoStdDataObj}
\alias{convertSeuratObjtoStdDataObj}
\title{convertSeuratObjtoStdDataObj}
\description{Convert a Seurat object into Standardized Data format files, if Seurat is available.}
\usage{
convertSeuratObjtoStdDataObj(theSeuratObj, theAssayToUse, theDataSlotToUse)
}
\arguments{
  \item{theSeuratObj}{Seurat object to convert to a Standardized Data object.}
  \item{theAssayToUse}{String: assay name to use, as per code seuObj@assays[[theAssayToUse]]. Not sure of possible total list other than: RNA, RNA3, RNA5}
  \item{theDataSlotToUse}{String: within assay, name of type of data to use. Also slot name in assay object. Not all types will be present. Known types are: counts, data, scale.data}
}
\details{Convert a Seurat object into Standardized Data format files, if Seurat is available.}
\value{
TRUE is process was performed and completed.
FALSE if Seurat package was not present, and process was not performed.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (theSeuratObj, theAssayToUse, theDataSlotToUse)
{
    logDebug("convertSeuratObjtoStdDataObj")
    collateOrigValue <- Sys.getlocale("LC_COLLATE")
    on.exit(Sys.setlocale("LC_COLLATE", collateOrigValue), add = TRUE)
    Sys.setlocale("LC_COLLATE", "C")
    logDebug("Changing LC_COLLATE to C for duration of run")
    checkPackageSettings()
    if (isNamespaceLoaded("Seurat")) {
        logDebug("convertSeuratObjtoStdDataObj converting")
        logDebug("convertSeuratObjtoStdDataObj theAssayToUse ",
            theAssayToUse)
        logDebug("convertSeuratObjtoStdDataObj theDataSlotToUse ",
            theDataSlotToUse)
        assayObj <- seuObj@assays[[theAssayToUse]]
        assayData <- slot(assayObj, theDataSlotToUse)
        assayMatrix <- as.matrix(assayData)
        metaData <- seuObj@meta.data
        metaDF <- as.data.frame(metaData)
        metaDF$Sample <- rownames(metaDF)
        metaDF <- metaDF[, c(ncol(metaDF), 1:(ncol(metaDF) -
            1))]
        rownames(metaDF) <- NULL
        myData <- new("BEA_DATA", sortMatrix(assayMatrix), sortDataframe(metaDF),
            data.frame())
        return(myData)
    }
    else {
        logDebug("convertSeuratObjtoStdDataObj Seurat package not available")
        return(NULL)
    }
  }
}
