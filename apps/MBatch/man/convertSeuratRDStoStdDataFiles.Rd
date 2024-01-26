\name{convertSeuratRDStoStdDataFiles}
\alias{convertSeuratRDStoStdDataFiles}
\title{convertSeuratRDStoStdDataFiles}
\description{Convert a Seurat RDS file into Standardized Data format files, if Seurat is available.}
\usage{
convertSeuratRDStoStdDataFiles(theSeuratRdsFilePath, theMatrixFilePath,
theBatchFilePath, theAssayToUse, theDataSlotToUse)
}
\arguments{
  \item{theSeuratRdsFilePath}{String: full path and filename to Seurat RDS file to process.}
  \item{theMatrixFilePath}{String: full path and filename to matrix TSV file to write.}
  \item{theBatchFilePath}{String: full path and filename to batches TSV file to write.}
  \item{theAssayToUse}{String: assay name to use, as per code seuObj@assays[[theAssayToUse]]. Not sure of possible total list other than: RNA, RNA3, RNA5}
  \item{theDataSlotToUse}{String: within assay, name of type of data to use. Also slot name in assay object. Not all types will be present. Known types are: counts, data, scale.data}
}
\details{Convert a Seurat RDS file into Standardized Data format files, if Seurat is available.}
\value{
TRUE is process was performed and completed.
FALSE if Seurat package was not present, and process was not performed.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (theSeuratRdsFilePath, theMatrixFilePath, theBatchFilePath,
    theAssayToUse, theDataSlotToUse)
{
    logDebug("convertSeuratRDStoStdDataFiles")
    collateOrigValue <- Sys.getlocale("LC_COLLATE")
    on.exit(Sys.setlocale("LC_COLLATE", collateOrigValue), add = TRUE)
    Sys.setlocale("LC_COLLATE", "C")
    logDebug("Changing LC_COLLATE to C for duration of run")
    checkPackageSettings()
    if (require(Seurat, warn.conflicts = FALSE)) {
        logDebug("convertSeuratRDStoStdDataFiles converting")
        logDebug("convertSeuratRDStoStdDataFiles theSeuratRdsFilePath ",
            theSeuratRdsFilePath)
        logDebug("convertSeuratRDStoStdDataFiles theMatrixFilePath ",
            theMatrixFilePath)
        logDebug("convertSeuratRDStoStdDataFiles theBatchFilePath ",
            theBatchFilePath)
        logDebug("convertSeuratRDStoStdDataFiles theAssayToUse ",
            theAssayToUse)
        logDebug("convertSeuratRDStoStdDataFiles theDataSlotToUse ",
            theDataSlotToUse)
        seuObj <- readRDS(theSeuratRdsFilePath)
        assayObj <- seuObj@assays[[theAssayToUse]]
        assayData <- slot(assayObj, theDataSlotToUse)
        assayMatrix <- as.matrix(assayData)
        metaData <- seuObj@meta.data
        metaDF <- as.data.frame(metaData)
        metaDF$Sample <- rownames(metaDF)
        metaDF <- metaDF[, c(ncol(metaDF), 1:(ncol(metaDF) -
            1))]
        rownames(metaDF) <- NULL
        writeAsGenericMatrix(theMatrixFilePath, assayMatrix)
        writeAsGenericDataframe(theBatchFilePath, metaDF)
        return(TRUE)
    }
    else {
        logDebug("convertSeuratRDStoStdDataFiles Seurat package not available")
        return(FALSE)
    }
  }
}
