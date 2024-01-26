\name{downloadDataMatrix}
\alias{downloadDataMatrix}
\title{downloadDataMatrix}
\description{
Download the given data matrix TSV file for the given dataset and version.
}
\usage{
downloadDataMatrix(theQueryPyObj, theDownloadFile, theDatasetId,
theVersion, theOriginalFlag)
}
\arguments{
  \item{theQueryPyObj}{Python object to query.}
  \item{theDownloadFile}{Path and filename to which to download.}
  \item{theDatasetId}{String for dataset id to download, from the available_data entry.}
  \item{theVersion}{String for dataset id to download, from the getDownloadableForDataset versions vector.}
  \item{theOriginalFlag}{If TRUE, download from original dataset. Otherwise, download from pipeline.}
}
\details{
Download the given data matrix TSV file for the given dataset and version.
}
\value{Nothing Meaningful.}
\examples{
## The function is currently defined as
function (theQueryPyObj, theDownloadFile, theDatasetId, theVersion,
    theOriginalFlag)
{
    theQueryPyObj$download_data_batches_to_file(theDownloadFile,
        theDatasetId, theVersion, theOriginalFlag)
  }
}

