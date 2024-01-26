\name{getDownloadableForDataset}
\alias{getDownloadableForDataset}
\title{getDownloadableForDataset}
\description{
Get list of available versions and the list of NGCHMs for a dataset.
}
\usage{
getDownloadableForDataset(theQueryPyObj, theDataEntry)
}
\arguments{
  \item{theQueryPyObj}{Current Python query object.}
  \item{theDataEntry}{The Data Entry object pulled from available_dataset vector in the Python query object.}
}
\details{
Get list of available versions and the list of NGCHMs for a dataset.
}
\value{
Return value[[1]] is vector of strings of data versions available.
Return value[[2]] is vector of strings of NGHCMs available.
}
\examples{
## The function is currently defined as
function (theQueryPyObj, theDataEntry)
{
    downloadableList <- theQueryPyObj$get_downloadable(theDataEntry)
    return(downloadableList)
  }
}

