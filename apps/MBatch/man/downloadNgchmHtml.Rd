\name{downloadNgchmHtml}
\alias{downloadNgchmHtml}
\title{downloadNgchmHtml}
\description{
Download the given NGHCM .html file for the given dataset.
}
\usage{
downloadNgchmHtml(theQueryPyObj, theDownloadFile, theDatasetId, theNgchmZipPath)
}
\arguments{
  \item{theQueryPyObj}{Python object to query.}
  \item{theDownloadFile}{Path and filename to which to download.}
  \item{theDatasetId}{String for dataset id to download, from the available_data entry.}
  \item{theNgchmZipPath}{String from getDownloadableForDataset describing NGCHM.html to download, from the getDownloadableForDataset NGCHM vector.}
}
\details{
Download the given NGHCM .html file for the given dataset.
}
\value{Nothing meaningful.}
\examples{
## The function is currently defined as
function (theQueryPyObj, theDownloadFile, theDatasetId, theNgchmZipPath)
{
    theQueryPyObj$download_ngchm_html_to_file(theDownloadFile,
        theDatasetId, theNgchmZipPath)
  }
}

