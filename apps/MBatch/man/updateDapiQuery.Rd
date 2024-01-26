\name{updateDapiQuery}
\alias{updateDapiQuery}
\title{updateDapiQuery}
\description{
Update available_datasets and other available members of the Python object based on the selected members.
}
\usage{
updateDapiQuery(theQueryPyObj)
}
\arguments{
  \item{theQueryPyObj}{Python query object to use in query and to update.}
}
\value{
Returns the update Python query object.
}
\examples{
## The function is currently defined as
function (theQueryPyObj)
{
    theQueryPyObj$update_from_selected()
    return(theQueryPyObj)
  }
}

