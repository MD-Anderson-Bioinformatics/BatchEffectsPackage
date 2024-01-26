\name{getDapiQuery}
\alias{getDapiQuery}
\title{getDapiQuery}
\description{
Creates the basic Python query object and returns it after updating it via the URL.
}
\usage{getDapiQuery(theUrl)}
\arguments{
  \item{theUrl}{string giving URL to web app to query for data. Should not have trailing slash.}
}
\details{
Creates the basic Python query object and returns it after updating it via the URL.
}
\value{
This is a Python object with the following members.

The available members are vectors of strings valid to use in the "selected" members, as search options. Note that you can build searches that return no results. See the website at https://bioinformatics.mdanderson.org/MQA/ for more details.
        self.available_jobtype = []
        self.available_platforms = []
        self.available_categories = []
        self.available_projects = []
        self.available_details = []
        self.available_dataversions = []
        self.available_testversions = []
        self.available_program = []
        self.available_data = []
        self.available_sources = []

The selected members are vectirs if strings describing what you are searching for.
        self.selected_jobtype = []
        self.selected_platforms = []
        self.selected_categories = []
        self.selected_projects = []
        self.selected_details = []
        self.selected_dataversions = []
        self.selected_testversions = []
        self.selected_program = []
        self.selected_data = []
        self.selected_sources = []

The available_datasets member gives a list of datasets that match your search request.
        self.available_datasets = []
}
\examples{
## The function is currently defined as
function (theUrl)
{
    dapiQueryImport <- import("mbatch.dapi.query")
    queryPyObj <- dapiQueryImport$DapiQuery(theUrl)
    queryPyObj$update_from_selected()
    return(queryPyObj)
  }
}

