\name{htmlMutationBatchEffects}
\alias{htmlMutationBatchEffects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{htmlMutationBatchEffects}
\description{Build HTML index files for each Mutation Batch output, to allow creation of a PDF "report".}
\usage{
htmlMutationBatchEffects(theBaseDir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theBaseDir}{Base directory to start looking for assessment results that need report files.}
}
\details{
Looks for callReference.tsv files and creates index.html files in the same directory for the data.
}
\value{
Does not return a meaningful value.
}
