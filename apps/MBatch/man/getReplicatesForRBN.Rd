\name{getReplicatesForRBN}
\alias{getReplicatesForRBN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getReplicatesForRBN}
\description{Function to find exact replicates in sample (column) names.}
\usage{
getReplicatesForRBN(matrix1, matrix2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix1}{First matrix, with sample ids in colnames() and
  feature names (like genes) in rownames().}
  \item{matrix2}{Second matrix, with sample ids in colnames() and
  feature names (like genes) in rownames().}
}
\details{
Find colnames that match between matrixes.
}
\value{
A vector of matching feature ids.
}
