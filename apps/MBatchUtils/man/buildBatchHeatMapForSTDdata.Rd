\name{buildBatchHeatMapForSTDdata}
\alias{buildBatchHeatMapForSTDdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{buildBatchHeatMapForSTDdata}
\description{Build an NGCHM for the given Standardized Data dataset.}
\usage{
buildBatchHeatMapForSTDdata(theFile, theBatchType, theOutputDir,
theShaidyMapGen,
theShaidyMapGenJava = "/usr/bin/java")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theFile}{Full path to the matrix_data.tsv file from Standardized Data.
  There should be a batches.tsv in the same directory.}
  \item{theBatchType}{Batch Type for which to create the NGCHM.
  (Should be a column in batches.tsv.)}
  \item{theOutputDir}{Output directory for NGCHM.}
  \item{theShaidyMapGen}{Full path to ShaidyMapGen.jar.
  Required by NGCHM package.}
  \item{theShaidyMapGenJava}{Full path to Java. Required by NGCHM package.}
}
\details{
Build an NGCHM named <theBatchType>_matrix.ngchm
}
\value{Does not return a meaningful value.}
