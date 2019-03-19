\name{buildBatchHeatMap_Files}
\alias{buildBatchHeatMap_Files}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{buildBatchHeatMap_Files}
\description{Build an NGCHM for the given files.}
\usage{
buildBatchHeatMap_Files(theMatrixFile, theBatchFile, theTitle,
theOutputFile, theSortByType, theRowType = "labels",
theColType = "bio.tcga.barcode.sample", theRowCluster = NULL,
theColCluster = NULL,
theShaidyMapGen,
theShaidyMapGenJava = "/usr/bin/java", theShaidyMapGenArgs = "-Xmx16G")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theMatrixFile}{Full path to the matrix_data.tsv formatted file. See the MBatch data format docs for more details.}
  \item{theBatchFile}{Full path to the batches.tsv formatted file. See the MBatch data format docs for more details.}
  \item{theTitle}{Title to use for the NGCHM.}
  \item{theOutputFile}{Full path to NGCHM file name.}
  \item{theSortByType}{Batch Type (from theBatchType) by which to sort the samples.}
  \item{theRowType}{See NGCHM docs for details, but since many Standardized Data features are customized, usually use "labels".}
  \item{theColType}{See NGCHM docs for details, but since Standardized Data uses TCGA barcodes, usually use "bio.tcga.barcode.sample".}
  \item{theRowCluster}{If NULL, do not cluster. Otherwise, this should be a vector of string, with the first entry being the method argument to distanceMatrix (like "pearson") and the second entry being the method argument to hclust (like "ward.D2").
  }
  \item{theColCluster}{If NULL, do not cluster. Otherwise, this should be a vector of string, with the first entry being the method argument to distanceMatrix (like "pearson") and the second entry being the method argument to hclust (like "ward.D2").}
  \item{theShaidyMapGen}{Full path to ShaidyMapGen.jar. Required by NGCHM package.}
  \item{theShaidyMapGenJava}{Full path to Java. Required by NGCHM package.}
  \item{theShaidyMapGenArgs}{Java run-time arguments, like -Xmx16G".}
}
\details{
Build an NGCHM from files.
}
\value{
Returns value from NGCHM function chmExportToFile.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
