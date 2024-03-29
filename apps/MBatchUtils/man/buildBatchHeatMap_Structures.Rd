\name{buildBatchHeatMap_Structures}
\alias{buildBatchHeatMap_Structures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{buildBatchHeatMap_Structures}
\description{Build an NGCHM for the given data structures.}
\usage{
buildBatchHeatMap_Structures(theMatrixData, theBatchData, theTitle,
theOutputFile, theSortByType, theRowType = "scholar",
theColType = "bio.tcga.barcode.sample", theRowCluster = NULL,
theColCluster = NULL,
theShaidyMapGen, theNgchmWidgetJs,
theShaidyMapGenJava = "/usr/bin/java", theShaidyMapGenArgs = "-Xmx16G",
theEncoding="en_US.UTF-8")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theMatrixData}{A matrix with samples in columns and features in rows. See the MBatch data structures docs for more information.}
  \item{theBatchData}{A dataframe with the first column being Sample and the rest being batch types. See the MBatch data structures docs for more information.}
  \item{theTitle}{Title to use for the NGCHM.}
  \item{theOutputFile}{Full path to NGCHM file name.}
  \item{theSortByType}{Batch Type (from theBatchType) by which to sort the samples.}
  \item{theRowType}{See NGCHM docs for details, but since many Standardized Data features are customized, usually use "scholar".}
  \item{theColType}{See NGCHM docs for details, but since Standardized Data uses TCGA barcodes, usually use "bio.tcga.barcode.sample".}
  \item{theRowCluster}{If NULL, do not cluster. Otherwise, this should be a vector of string, with the first entry being the method argument to distanceMatrix (like "pearson") and the second entry being the method argument to hclust (like "ward.D2").
  }
  \item{theColCluster}{If NULL, do not cluster. Otherwise, this should be a vector of string, with the first entry being the method argument to distanceMatrix (like "pearson") and the second entry being the method argument to hclust (like "ward.D2").}
  \item{theShaidyMapGen}{Full path to ShaidyMapGen.jar. Required by NGCHM package.}
  \item{theNgchmWidgetJs}{Full path to ngchmWidget-min.js. Required by NGCHM package.}
  \item{theShaidyMapGenJava}{Full path to Java. Required by NGCHM package.}
  \item{theShaidyMapGenArgs}{Java run-time arguments, like -Xmx16G".}
  \item{theEncoding}{Text encoding, defaults to "en_US.UTF-8". Other encodings may not work properly.}
}
\details{
Build an NGCHM from data structures.
}
\value{
Returns value from NGCHM function chmExportToFile.
}
