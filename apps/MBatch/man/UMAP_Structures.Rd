\name{UMAP_Structures}
\alias{UMAP_Structures}
\title{UMAP_Structures}
\description{Generate UMAP image and data files}
\usage{
UMAP_Structures(theData, theTitle, theOutputDir, theDataVersion,
theTestVersion, theBatchTypeAndValuePairsToRemove=list(),
theBatchTypeAndValuePairsToKeep=list(),
theDoDSCFlag=FALSE, theDSCPermutations=0,
theDSCThreads=0, theDoDscPermsFileFlag=FALSE, theSeed=314)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theData}{An instance of BEA_DATA \code{\link{BEA_DATA-class}}}
  \item{theTitle}{Object of class \code{"character"} Title to use in PNG files. }
  \item{theOutputDir}{Location for output.}
	\item{theDataVersion}{Object of class \code{"character"} Version of the data prefaced
	with DATA_, such as DATA_2022-09-09-1600. Use empty string to not include. }
	\item{theTestVersion}{Object of class \code{"character"} Version of the test prefaced
	with TEST_, such as TEST_2022-10-10-1300. Use empty string to not include. }
  \item{theBatchTypeAndValuePairsToRemove}{Object of class \code{"list"} A list of vectors
  containing the batch type (or * for all types) and the value to remove. list() indicates
  none while NULL will cause an error. }
  \item{theBatchTypeAndValuePairsToKeep}{Object of class \code{"list"} A list of vectors
  containing the batch type (or * for all types) and a vector of the the value(s) to keep.
  list() indicates none while NULL will cause an error. }
  \item{theDoDSCFlag}{Object of class \code{"logical"}
  A flag indicating whether or not to
  perform the DSC calculation for PCA-Plus plots. }
  \item{theDSCPermutations}{Object of class \code{"numeric"}
  The number of permutations to perform for DSC calculations.
  Defaults to 2000}
  \item{theDSCThreads}{Object of class \code{"numeric"}
  The number of threads (or CPUs or cores)
  to use to perform DSC calculations.
  Defaults to 1}
  \item{theDoDscPermsFileFlag}{Object of class \code{"logical"}
  A flag indicating whether
  or not to perform the DSC calculation for PCA-Plus plots. }
  \item{theSeed}{Default to 314.}
}
\details{Generate UMAP output and data files}
\value{Vector of strings being full paths to output files.}
