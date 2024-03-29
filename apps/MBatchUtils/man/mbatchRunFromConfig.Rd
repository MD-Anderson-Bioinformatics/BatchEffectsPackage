\name{mbatchRunFromConfig}
\alias{mbatchRunFromConfig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mbatchRunFromConfig}
\description{Run MBatch from a configuration file. This is for
reproduciblity or for use with Batch Effects Interface}
\usage{
mbatchRunFromConfig(theConfigFile, theMatrixFile, theBatchesFile,
  theZipDataDir, theZipResultsDir, theNaStrings,
  theShaidyMapGen, theNgchmWidgetJs, theShaidyMapGenJava,
  theNGCHMShaidyMem="16G", theRunPostFlag=FALSE,
  theMatrixFile2=NULL, theBatchesFile2=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theConfigFile}{Full path to the configuration file.
  Documentation will be in the MBatch release documentation.}
  \item{theMatrixFile}{Full path to data file.}
  \item{theBatchesFile}{Full path to batch file.}
  \item{theZipDataDir}{Full path to directory for DATA directory.}
  \item{theZipResultsDir}{Full path to directory for RESULTS directory.}
  \item{theNaStrings}{This is the argument to mbatchLoadFiles,
  and should generally be "NA".}
  \item{theShaidyMapGen}{Full path to ShaidyMapGen.jar.
  Required by NGCHM package.}
  \item{theNgchmWidgetJs}{Full path to ngchmWidget-min.js.
  Required by NGCHM package.}
  \item{theShaidyMapGenJava}{Full path to Java. Required by NGCHM package.}
  \item{theNGCHMShaidyMem}{16G default}
  \item{theRunPostFlag}{FALSE default. If TRUE, build DSC Overview file.}
  \item{theMatrixFile2}{Full path to secondary data file.}
  \item{theBatchesFile2}{Full path to secondary batch file.}
}
\details{
Does not return any meaningful value, but creates an empty file
MBATCH_SUCCESS.txt in the output directory, if successful.
}
\value{
Does not return any meaningful value, but creates an empty file
MBATCH_SUCCESS.txt in the output directory, if successful.
}
