\name{mbatchRunFromConfig}
\alias{mbatchRunFromConfig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mbatchRunFromConfig}
\description{Run MBatch from a configuration file. This is for
reproduciblity or for use with Batch Effects Interface}
\usage{
mbatchRunFromConfig(theConfigFile, theOutputDir, theNaStrings,
theShaidyMapGen, theShaidyMapGenJava,
                 theNGCHMShaidyMem = "16G", thePCAMem = "4800m",
                 theBoxplotMem = "16G", theRunPostFlag=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theConfigFile}{Full path to the configuration file.
  Documentation will be in the MBatch release documentation.}
  \item{theOutputDir}{Full path to directory to create output, and the
  function attempts to create the directory.}
  \item{theNaStrings}{This is the argument to mbatchLoadFiles,
  and should generally be "NA".}
  \item{theShaidyMapGen}{Full path to ShaidyMapGen.jar.
  Required by NGCHM package.}
  \item{theShaidyMapGenJava}{Full path to Java. Required by NGCHM package.}
  \item{theNGCHMShaidyMem}{16G default}
  \item{thePCAMem}{4800m default.}
  \item{theBoxplotMem}{16G default.}
  \item{theRunPostFlag}{FALSE default. If TRUE, build DSC Overview file.}
}
\details{
Does not return any meaningful value, but creates an empty file
MBATCH_SUCCESS.txt in the output directory, if successful.
}
\value{
Does not return any meaningful value, but creates an empty file
MBATCH_SUCCESS.txt in the output directory, if successful.
}
