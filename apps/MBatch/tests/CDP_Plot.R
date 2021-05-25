# MBatch Copyright (c) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021 University of Texas MD Anderson Cancer Center
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
# MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>

require(MBatch)

inputDir <- getTestInputDir()
outputDir <- getTestOutputDir()
compareDir <- getTestCompareDir()

theGeneFile1=file.path(inputDir, "CDP_norep_data1.tsv")
theGeneFile2=file.path(inputDir, "CDP_norep_data2.tsv")
theOutputDir=file.path(outputDir, "CDP_Plot")
theCompareFile=file.path(compareDir, "CDP_Plot.tsv")

theRandomSeed=314

if (!is.null(inputDir))
{
  warnLevel<-getOption("warn")
  on.exit(options(warn=warnLevel))
  # warnings are errors
  options(warn=3)
  # if there is a warning, show the calls leading up to it
  options(showWarnCalls=TRUE)
  # if there is an error, show the calls leading up to it
  options(showErrorCalls=TRUE)
  #
  unlink(theOutputDir, recursive=TRUE, force=TRUE)
  dir.create(theOutputDir, showWarnings=FALSE, recursive=TRUE)
  # load data
  theData1 <- readAsGenericMatrix(theGeneFile1)
  theData2 <- readAsGenericMatrix(theGeneFile2)
  ##############################################################################
  theUseReplicatesUnpaired <- FALSE
  theUnmatchedCount <- 1000
  # get list of natural replicates
  pairedSamples <- colnames(theData1)[colnames(theData1) %in% colnames(theData2)]
  # get list of unmatched replicates
  unpairedSamples1 <- NULL
  unpairedSamples2 <- NULL
  if (TRUE==theUseReplicatesUnpaired)
  {
    unpairedSamples1 <- sample(colnames(theData1), theUnmatchedCount, replace=TRUE)
    unpairedSamples2 <- sample(colnames(theData2), theUnmatchedCount, replace=TRUE)
  }
  else
  {
    data1samples <- colnames(theData1)[!colnames(theData1) %in% pairedSamples]
    data2samples <- colnames(theData2)[!colnames(theData2) %in% pairedSamples]
    if (is.null(pairedSamples))
    {
      data1samples <- colnames(theData1)
      data2samples <- colnames(theData2)
    }
    if ((0==length(data1samples))||(0==length(data2samples)))
    {
      unpairedSamples1 <- c()
      unpairedSamples2 <- c()
    }
    else
    {
      unpairedSamples1 <- sample(data1samples, theUnmatchedCount, replace=TRUE)
      unpairedSamples2 <- sample(data2samples, theUnmatchedCount, replace=TRUE)
    }
  }
  CDP_Plot(file.path(theOutputDir, "CDP_Plot.png"), theData1, theData2, pairedSamples, pairedSamples, unpairedSamples1, unpairedSamples2,
           theSubTitle="no replicates", theMethod="pearson", theUse="pairwise.complete.obs", theSeed=theRandomSeed,
           theLinePlot=TRUE, theHistPlot=TRUE, theBinWidth=NULL)
  message("No error means test was OK.")
  TRUE
} else {
  message("No test data. Skip test.")
  TRUE
}
