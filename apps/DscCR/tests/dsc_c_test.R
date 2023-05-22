#DscCR Copyright 2023 University of Texas MD Anderson Cancer Center
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(DscCR)


inputDir <- getTestInputDirDscCR()
# outputDir <- getTestOutputDirDscCR()

myScores <- file.path(inputDir, "PCAScores.tsv")
myBatches <- file.path(inputDir, "BatchData.tsv")

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
  # do not unlink or create, is not used
  # also, do not want to delete results from compare_dsc_times
  # unlink(outputDir, recursive=TRUE, force=TRUE)
  # dir.create(outputDir, showWarnings=FALSE, recursive=TRUE)

  # build the data to test
  pcascores <- t(read.table(myScores, row.names=1, sep="\t", na.strings="NA", header=TRUE, stringsAsFactors=FALSE, quote="", check.names=FALSE))
	batchInfo <- read.table(myBatches, sep="\t", na.strings="NA", header=TRUE, stringsAsFactors=FALSE, quote="", colClasses=c("character"), as.is=TRUE, check.names=FALSE)
	indexList <- !is.na(match(batchInfo$Sample, rownames(pcascores)))
	batchIdsForSamples <- batchInfo$BatchId[indexList]
	myDSCPermutations <- 1000
	myDSCThreads <- 5
	dscAllResults <- NULL
	time <- system.time(dscAllResults <- pvalueDSC_C(pcascores, batchIdsForSamples, myDSCPermutations, 0, 0, myDSCThreads))
	message("mPvalue=",dscAllResults@mPvalue, " should be ", "0")
	message("mDSC=",dscAllResults@mDSC, " should be ", "0.523")
	message("mDB=",dscAllResults@mDB, " should be ", "33.837")
	message("mDW=",dscAllResults@mDW, " should be ", "64.659")
	message(time)
	matches = TRUE
	if (!startsWith(toString(dscAllResults@mPvalue), "0"))
	{
	  matches = FALSE
	  message("pValue does not match")
	}
	if (!startsWith(toString(dscAllResults@mDSC), "0.523"))
	{
	  matches = FALSE
	  message("DSC does not match")
	}
	if (!startsWith(toString(dscAllResults@mDB), "33.837"))
	{
	  matches = FALSE
	  message("DB does not match")
	}
	if (!startsWith(toString(dscAllResults@mDW), "64.659"))
	{
	  matches = FALSE
	  message("DW does not match")
	}
	# return TRUE or FALSE for testing
	matches
}
