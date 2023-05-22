#DscCR Copyright 2023 University of Texas MD Anderson Cancer Center
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(DscCR)


inputDir <- getTestInputDirDscCR()
outputDir <- getTestOutputDirDscCR()

myScores <- file.path(inputDir, "PCAScores.tsv")
myBatches <- file.path(inputDir, "BatchData.tsv")
timefile <- file.path(outputDir, "times.tsv")

message(myScores)
message(myBatches)
message(timefile)

theRandomSeed=314


writeTimeFile <- function(theFile, theCTime, theRTime)
{
  message("writeTimeFile theFile=",theFile)
  myFile <- file(theFile, "w+")
  on.exit(close(myFile))
  # ELAPSED is time perceived by user/waiting process
  # USER + SYSTEM for single threaded or C code is system resource time spent
  # CHILD_USER + CHILD_SYSTEM for parallel R/Python is system resource time spent
  cat("PROCESS\tUSER\tSYSTEM\tELAPSED\tCHILD_USER\tCHILD_SYSTEM\n", file=myFile, append=TRUE)
  cat(paste("CTIME\t", paste(theCTime, collapse = "\t"), "\n", sep="", collapse=""), file=myFile, append=TRUE)
  cat(paste("RTIME\t", paste(theRTime, collapse = "\t"), "\n", sep="", collapse=""), file=myFile, append=TRUE)
}

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
  # other tests should not unlink or create
  # do not want to delete results from compare_dsc_times
  unlink(outputDir, recursive=TRUE, force=TRUE)
  dir.create(outputDir, showWarnings=FALSE, recursive=TRUE)

  # build the data to test
  pcascores <- t(read.table(myScores, row.names=1, sep="\t", na.strings="NA", header=TRUE, stringsAsFactors=FALSE, quote="", check.names=FALSE))
	batchInfo <- read.table(myBatches, sep="\t", na.strings="NA", header=TRUE, stringsAsFactors=FALSE, quote="", colClasses=c("character"), as.is=TRUE, check.names=FALSE)
	indexList <- !is.na(match(batchInfo$Sample, rownames(pcascores)))
	batchIdsForSamples <- batchInfo$BatchId[indexList]
	myDSCPermutations <- 1000
	myDSCThreads <- 5
	dscAllResults <- NULL
	# call C code
	# TODO: add passing in SEED
	timeC <- system.time(dscAllResults <- pvalueDSC_C(pcascores, batchIdsForSamples, myDSCPermutations, 0, 0, myDSCThreads))
	message(paste(timeC, collapse = "\t"))
	# call R/Python code
	timeRP <- system.time(dscAllResults <- pvalueDSC_R(pcascores, batchIdsForSamples, myDSCPermutations, 0, 0, myDSCThreads, theRandomSeed))
	message(paste(timeRP, collapse = "\t"))
	# write results to timefile
	print(Sys.getenv())
	message("compare_dsc_times outputDir=",outputDir)
	writeTimeFile(timefile, timeC, timeRP)
	message(myScores)
	message(myBatches)
	message(timefile)
	print(timeC)
	message(paste(timeC, collapse = "\t"))
	print(timeRP)
	message(paste(timeRP, collapse = "\t"))
	# return TRUE or FALSE for testing
	TRUE
}


