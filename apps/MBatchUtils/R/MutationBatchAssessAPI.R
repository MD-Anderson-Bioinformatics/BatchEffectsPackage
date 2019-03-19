#MBatchUtils Copyright ? 2018 University of Texas MD Anderson Cancer Center
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(dunn.test)
library(MBatch)
#source("GSColors.R")
#source("MutationBatchExtractAPI.R")

# used internally for marking ShipDate for MBatch
isTrendBatch<-function(theBatchTypeName, theListOfBatchIds)
{
	return(is.element(theBatchTypeName, c("ShipDate")))
}

mutationBatchAssess <- function(theTypeCountDir, theOutputDir, theJavaArgs="-Xms24000m", theThreads=5,
																thePvalueCutoff=.00001, theZScoreCutoff=1.96, thePCAflag=FALSE,
																theBatchTypes=c("BatchId", "PlateId", "ShipDate", "TSS"),
																theMutationTypes=NULL)
{
	# START HERE
	# theTypeCountDir - full path to top level of directory containing disease sub-directories and files created by mutationBatchExtract
	# theOutputDir - full path to output directory for this function
	# theJavaArgs="-Xms24000m" - argument to Java for Boxplot and DSC computations.
	# theThreads=5 - argument to Java for # cores (threads) to use for DSC computations
	# thePvalueCutoff=.00001 - P-value to use for cutoff for "significant" difference in Kruskal-Wallis and Dunn's tests
	# theZScoreCutoff=1.96 - cut off for Dunn's test for "significant" confidence interval
	# 			95% Confidence Internal (1.96 Standard Deviations).
	# 			99% Confidence Internal (2.576 Standard Deviations).
	# 			99.9999% (4.892 Standard Deviations)
	# thePCAflag=FALSE - TRUE means to PCA with DSC (which is slow), FALSE means don't do PCA
	# theBatchTypes=c("BatchId", "PlateId", "ShipDate", "TSS") - Batch types to check (correspond to columns for batch files)
	#
	# This function processes the directories and files created by mutationBatchExtract
	# and does batch assessment processing on mutation data. Batch comparisons are both
	# between platforms, within platforms, and between disease types.
	# Output is divided into <BatchType>_<MutationType> subdirectories. Datasets without Kruskal-Wallis-called significant results
	# will ony have a PNG (FullMutCounts_<BatchType>_<MutationType>_Diagram.PNG) giving the Kruskal-Wallis/Dunn's Test results.
	# For datasets with significant calls, the aforementioned PNG will exist, as will files using the patterns
	# FullMutCounts_<BatchType>_<MutationType>_<DiseaseType>_MutDots_Diagram.PNG
	#			which contains all samples sorted by batch with mutations for all platforms plotted
	# NarrowBoxplot_<BatchType>_<MutationType>_<DiseaseType>_<Log10|ZScore>_<Platform_GenomeReference_MutationType>_Diagram.PNG
	#			which contains mutations for the given platform plotted by batches (with more than the average number of samples)
	#			and plots are provided both in terms of Log10 and ZScores.
	# WideBoxplot_<BatchType>_<MutationType>_<DiseaseType>_<Log10|ZScore>_Diagram.PNG
	#			which contains mutations for all platforms plotted by batches and plots are provided both in terms of Log10 and ZScores.
	# callReference.tsv
	#			is a tab delimited file with a header of MutationType, BatchTuype, MutationFile, and called Batches.
	#     Batches are comma delimited within parenthesis.
	message(mbatchUtilVersion())
	message("mutationBatchAssess")
	message("find matrix files")
	# make sure output directory exists
	dir.create(theOutputDir, showWarnings=FALSE, recursive = TRUE)
	# process file names to get mutation types that are available
	dataTypes <- theMutationTypes
	if (is.null(dataTypes))
	{
		dataTypes <- collectCountTypes(theTypeCountDir)
	}
	print(dataTypes)
	for(myCountType in dataTypes)
	{
		# each mutation type
		message("Find all data files for ", myCountType)
		# find all the data files for that mutation type
		matrixFiles <- list.files(theTypeCountDir, pattern=paste(".*", myCountType, ".tsv", sep=""), full.names=TRUE, recursive=TRUE)
		message("found files ", length(matrixFiles))
		# list (not vector) of Kruskal results
		# list of sub-list. The sub-list has names() corresponding to batch types
		# the values for each batch type will either be a single p-value or
		# a sub-list with the p-value and a vector of batches found as significant by Dunn's test
		resultList <- list()
		# vector of titles for each set of p-value results
		titleVector <- c()
		# checkpoint is to complete all entries for a batch type
		dirs <- paste(theBatchTypes, myCountType, sep="_")
		if (!all(dir.exists(file.path(theOutputDir, dirs))))
		{
			for (myDataFile in matrixFiles)
			{
				# for each data file in this set of mutation type files
				message("Do Kruskal calculations for: ", myDataFile)
				# get the corresponding batch file (see MutationBatchExtractAPI.R for naming conventions)
				myBatchFile <- file.path(dirname(myDataFile), gsub(x=basename(myDataFile), pattern=paste("HG.*", myCountType, sep=""), replacement="batches"))
				# build the title based on file names
				myTitle <- gsub(x=basename(myDataFile), pattern=paste("_", myCountType, ".tsv", sep=""), fixed=TRUE, replacement="")
				# get vector of p-values (and possible Dunn's batch values) for this data (see resultList description for details)
				pValueVector <- loadAndRunForKruskal(myDataFile, myBatchFile, theBatchTypes, thePvalueCutoff, theZScoreCutoff)
				# add p-values to result list
				resultList <- c(resultList, list(pValueVector))
				# add title to title vector
				titleVector <- c(titleVector, myTitle)
			}
			message("get batch files")
			# get vector of all batch files
			batchFiles <- list.files(theTypeCountDir, pattern=".*_batches.tsv", full.names=TRUE, recursive=TRUE)
			for (myBatchType in theBatchTypes)
			{
				if(!dir.exists(file.path(theOutputDir, paste(myBatchType, myCountType, sep="_"))))
				{
					# for each batch type
					message("Do Kruskal output for: ", myBatchType)
					# extract and convert p-values into logged p-values (-log10) for the given batch type
					loggedPvalues <- collectLoggedPvalues(resultList, myBatchType)
					# if the values for a batch type has batch names (in addition to the p-value)
					# collect the batch names found as significant by a Dunn's test
					collectBatches <- collectLoggedBatches(resultList, myBatchType)
					# subdir named for batch type and count type
					outdir <- file.path(theOutputDir, paste(myBatchType, "_", myCountType, sep=""))
					# make sure sub-dir exists
					dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
					# output for Kruval - name the file "FullMutCounts_<batch-type>_<mutation-type>_Diagram.PNG"
					outfile <- file.path(outdir, paste("FullMutCounts_", myBatchType, "_", myCountType, "_Diagram.PNG", sep=""))
					outfileTSV <- file.path(outdir, paste("FullMutCounts_", myBatchType, "_", myCountType, ".tsv", sep=""))
					# add names from title vector to logged p-values
					names(loggedPvalues) <- titleVector
					# if Dunn's test returned batch names, put them into parens to put in diagram otherwise, put an empty string
					batchAddOns <- getBatchAddOn(collectBatches)
					# write the PNG for Kruskal's output
					makeKruskalPng(outfile, loggedPvalues, collectBatches, paste(myBatchType, "for", myCountType, sep=" "), batchAddOns, thePvalueCutoff)
					makeKruskalTsv(outfileTSV, loggedPvalues, collectBatches, paste(myBatchType, "for", myCountType, sep=" "), batchAddOns, thePvalueCutoff)
					# if any batches were called by Dunn't test, get the names, and write strip charts and MBatch for their batch types
					if (sum(!(batchAddOns==""))>0)
					{
						# get disease for which batches were called by Dunn's
						mySBdiseases <- extractDiseasesFromNames(titleVector[!(batchAddOns=="")])
						# find matrix (data) files for diseases with batches called by Dunn's
						matrixFilesSB <- matrixFiles[grepl(paste(".*", paste(mySBdiseases, collapse="|", sep="") ,".*", sep=""), matrixFiles)]
						# find batch files for diseases with batches called by Dunn's
						batchFilesSB  <- batchFiles[grepl(paste(".*", paste(mySBdiseases, collapse="|", sep="") ,".*", sep=""), batchFiles)]
						message("Do combined dotplots for: ", myBatchType, " and ", myCountType)
						# new subdir based on batch type and mutation type
						outdirSandB <- file.path(theOutputDir, paste(myBatchType, "_", myCountType, sep=""))
						# make/ensure subdir exists
						dir.create(outdirSandB, showWarnings=FALSE, recursive=TRUE)
						# collect output for reference
						refLines <- c()
						for(myIndex in 1:length(batchAddOns))
						{
							refBatchAddOn <- batchAddOns[myIndex]
							if (""!=refBatchAddOn)
							{
								refMutationFile <- names(loggedPvalues)[myIndex]
								refLines <- c(refLines, paste(myCountType, myBatchType, refMutationFile, refBatchAddOn, sep="\t"))
							}
						}
						# write output for reference
						if (length(refLines)>0)
						{
							writeReferenceOutput(refLines, outdirSandB)
						}
						# output for dot plot and MBatch
						#TODO: Removed until revisit for better plots
						#writeStripchartsAndMBatch(outdirSandB, matrixFilesSB, myBatchType, myCountType, batchFilesSB, theJavaArgs, theThreads, thePCAflag)
					}
				}
			}
		}
	}
}

writeReferenceOutput <- function(theLines, theOutdir)
{
	filename <- file.path(theOutdir, "callReference.tsv")
	writeLines(c("MutationType\tBatchType\tMutationFile\tBatches", theLines), con=filename)
}

extractDiseasesFromNames <- function(theMatrixFiles)
{
	# get disease for which batches were called by Dunn's
	diseases <- c()
	for (myFile in theMatrixFiles)
	{
		#prefix will be "TCGA-ACC" or similar patter for GDC and "acc" for TCGA
		diseases <- unique(sort(c(diseases, strsplit(basename(myFile), split=".", fixed=TRUE)[[1]][1])))
	}
	diseases
}

writeStripchartsAndMBatch <- function(theOutdir, theMatrixFiles, theBatchType, theCountType, theBatchFiles,
																				theJavaArgs, theThreads, thePCAflag)
{
	message("get diseases")
	# pull diseases from file names
	diseases <- extractDiseasesFromNames(theMatrixFiles)
	for (myDis in diseases)
	{
		message("process disease ", myDis, " for ", theBatchType, " in ", theCountType)
		matrixFiles <- theMatrixFiles[startsWith(basename(theMatrixFiles), myDis)]
		message("get max value for mutations")
		maxValue <- findMaxValue(matrixFiles)
		message("collect all sample ids")
		allSampleIds <- collectNonZeroSampleIds(matrixFiles)
		if (length(allSampleIds)>0)
		{
			# build batch dataframe for all samples (includes those with Unknown batch values form allSampleIds) and remove unused batch data
			message("batch data")
			# build combined dataframe for batch files
			allBatchesDataframe <- buildCombinedDataframe(theBatchFiles, theCountType, allSampleIds)
			# sort batch dataframe by samples and then...
			allBatchesDataframe <- allBatchesDataframe[order(allBatchesDataframe$Sample),]
			# sort batch dataframe by selected batch type batch values
			allBatchesDataframe <- allBatchesDataframe[order(allBatchesDataframe[theBatchType]),]
			# overwrite allSampleIds with dataframe order
			allSampleIds <- allBatchesDataframe$Sample
			# size checking
			message("allSampleIds")
			print(length(allSampleIds))
			message("allBatchesDataframe")
			print(dim(allBatchesDataframe))
			if (length(unique(sort(as.vector(unlist(allBatchesDataframe[theBatchType])))))>1)
			{
				#message("write Mbatch")
				#for (myMatrix in matrixFiles)
				#{
					# make output directory based on matrix file name
					#myOutdir <- file.path(theOutdir, gsub(x=basename(myMatrix), pattern=".tsv", replacement="", fixed=TRUE))
					# make output direcotry
					#dir.create(myOutdir, showWarnings=FALSE, recursive=TRUE)
					# call and plot MBatch data
					# TODO: finish plotMBatch docs
					#plotMBatch(myOutdir, myMatrix, allBatchesDataframe, paste(myDis, theBatchType, theCountType, sep=" - "),
					#					 theJavaArgs, theThreads, thePCAflag, theBatchType)
				#}
				message("write strips")
				# TODO: finish writeMutationSamplesDotPlot docs
				message("call writeMutationSamplesDotPlot")
				writeMutationSamplesDotPlot(theOutdir, matrixFiles, allBatchesDataframe, allSampleIds, theBatchType, theCountType, myDis, maxValue)
				# TODO: finish writeMutationStripchart docs
				message("call writeMutationBoxplot_Wide 1")
				writeMutationBoxplot_Wide(theOutdir, matrixFiles, allBatchesDataframe, allSampleIds, theBatchType, theCountType, myDis, maxValue, theZscoreFlag=FALSE)
				message("call writeMutationBoxplot_Wide 2")
				writeMutationBoxplot_Wide(theOutdir, matrixFiles, allBatchesDataframe, allSampleIds, theBatchType, theCountType, myDis, maxValue, theZscoreFlag=TRUE)
				message("call writeMutationBoxplot_Narrow 1")
				writeMutationBoxplot_Narrow(theOutdir, matrixFiles, allBatchesDataframe, allSampleIds, theBatchType, theCountType, myDis, maxValue, theZscoreFlag=FALSE)
				message("call writeMutationBoxplot_Narrow 2")
				writeMutationBoxplot_Narrow(theOutdir, matrixFiles, allBatchesDataframe, allSampleIds, theBatchType, theCountType, myDis, maxValue, theZscoreFlag=TRUE)
				#TDC#writeMutationStripchart(theOutdir, matrixFiles, allBatchesDataframe, allSampleIds, theBatchType, theCountType, myDis, maxValue)
			}
		}
	}
}

collectLoggedPvalues <- function(theResultList, theBatchType)
{
	# extract and convert p-values into logged p-values (-log10) for the given batch type
	pvalues <- lapply(theResultList, function(thePvalues)
	{
		# If there are Dunn's test results, strip the batches out and return the p-value only
		if (length(thePvalues[theBatchType][[1]])>1)
		{
			return(thePvalues[theBatchType][[1]][1])
		}
		else
		{
			# without Dunn's test results, just use the p-values
			return(thePvalues[theBatchType])
		}
	})
	# -log10 all values
	pvalues <- -log10(as.vector(unlist(pvalues)))
	pvalues
}

writeMutationSamplesDotPlot <- function(theOutdir, theMatrixFiles, theAllBatchesDataframe, theAllSampleIds,
																				theBatchType, theCountType, theDisease, theMax, theZscoreFlag=FALSE)
{
	message("writeMutationSamplesDotPlot")
	outfile <- file.path(theOutdir, paste("FullMutCounts_", theBatchType, "_", theCountType, "_", theDisease, "_MutDots_Diagram.PNG", sep=""))
	message("write ", basename(outfile))
	symbol <- 0
	# build matrix of values
	bigDataVector <- c()
	for (myDataFile in theMatrixFiles)
	{
		dataMatrix <- readMatrixScore(myDataFile)
		fileDataVector <- rowSums(t(dataMatrix), na.rm=TRUE)
		plotDataVector <- rep(0, length.out=length(theAllSampleIds))
		names(plotDataVector) <- theAllSampleIds
		plotDataVector[names(fileDataVector)] <- fileDataVector
		bigDataVector <- c(bigDataVector, as.vector(unlist(plotDataVector)))
	}
	myBigDataMatrix <- matrix(data=as.vector(unlist(bigDataVector)),
														nrow=length(theAllSampleIds),
														ncol=length(theMatrixFiles),
														dimnames=list(theAllSampleIds, theMatrixFiles),
														byrow=TRUE)
	##############################################################################
	if (isTRUE(theZscoreFlag))
	{
		# z-score max value
		theMax <- max(as.vector(unlist(apply(myBigDataMatrix, 2, function(theColumn)
		{
			max((abs(theColumn - mean(theColumn))) / sd(theColumn))
		}))), na.rm = TRUE)
	}
	else
	{
		# log10(+1) scores
		myBigDataMatrix <- log10(myBigDataMatrix+1)
		theMax <- log10(theMax+1)
	}
	##############################################################################
	myWidth <- length(rownames(myBigDataMatrix))*15
	if (myWidth<100)
	{
	  myWidth <- 100
	}
	if (myWidth>5000)
	{
	  myWidth <- 5000
	}
	#batchesAsLevels <- as.factor(as.vector(unlist(theAllBatchesDataframe[theBatchType]))[keepVector])
	batchesAsLevels <- as.factor(as.vector(unlist(theAllBatchesDataframe[theBatchType])))
	colors <- getColorListAny(length(unique(sort(batchesAsLevels))))
	######
	message("writeMutationSamplesDotPlot outfile=", outfile)
	message("writeMutationSamplesDotPlot myWidth=", myWidth)
	CairoPNG(filename=outfile, width = myWidth, height = 1000, pointsize=12)
	on.exit(dev.off(), add = TRUE)
	# c(bottom, left, top, right)
	par(mar=c(30,10,2,2))
	for(dataSet in colnames(myBigDataMatrix))
	{
		thisData <- myBigDataMatrix[, dataSet]
		##############################################################################
		if (isTRUE(theZscoreFlag))
		{
			# z-scores
			thisData <- ((abs(thisData - mean(thisData))) / sd(thisData))
		}
		##############################################################################
		useSampleIds <- rownames(myBigDataMatrix)
		listOfData <- lapply(thisData, function(theData)
		{
			as.vector(unlist(theData))
		})
		names(listOfData) <- useSampleIds
		titlePrefix <- "log10(+1)"
		if (isTRUE(theZscoreFlag))
		{
			titlePrefix <- "z-scores"
		}
		names(listOfData) <- paste(names(listOfData), as.vector(unlist(theAllBatchesDataframe[theBatchType])), sep=" : ")
		if (0==symbol)
		{
			stripchart(x=listOfData, vertical=TRUE, pch=symbol, las=2, ylim=c(0, theMax),
								 col=colors[batchesAsLevels],
								 main=paste(titlePrefix, theBatchType, theCountType, theDisease, sep=" "))
		}
		else
		{
			stripchart(x=listOfData, vertical=TRUE, pch=symbol, las=2, ylim=c(0, theMax),
								 col=colors[batchesAsLevels],
								 main=paste(titlePrefix, theBatchType, theCountType, theDisease, sep=" "), add=TRUE)
		}
		symbol <- symbol + 1
	}
	legend(x="bottomleft", legend=basename(theMatrixFiles), pch=0:length(theMatrixFiles),
				 xpd=TRUE, inset=c(0,-.75))
	legend(x="bottomright", legend=levels(batchesAsLevels), text.col=colors[unique(batchesAsLevels)],
				 xpd=TRUE, inset=c(0,-.75))
}

readDataSets <- function(theMatrixFiles)
{
	datasets <- lapply(theMatrixFiles,
										 function(theDataFile)
	{
		dataMatrix <- readMatrixScore(theDataFile)
		fileDataVector <- rowSums(t(dataMatrix), na.rm=TRUE)
		fileDataVector
	})
	names(datasets) <- basename(theMatrixFiles)
	datasets
}

writeMutationBoxplot_Wide <- function(theOutdir, theMatrixFiles, theAllBatchesDataframe, theAllSampleIds,
																			theBatchType, theCountType, theDisease, theMax, theZscoreFlag=FALSE)
{
	message("writeMutationStripchart")
	filePost <- "Log10"
	if (isTRUE(theZscoreFlag))
	{
		filePost <- "ZScore"
	}
	outfile <- file.path(theOutdir, paste("WideBoxplot_", theBatchType, "_", theCountType, "_", theDisease, "_", filePost, "_Diagram.PNG", sep=""))
	message("write ", basename(outfile))
	##TODO: maybe write out TSV file?
	outfileTSV <- file.path(theOutdir, paste("WideBoxplot_", theBatchType, "_", theCountType, "_", theDisease, "_", filePost, ".tsv", sep=""))
	message("write ", basename(outfileTSV))
	##############################################################################
	# build list of values to plot, where each element of list is combination of a batch and a data file
	theMatrixFiles <- sort(theMatrixFiles)
	plotList <- list()
	dataSetList <- readDataSets(theMatrixFiles)
	message("writeMutationBoxplot_Wide reading data files for stripcharts")
	for (myBatch in unique(sort(as.vector(unlist(theAllBatchesDataframe[theBatchType])))))
	{
		samplesInBatch <- theAllBatchesDataframe$Sample[theAllBatchesDataframe[theBatchType]==myBatch]
		for (myIndex in 1:length(dataSetList))
		{
			myDataSet <- dataSetList[[myIndex]]
			curFilename <- names(dataSetList)[myIndex]
			curFilename <- gsub(".tsv", "", curFilename, fixed=TRUE)
			message(curFilename, " and ", myBatch)
			dataEntry <- myDataSet[samplesInBatch]
			if (!is.na(sum(dataEntry)))
			{
			  if (sum(dataEntry)>0)
			  {
			    nameOfEntry <- paste(curFilename, " - ", myBatch, " (", length(dataEntry), ")", sep="")
			    plotList[[nameOfEntry]] <- dataEntry
			  }
			}
		}
	}
	##############################################################################
	# get max value for plotting
	message("get max")
	message(theMax)
	if (isTRUE(theZscoreFlag))
	{
		# z-scores max value
		theMax <- max(as.vector(unlist(lapply(plotList, function(theColumn)
		{
			max((abs(theColumn - mean(theColumn))) / sd(theColumn))
		}))), na.rm = TRUE)
	}
	else
	{
		# log10(+1) max value
		theMax <- max(log10(theMax+1), na.rm = TRUE)
	}
	message(theMax)
	##############################################################################
	# adjust data for z-scores or log10
	message("check for z-scores or log10")
	if (isTRUE(theZscoreFlag))
	{
		plotList <- lapply(plotList, function(theBatchData) { ((abs(theBatchData - mean(theBatchData))) / sd(theBatchData)) })
	}
	else
	{
		plotList <- lapply(plotList, function(theColumn)
		{
			log10(theColumn+1)
		})
	}
	##############################################################################
	#
	symbol <- 1
	colors <- getColorListAny(length(theMatrixFiles))
	##
	message("plot")
	myWidth <- 100*length(plotList)
	if (myWidth<1000)
	{
		myWidth <- 1000
	}
	if (myWidth>5000)
	{
	  myWidth <- 5000
	}
	##
	if (!is.finite(theMax))
	{
		message("Max was not a good number, so don't try plotting")
	}
	else
	{
		message("plot")
		CairoPNG(filename=outfile, width=myWidth, height = 2000, pointsize=12)
		on.exit(dev.off(), add = TRUE)
		# c(bottom, left, top, right)
		par(mar=c(60,10,2,2))
		#par(mar=c(10,60,2,2))
		#par(mar=c(60,60,2,2))
		titlePrefix <- "log10(+1) (zeros removed)"
		if (isTRUE(theZscoreFlag))
		{
			titlePrefix <- "z-scores (zeros removed)"
		}
		boxplot(x=plotList, range=c(0, theMax), horizontal=FALSE,
						las=2,
						col=colors,
						main=paste(titlePrefix, theBatchType, theCountType, theDisease, sep=" "))
		stripchart(x=plotList, vertical=TRUE, pch=symbol, add=TRUE,
							 method="jitter", jitter=0.3, las=2, ylim=c(0, theMax),
							 col="black")
		#browser()
		####lapply(plotList, function(theCell) { print(length(as.vector(unlist(theCell))))})
		# add header line
		#paste(c("DataSet", names(plotList)), sep="\t", collapse="\t"),
		#lapply(plotList, function(theCell)
		#{
		#
		#})
		#legend(x="bottomleft", legend=basename(theMatrixFiles), pch=0:(length(theMatrixFiles)-1),
		#			 xpd=TRUE, inset=c(0,-.75), text.col=colors[1:length(theMatrixFiles)], col=colors[1:length(theMatrixFiles)])
		#legend(x="bottomright", legend=levels(batchesAsLevels), text.col=colors[unique(batchesAsLevels)],
		#			 xpd=TRUE, inset=c(0,-.75))
	}
}

writeMutationBoxplot_Narrow <- function(theOutdir, theMatrixFiles, theAllBatchesDataframe, theAllSampleIds,
																			theBatchType, theCountType, theDisease, theMax, theZscoreFlag=FALSE)
{
	message("writeMutationStripchart")
	filePost <- "Log10"
	if (isTRUE(theZscoreFlag))
	{
		filePost <- "ZScore"
	}
	theMatrixFiles <- sort(theMatrixFiles)
	for (myMatrixFile in theMatrixFiles)
	{
		subFileName <- strsplit(x=myMatrixFile, split=".", fixed=TRUE)[[1]][2]
		#subFileName <- gsub(x=subFileName, pattern=paste("_", theCountType), replacement="", fixed=TRUE)
		# only applies to GDC, but won't hurt DCC
		subFileName <- gsub(x=subFileName, pattern="AggregationandMasking", replacement="", fixed=TRUE)
		outfile <- file.path(theOutdir, paste("NarrowBoxplot_", theBatchType, "_", theCountType, "_", theDisease, "_", filePost,
		                                      "_", subFileName, "_Diagram.PNG", sep=""))
		outfileTSV <- file.path(theOutdir, paste("NarrowBoxplot_", theBatchType, "_", theCountType, "_", theDisease, "_", filePost,
		                                      "_", subFileName, ".tsv", sep=""))
		writeMutationBoxplot_Narrow_Internal(outfile, outfileTSV, myMatrixFile, theAllBatchesDataframe, theAllSampleIds,
																						theBatchType, theCountType, theDisease, theMax, theZscoreFlag)
	}
}

writeMutationBoxplot_Narrow_Internal <- function(theOutFile, theTsvFile, theMatrixFiles, theAllBatchesDataframe, theAllSampleIds,
																				theBatchType, theCountType, theDisease, theMax, theZscoreFlag)
{
  ##TODO: maybe write out TSV file?
	message("write ", basename(theOutFile))
	##############################################################################
	# build list of values to plot, where each element of list is combination of a batch and a data file
	plotList <- list()
	dataSetList <- readDataSets(theMatrixFiles)
	message("writeMutationBoxplot_Narrow_Internal reading data files for stripcharts")
	for (myBatch in unique(sort(as.vector(unlist(theAllBatchesDataframe[theBatchType])))))
	{
		samplesInBatch <- theAllBatchesDataframe$Sample[theAllBatchesDataframe[theBatchType]==myBatch]
		for (myIndex in 1:length(dataSetList))
		{
			myDataSet <- dataSetList[[myIndex]]
			curFilename <- names(dataSetList)[myIndex]
			curFilename <- gsub(".tsv", "", curFilename, fixed=TRUE)
			message(curFilename, " and ", myBatch)
			dataEntry <- myDataSet[samplesInBatch]
			dataEntry <- dataEntry[!is.na(dataEntry)]
			if (!is.na(sum(dataEntry)))
			{
			  if (sum(dataEntry)>0)
			  {
			    nameOfEntry <- paste(curFilename, " - ", myBatch, " (", length(dataEntry), ")", sep="")
			    plotList[[nameOfEntry]] <- dataEntry
			  }
			}
		}
	}
	##############################################################################
	# get average group size and remove anything below that
	aveSize <- mean(as.vector(unlist(lapply(plotList, function(theColumn)
	{
		length(theColumn)
	}))))
	newPlotList <- list()
	newNameVector <- c()
	for(myIndex in 1:length(plotList))
	{
		if (length(plotList[[myIndex]])>aveSize)
		{
			newPlotList <- c(newPlotList, list(plotList[[myIndex]]))
			newNameVector <- c(newNameVector, names(plotList)[myIndex])
		}
	}
	if (length(newPlotList)>3)
	{
		plotList <- newPlotList
		names(plotList) <- newNameVector
	}
	newPlotList <- NULL
	##############################################################################
	# get max value for plotting
	message("get max")
	message(theMax)
	if (isTRUE(theZscoreFlag))
	{
		message("Do Zscore max")
		# z-scores max value
		theMax <- max(as.vector(unlist(lapply(plotList, function(theColumn)
		{
			max((abs(theColumn - mean(theColumn))) / sd(theColumn))
		}))), na.rm = TRUE)
	}
	else
	{
		message("Do log10 max")
		# log10(+1) max value
		theMax <- max(log10(theMax+1), na.rm = TRUE)
	}
	message(theMax)
	##############################################################################
	# adjust data for z-scores or log10
	message("check for z-scores or log10")
	if (isTRUE(theZscoreFlag))
	{
		message("Do Zscore flag")
		plotList <- lapply(plotList, function(theBatchData) { ((abs(theBatchData - mean(theBatchData))) / sd(theBatchData)) })
	}
	else
	{
		message("Do log10 flag")
		plotList <- lapply(plotList, function(theColumn)
		{
			log10(theColumn+1)
		})
	}
	##############################################################################
	#
	symbol <- 1
	colors <- getColorListAny(length(theMatrixFiles))
	##
	if (!is.finite(theMax))
	{
		message("Max was not a good number, so don't try plotting")
	}
	else
	{
		message("plot")
		print(length(plotList))
		myWidth <- 100*length(plotList)
		if (myWidth<1000)
		{
		  myWidth <- 1000
		}
		if (myWidth>5000)
		{
		  myWidth <- 5000
		}
		CairoPNG(filename=theOutFile, width = myWidth, height = 2000, pointsize=12)
		on.exit(dev.off(), add = TRUE)
		# c(bottom, left, top, right)
		#TDC#par(mar=c(60,10,2,2))
		#par(mar=c(10,60,2,2))
		par(mar=c(60,2,2,2))
		titlePrefix <- "log10(+1) (zeros removed)"
		if (isTRUE(theZscoreFlag))
		{
			titlePrefix <- "z-scores (zeros removed)"
		}
		#browser()
		boxplot(x=plotList, range=c(0, theMax), horizontal=FALSE,
						las=2,
						col=colors,
						main=paste(titlePrefix, theBatchType, theCountType, theDisease, sep=" "))
		stripchart(x=plotList, vertical=TRUE, pch=symbol, add=TRUE,
							 method="jitter", jitter=0.3, las=2, ylim=c(0, theMax),
							 col="black")
		#legend(x="bottomleft", legend=basename(theMatrixFiles), pch=0:(length(theMatrixFiles)-1),
		#			 xpd=TRUE, inset=c(0,-.75), text.col=colors[1:length(theMatrixFiles)], col=colors[1:length(theMatrixFiles)])
		#legend(x="bottomright", legend=levels(batchesAsLevels), text.col=colors[unique(batchesAsLevels)],
		#			 xpd=TRUE, inset=c(0,-.75))
	}
}

writeMutationStripchart <- function(theOutdir, theMatrixFiles, theAllBatchesDataframe, theAllSampleIds,
																		theBatchType, theCountType, theDisease, theMax, theZscoreFlag=FALSE)
{
	message("writeMutationStripchart")
	outfile <- file.path(theOutdir, paste("FullMutCounts_", theBatchType, "_", theCountType, "_", theDisease, "_MutStrip_Diagram.PNG", sep=""))
	message("write ", basename(outfile))
	##############################################################################
	# build list of values to plot, where each element of list is combination of a batch and a data file
	theMatrixFiles <- sort(theMatrixFiles)
	plotList <- list()
	dataSetList <- readDataSets(theMatrixFiles)
	message("writeMutationStripchart reading data files for stripcharts")
	for (myBatch in unique(sort(as.vector(unlist(theAllBatchesDataframe[theBatchType])))))
	{
		samplesInBatch <- theAllBatchesDataframe$Sample[theAllBatchesDataframe[theBatchType]==myBatch]
		for (myIndex in 1:length(dataSetList))
		{
			myDataSet <- dataSetList[[myIndex]]
			curFilename <- names(dataSetList)[myIndex]
			curFilename <- gsub(".tsv", "", curFilename, fixed=TRUE)
			message(curFilename, " and ", myBatch)
			dataEntry <- myDataSet[samplesInBatch]
			if (!is.na(sum(dataEntry)))
			{
			  if (sum(dataEntry)>0)
			  {
			    nameOfEntry <- paste(curFilename, " - ", myBatch, " (", length(dataEntry), ")", sep="")
    			plotList[[nameOfEntry]] <- dataEntry
			  }
			}
		}
	}
	##############################################################################
	# get max value for plotting
	message("get max")
	message(theMax)
	if (isTRUE(theZscoreFlag))
	{
		# z-scores max value
		theMax <- max(as.vector(unlist(apply(plotList, 2, function(theColumn)
		{
			max((abs(theColumn - mean(theColumn))) / sd(theColumn))
		}))))
	}
	else
	{
		# log10(+1) scores
		theMax <- max(log10(theMax+1))
	}
	message(theMax)
	##############################################################################
	# adjust data for z-scores
	message("check for z-scores")
	if (isTRUE(theZscoreFlag))
	{
		plotList <- lapply(plotList, function(theBatchData) { ((abs(theBatchData - mean(theBatchData))) / sd(theBatchData)) })
	}
	##############################################################################
	#
	symbol <- 1
	colors <- getColorListAny(length(theMatrixFiles))
	##
	message("plot")
	myWidth <- length(plotList)*20
	if (myWidth<100)
	{
	  myWidth <- 100
	}
	if (myWidth>5000)
	{
	  myWidth <- 5000
	}
	CairoPNG(filename=outfile, width = myWidth, height = 2000, pointsize=12)
	on.exit(dev.off(), add = TRUE)
	# c(bottom, left, top, right)
	par(mar=c(60,10,2,2))
	#par(mar=c(10,60,2,2))
	titlePrefix <- "log10(+1) (zeros removed)"
	if (isTRUE(theZscoreFlag))
	{
		titlePrefix <- "z-scores (zeros removed)"
	}
	#TDC#stripChart?? EnvStats?
	stripchart(x=plotList, vertical=TRUE, pch=symbol,
						 method="jitter", jitter=0.3, las=2, ylim=c(0, theMax),
						 col=colors,
						 main=paste(titlePrefix, theBatchType, theCountType, theDisease, sep=" "),
						 n.text="bottom"
						 )
	#stripchart(x=plotList, vertical=TRUE, pch=symbol,
	#					 method="jitter", jitter=0.3, las=2, ylim=c(0, theMax),
	#					 col=colors,
	#					 main=paste(titlePrefix, theBatchType, theCountType, theDisease, sep=" "))
	#legend(x="bottomleft", legend=basename(theMatrixFiles), pch=0:(length(theMatrixFiles)-1),
	#			 xpd=TRUE, inset=c(0,-.75), text.col=colors[1:length(theMatrixFiles)], col=colors[1:length(theMatrixFiles)])
	#legend(x="bottomright", legend=levels(batchesAsLevels), text.col=colors[unique(batchesAsLevels)],
	#			 xpd=TRUE, inset=c(0,-.75))
}

#
# OLD_writeMutationStripchart_OLD <- function(theOutdir, theMatrixFiles, theAllBatchesDataframe, theAllSampleIds,
# 																		theBatchType, theCountType, theDisease, theMax, theZscoreFlag=FALSE)
# {
# 	message("writeMutationStripchart")
# 	outfile <- file.path(theOutdir, paste("FullMutCounts_", theBatchType, "_", theCountType, "_", theDisease, "_MutStrip_Diagram.PNG", sep=""))
# 	message("write ", basename(outfile))
# 	symbol <- 0
# 	batchesAsLevels <- as.factor(as.vector(unlist(theAllBatchesDataframe[theBatchType])))
# 	colors <- getColorListAny(length(theMatrixFiles))
# 	# build matrix of values
# 	bigDataVector <- c()
# 	for (myDataFile in theMatrixFiles)
# 	{
# 		dataMatrix <- readMatrixScore(myDataFile)
# 		fileDataVector <- rowSums(t(dataMatrix), na.rm=TRUE)
# 		plotDataVector <- rep(0, length.out=length(theAllSampleIds))
# 		names(plotDataVector) <- theAllSampleIds
# 		plotDataVector[names(fileDataVector)] <- fileDataVector
# 		bigDataVector <- c(bigDataVector, as.vector(unlist(plotDataVector)))
# 	}
# 	myBigDataMatrix <- matrix(data=as.vector(unlist(bigDataVector)),
# 														nrow=length(theAllSampleIds),
# 														ncol=length(theMatrixFiles),
# 														dimnames=list(theAllSampleIds, theMatrixFiles),
# 														byrow=TRUE)
# 	##############################################################################
# 	if (isTRUE(theZscoreFlag))
# 	{
# 		# z-scores max value
# 		theMax <- max(as.vector(unlist(apply(myBigDataMatrix, 2, function(theColumn)
# 		{
# 			max((abs(theColumn - mean(theColumn))) / sd(theColumn))
# 		}))))
# 	}
# 	else
# 	{
# 		# log10(+1) scores
# 		myBigDataMatrix <- log10(myBigDataMatrix+1)
# 		theMax <- log10(theMax+1)
# 	}
# 	##############################################################################
# 	# get groups
# 	groups <- batchesAsLevels
# 	groupNames <- as.vector(unlist(theAllBatchesDataframe[theBatchType]))
# 	#####
# 	CairoPNG(filename=outfile, width = 2000, height = 1000, pointsize=12)
# 	on.exit(dev.off(), add = TRUE)
# 	# c(bottom, left, top, right)
# 	par(mar=c(30,10,2,2))
# 	titlePrefix <- "log10(+1)"
# 	if (isTRUE(theZscoreFlag))
# 	{
# 		titlePrefix <- "z-scores"
# 	}
# 	for(dataSet in colnames(myBigDataMatrix))
# 	{
# 		thisData <- myBigDataMatrix[, dataSet]
# 		##############################################################################
# 		# z-scores
# 		if (isTRUE(theZscoreFlag))
# 		{
# 			thisData <- ((abs(thisData - mean(thisData))) / sd(thisData))
# 		}
# 		##############################################################################
# 		if (0==symbol)
# 		{
# 			stripchart(x=thisData ~ groups, vertical=TRUE, pch=symbol,
# 								 method="jitter", jitter=0.3, las=2, ylim=c(0, theMax),
# 								 col=colors[symbol+1],
# 								 main=paste(titlePrefix, theBatchType, theCountType, theDisease, sep=" "))
# 		}
# 		else
# 		{
# 			stripchart(x=thisData ~ groups, vertical=TRUE, pch=symbol,
# 								 method="jitter", jitter=0.3, las=2, ylim=c(0, theMax),
# 								 col=colors[symbol+1],
# 								 main=paste(titlePrefix, theBatchType, theCountType, theDisease, sep=" "),
# 								 add=TRUE)
# 		}
# 		symbol <- symbol + 1
# 	}
# 	legend(x="bottomleft", legend=basename(theMatrixFiles), pch=0:(length(theMatrixFiles)-1),
# 				 xpd=TRUE, inset=c(0,-.75), text.col=colors[1:length(theMatrixFiles)], col=colors[1:length(theMatrixFiles)])
# 	#legend(x="bottomright", legend=levels(batchesAsLevels), text.col=colors[unique(batchesAsLevels)],
# 	#			 xpd=TRUE, inset=c(0,-.75))
# }

buildCombinedDataframe <- function(theBatchFiles, theCountType, theAllSampleIds)
{
	# build combined dataframe for batch files
	combined <- NULL
	for (myBatchFile in theBatchFiles)
	{
		# read dataframe
		myDF <- readAsGenericDataframe(myBatchFile, theNaString="NULL")
		##TDC myDF$MBABatchFile <- myBatchFile
		if (is.null(combined))
		{
			# if first dataframe, assign value
			combined <- myDF
		}
		else
		{
			# if not first, combined using merge
			combined <- merge(combined, myDF, by=names(combined), all=TRUE)
		}
	}
	# add in unknown batch ids
	# first find sample ids not in combined batches
	notIncluded <- theAllSampleIds[!(theAllSampleIds %in% combined$Sample)]
	for(newId in notIncluded)
	{
		# assign Unknown to combined batches
		combined <- rbind(combined, c(newId, rep("Unknown", length.out=length(names(combined))-1)))
	}
	# remove ids in batch files that aren't in data files
	# by finding indexes to keep
	keepIndexes <- which(combined$Sample %in% theAllSampleIds)
	# and filtering based on that
	combined <- combined[keepIndexes,]
	return(combined)
}

findMaxValue <- function(theMatrixFiles)
{
	maxValue <- 0
	for (myDataFile in theMatrixFiles)
	{
		dataMatrix <- readMatrixScore(myDataFile)
		dataVector <- rowSums(t(dataMatrix), na.rm=TRUE)
		fileMax <- max(dataVector)
		if (fileMax>maxValue)
		{
			maxValue <- fileMax
		}
	}
	return(maxValue)
}

collectNonZeroSampleIds <- function(theMatrixFiles)
{
	sampleVector <- c()
	for (myDataFile in theMatrixFiles)
	{
		sampleLine <- readLines(myDataFile, n=1)
		fileSamples <- strsplit(x=sampleLine, split="\t", fixed=TRUE)[[1]]
		fileSamples <- fileSamples[fileSamples!=""]
		sampleVector <- unique(sort(c(sampleVector, fileSamples)))
	}
	return(sampleVector)
}

collectSampleIds <- function(theMatrixFiles)
{
	sampleVector <- c()
	for (myDataFile in theMatrixFiles)
	{
		sampleLine <- readLines(myDataFile, n=1)
		fileSamples <- strsplit(x=sampleLine, split="\t", fixed=TRUE)[[1]]
		fileSamples <- fileSamples[fileSamples!=""]
		sampleVector <- unique(sort(c(sampleVector, fileSamples)))
	}
	return(sampleVector)
}

collectLoggedBatches <- function(theResultList, theBatchType)
{
	# if the values for a batch type has batch names (in addition to the p-value)
	# collect the batch names found as significant by a Dunn's test	batches <- lapply(theResultList, function(thePvalues)
	batches <- lapply(theResultList, function(thePvalues)
	{
		if (length(thePvalues[theBatchType][[1]])>1)
		{
			# if there is more than one value for a batch type,
			# return the second value, which is a vector of batches
			return(thePvalues[theBatchType][[1]][2])
		}
		else
		{
			# if no Dunn's results, return NA
			return(NA)
		}
	})
	batches
}

getBatchAddOn <- function(theBatches)
{
	# if Dunn's test returned batch names, put them into parens to put in diagram otherwise, put an empty string
	batchAddOn <- as.vector(unlist(lapply(theBatches, function(theBatchSet)
	{
		if (is.na(theBatchSet))
		{
			# if NA, use empty string
			return("")
		}
		else
		{
			# otherwise make a "nice" display string
			paste("(",paste(as.vector(unlist(theBatchSet)), sep="", collapse=", "), ")", sep="")
		}
	})))
	batchAddOn
}

makeKruskalPng <- function(theFile, theValues, theBatches, theMain, batchAddOn, thePvalueCutoff)
{
  # PNG includes corresponding non-NA values (which may be a vector) from theBatches
  myWidth <- ((length(theValues)*25)+100)
  if (myWidth<100)
  {
    myWidth <- 100
  }
  if (myWidth>5000)
  {
    myWidth <- 5000
  }
  CairoPNG(filename=theFile, width = 1000, height = myWidth, pointsize=12)
  # not the dev.off on exit from this function, which also covers most errors
  on.exit(dev.off(), add = TRUE)
  # c(bottom, left, top, right)
  par(mar=c(5,28,2,2))
  # make a barplot
  #note: "max" value in xlim
  #note: horizontal plot
  bb <- barplot(theValues, xlab="-log10(pvalue)", ylab="", horiz=TRUE, las=2,
                main=theMain, xlim=c(0,max(c(10, theValues))))
  # add labels for Dunn's Test
  text(0, y=bb, labels=batchAddOn, pos=4, adj=0, cex=0.8)
  # add a red line at the cutoff value
  abline(v=-log10(thePvalueCutoff),  col="red")
}

makeKruskalTsv <- function(theFile, theValues, theBatches, theMain, batchAddOn, thePvalueCutoff)
{
  outDF <- data.frame(dataset=names(theValues),
                      negLog10PValue=as.vector(unlist(theValues)),
                      batchesCalled=batchAddOn,
                      negLog10Cutoff=rep(-log10(thePvalueCutoff), length(theValues)))
  writeAsDataframe(theFile, outDF)
}

dunnTestBatches <- function(dataVector, batchValues, thePvalueCutoff, theZScoreCutoff)
{
	noSpaceValues <- as.factor(gsub("[^[:alnum:]=\\.]", "", batchValues))
	names(noSpaceValues) <- names(batchValues)
	dunnResults <- dunn.test(x=dataVector, g=noSpaceValues, table=FALSE)
	dunnPairs <- dunnResults$comparisons
	dunnPVadj <- dunnResults$P.adjusted
	counts <- c()
	for(index in 1:length(dunnPairs))
	{
		batch1 <- strsplit(dunnPairs[index], " - ", fixed=TRUE)[[1]][1]
		batch2 <- strsplit(dunnPairs[index], " - ", fixed=TRUE)[[1]][2]
		myVal1 <- counts[batch1]
		myVal2 <- counts[batch2]
		if ((is.null(myVal1))||(is.na(myVal1)))
		{
			myVal1 <- 0
		}
		if ((is.null(myVal2))||(is.na(myVal2)))
		{
			myVal2 <- 0
		}
		if (dunnPVadj[index]<=thePvalueCutoff)
		{
			myVal1 <- myVal1 + 1
			myVal2 <- myVal2 + 1
		}
		counts[batch1] <- myVal1
		counts[batch2] <- myVal2
	}
	zscores <- ( (abs(counts - mean(counts))) / sd(counts))
	# the 95% Confidence Internal (1.96 Standard Deviations).
	# the 99% Confidence Internal (2.576 Standard Deviations).
	# 99.9999% (4.892 Standard Deviations)
	batchesWithEffects <- NULL
	batchesWithEffects <- names(which(zscores>theZScoreCutoff))
	if (0==length(batchesWithEffects))
	{
		batchesWithEffects <- names(which(zscores==max(zscores)))
	}
	batchesWithEffects <- unique(sort(as.vector(batchValues[which(noSpaceValues %in% batchesWithEffects)])))
	batchesWithEffects
}

loadAndRunForKruskal <- function(theMatrixFile, theBatchesFile, theBatchTypes, thePvalueCutoff, theZScoreCutoff)
{
	# theMatrixFile - full path to matrix data file
	# theBatchesFile - full path to corresponding batch file
	# theBatchTypes=c("BatchId", "PlateId", "ShipDate", "TSS") - Batch types to check (correspond to columns for batch files)
	# thePvalueCutoff=.00001 - P-value to use for cutoff for "significant" difference in Kruskal-Wallis and Dunn's tests
	# theZScoreCutoff=1.96 - cut off for Dunn's test for "significant" confidence interval
	# 			95% Confidence Internal (1.96 Standard Deviations).
	# 			99% Confidence Internal (2.576 Standard Deviations).
	# 			99.9999% (4.892 Standard Deviations)
	#
	# returns: a list of either p-values or a sub-list of p-value and batches called by Dunn's Test
	#
	# This function does Kruskal-Wallis Test by Rank for each batch type on the given data files.
	# If Kruskal-Wallis gives a p-value less than the given cutoff, a Dunn's Test is performed to find signifant batches
	message("read data")
	message(theMatrixFile)
	message(theBatchesFile)
	################################
	# read the Matrix Data File
	dataMatrix <- readMatrixScore(theMatrixFile)
	# read the Dataframe Batch File
	batchDF <- readDataframeScore(theBatchesFile)
	# print original matrix and df sizes
	print(dim(dataMatrix))
	print(dim(batchDF))
	################################
	# get sample names for which we have batch information
	sampleNames <- intersect(batchDF$Sample, colnames(dataMatrix))
	# keep only samples for which we have batch information
	dataMatrix <- dataMatrix[,sampleNames]
	# print reduced matrix and df sizes
	print(dim(dataMatrix))
	print(dim(batchDF))
	################################
	# for each batch type, do Kruskal-Wallis test by rank
	# if Kruskal-Wallis makes significance test, do Dunn's Test
	pValueVector <- lapply(theBatchTypes,
												 function(myBatchType)
												 {
												 	message(myBatchType)
												 	# do rowsums on transposed matrix to get totals for each sample across all genes
												 	dataVector <- rowSums(t(dataMatrix), na.rm=TRUE)
												 	message("length(dataVector)=")
												 	print(length(dataVector))
												 	# pull out batch information for this batch type
												 	batchValues <- batchDF[,myBatchType]
												 	# apply sample names to batch information
												 	names(batchValues) <- batchDF$Sample
												 	# get batch values for samples in this data matrix
												 	batchValues <- batchValues[names(dataVector)]
												 	# test # of batches
												 	if (length(unique(sort(batchValues)))>1)
												 	{
												 		# if there is more than one batch (can't really test one batch)
												 		#
												 		# convert batch values to factors
												 		# Kruskal-Wallis and Dunn's don't require this in docs, but do not work properly on strings
												 		batchFactors <- as.factor(batchValues)
												 		# perform Kruskal-Wallis test by rank for data and batches
												 		testResult <- kruskal.test(x=dataVector, g=batchFactors)
												 		message(myBatchType, " results are ")
												 		print(testResult)
												 		# check result values
												 		if (is.nan(testResult$p.value))
												 		{
												 			# for an NA (untestable) value, change to a 1 to say "nothing interesting here"
												 			message(myBatchType, " changing from NA to 1")
												 			testResult$p.value <- 1
												 		}
												 		else if (-log10(testResult$p.value)>=-log10(thePvalueCutoff))
												 		{
												 			# if the p-balue result is greater than or equal to the cutoff, do a Dunn's test
												 			batchesWithEffects <- dunnTestBatches(dataVector, batchValues, thePvalueCutoff, theZScoreCutoff)
												 			# if Dunn's test returned interesting batches
												 			if (length(batchesWithEffects)>0)
												 			{
												 				# return the p-value and interesting batches in a list
												 				return(list(testResult$p.value, batchesWithEffects))
												 			}
												 		}
												 		# otherwise just return the p-value
												 		return(testResult$p.value)
												 	}
												 	else
												 	{
												 		# if there is one batch, assign p-value of 1, meaning "nothing interesting here"
														message("only one batch")
												 		return(1)
												 	}
												 })
	# add the batch types to the values, for latter reference
	names(pValueVector) <- theBatchTypes
	pValueVector
}

readMatrixScore <- function(file)
{
	# read the matrix, do not modify names or convert to factors. First column is row names
	as.matrix(read.delim(file, header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE, row.names=1))
}

readDataframeScore <- function(file)
{
	# read the dataframe and strings, do not modify names or convert to factors. First column is row names
	read.delim(file, sep="\t", as.is=TRUE, check.names=FALSE, stringsAsFactors=FALSE, colClasses=c("character"))
}

collectCountTypes <- function(theTypeCountDir)
{
	# give a directory with output from MutationBatchExtractAPI.R
	# find all files that end in .tsv
	dataTypes <- list.files(theTypeCountDir, pattern=".*.tsv", full.names=TRUE, recursive=TRUE)
	# remove the .tsv from all names
	dataTypes <- gsub(x=dataTypes, pattern=".tsv", replacement="", fixed=TRUE)
	# for each shortened name, split on underscore, and return the last element
	dataTypes <- unique(sort(as.vector(unlist(lapply(strsplit(x=dataTypes, split="_", fixed=TRUE),
																									 function(theVector)
																									 {
																									 	tail(theVector, n=1)
																									 })))))
	# remove "batches" from mutation types, since batches aren't a mutation type
	dataTypes <- dataTypes[dataTypes != "batches"]
	# return values
	dataTypes
}

################################################################################
################################################################################
################################################################################
################################################################################

plotMBatch <- function(theOutputDir, theDataFile, theBatchDF, theTitle,
											 theJavaArgs, theThreads, thePCAflag, theBatchType)
{
	message("plotMBatch")
	pipelineMean <- function(x)
	{
		mean(x, na.rm=TRUE)
	}
	message("read")
	myMatrix <- readMatrixScore(theDataFile)
	print(dim(myMatrix))
	message("update batches")
	keepIndexes <- which(theBatchDF$Sample %in% colnames(myMatrix))
	theBatchDF <- theBatchDF[keepIndexes,]
	print(dim(theBatchDF))
	# remove unused batch types
	theBatchDF <- theBatchDF[,c("Sample", theBatchType)]
	print(dim(theBatchDF))
	#theDataObject <- mbatchLoadStructures(myMatrix, theBatchDF)
	runMBatchForMutCounts(myMatrix, theBatchDF, theOutputDir, theTitle,
												theJavaArgs, theThreads, thePCAflag)
}


############################################################################################

runMBatchForMutCounts <- function(theDataMatrix, theBatchesDf, theOutputDir, theTitle,
																	theJavaArgs, theThreads, thePCAflag)
{
	thePermutations <- 2000
	thePermutationThreads <- theThreads
	message("runMBatchOnStructs")
	print(dim(theDataMatrix))
	print(dim(theBatchesDf))
	sampleNames <- intersect(theBatchesDf$Sample, colnames(theDataMatrix))
	theDataMatrix <- theDataMatrix[,sampleNames]
	print(dim(theDataMatrix))
	print(dim(theBatchesDf))
	theDataObject <- mbatchLoadStructures(theGeneMatrix=theDataMatrix, theBatchDataframe=theBatchesDf, theCovariatedataframe=NULL)
	#tmp <- as.vector(unlist(apply(theDataObject@mData,1, sd, na.rm=TRUE)))
	#mySD <- 0.005
	#while((sum(tmp>=mySD)>11000)&&(mySD<1))
	#{
	#	mySD <- mySD + 0.005
	#}
	#theTitle <- paste(theTitle, " SD>=", mySD, sep="")
	theDataObject <- mbatchFilterData(theBeaData=theDataObject,
																		theBatchTypeAndValuePairsToRemove=list(),
																		theBatchTypeAndValuePairsToKeep=list(),
																		theBatchTypesToRemove=NULL,
																		theMinIqr=0, theMinSd=0, theMinMad=0)
	#theMinIqr=0, theMinSd=mySD, theMinMad=0)
	callMBatch_SupervisedClustering_Structures_MBA(theOutputDir, theDataObject, theTitle)
	####
	#### HierarchicalClustering
	####
	callMBatch_HierarchicalClustering_Structures_MBA(theOutputDir, theDataObject, theTitle)
	####
	#### PCAPlus
	####
	if (isTRUE(thePCAflag))
	{
		callMBatch_PCA_Structures_MBA(theOutputDir, theDataObject, theTitle,
															theIsPcaTrendFunction=isTrendBatch,
															theDSCPermutations=thePermutations,
															theDSCThreads=thePermutationThreads,
															theMinBatchSize=5,
															theJavaParameters=theJavaArgs,
															theSeed=314,
															theMaxGeneCount=nrow(theDataObject@mData))
	}
	####
	#### BoxPlot
	###
	pipelineMean <- function(x)
	{
		mean(x, na.rm=TRUE)
	}
	bpGene <- nrow(theDataObject@mData)
	maxBoxPlotGenes <- bpGene
	print(maxBoxPlotGenes)
	#callMBatch_BoxplotAllSamplesData_Structures(theOutputDir, theDataObject, theTitle,
	#																						theJavaParameters="-Xms12000m", theMaxGeneCount=maxBoxPlotGenes)
	#callMBatch_BoxplotAllSamplesRLE_Structures(theOutputDir, theDataObject, theTitle,
	callMBatch_BoxplotGroup_Structures_MBA(theOutputDir, theDataObject, theTitle,
																		 theJavaParameters=theJavaArgs, theMaxGeneCount=maxBoxPlotGenes,
																		 theFunction=c(pipelineMean), theFunctionName=c("Mean"))
}

################################################################################
################################################################################
################################################################################
################################################################################

callMBatch_SupervisedClustering_Structures_MBA <- function(theOutputDir, theDataObject, theTitle)
{
	# output directory
	outdir <- file.path(theOutputDir, "SupervisedClustering")
	dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
	# here, we call supervised clustering, passing a title and an output path,
	# telling it to generate a heatmap and to use the data without removing any type/values
	# and keeping any other defaults
	SupervisedClustering_Batches_Structures(theData=theDataObject,
																					theTitle=theTitle,
																					theOutputPath=outdir,
																					theDoHeatmapFlag=TRUE,
																					theBatchTypeAndValuePairsToRemove=NULL,
																					theBatchTypeAndValuePairsToKeep=NULL)
}

callMBatch_HierarchicalClustering_Structures_MBA <- function(theOutputDir, theDataObject, theTitle)
{
	# output directory
	outdir <- file.path(theOutputDir, "HierarchicalClustering")
	dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
	# here, we take all the defaults to hierarchical clustering, passing a title and an output path
	HierarchicalClustering_Structures(theData=theDataObject,
																		theTitle=theTitle,
																		theOutputPath=outdir)
}

callMBatch_PCA_Structures_MBA <- function(theOutputDir, theDataObject, theTitle,
																			theIsPcaTrendFunction,
																			theDSCPermutations=1000,
																			theDSCThreads=1,
																			theMinBatchSize=2,
																			theJavaParameters="-Xms2000m",
																			theSeed=314,
																			theMaxGeneCount=10000)
{
	# output directory
	outdir <- file.path(theOutputDir, "PCA")
	dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
	# here, we call PCA, passing a title and an output path,
	# and to use the data without removing any type/values
	# and keeping any other defaults
	PCA_Regular_Structures(theData=theDataObject,
												 theTitle=theTitle,
												 theOutputPath=outdir,
												 theBatchTypeAndValuePairsToRemove=NULL,
												 theBatchTypeAndValuePairsToKeep=NULL,
												 theDoDscPermsFileFlag = TRUE,
												 theIsPcaTrendFunction=theIsPcaTrendFunction,
												 theDSCPermutations=theDSCPermutations,
												 theDSCThreads=theDSCThreads,
												 theMinBatchSize=theMinBatchSize,
												 theJavaParameters=theJavaParameters,
												 theSeed=theSeed,
												 theMaxGeneCount=theMaxGeneCount)
}

callMBatch_BoxplotGroup_Structures_MBA <- function(theOutputDir, theDataObject, theTitle,
																							 theJavaParameters="-Xms8000m", theMaxGeneCount = 10000,
																							 theFunction=c(mean), theFunctionName=c("Mean"))
{
	# output directory
	outdir <- file.path(theOutputDir, "BoxPlot")
	dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
	# here, we call boxplot group, passing a title and an output path,
	# and to use the data without removing any type/values
	# and tu use the mean function and the label "mean"
	# and tell it to use 8G of memory
	# and keeping any other defaults
	Boxplot_Group_Structures(theData=theDataObject,
													 theTitle=theTitle,
													 theOutputPath=outdir,
													 theBatchTypeAndValuePairsToRemove=NULL,
													 theBatchTypeAndValuePairsToKeep=NULL,
													 theListOfGroupBoxFunction=theFunction,
													 theListOfGroupBoxLabels=theFunctionName,
													 theJavaParameters=theJavaParameters,
													 theMaxGeneCount=theMaxGeneCount)
}

################################################################################
################################################################################
################################################################################
################################################################################
