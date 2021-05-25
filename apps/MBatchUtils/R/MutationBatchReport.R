# MBatchUtils Copyright (c) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021 University of Texas MD Anderson Cancer Center
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
# MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>

########################################
########################################
########################################
########################################

htmlMutationBatchEffects <- function(theBaseDir)
{
  message(mbatchUtilVersion())
	message("htmlMutationBatchEffects")
	baseName <- basename(theBaseDir)
	refDirs <- list.dirs(theBaseDir, full.names=TRUE, recursive=FALSE)
	refFiles <- c()
	for (myRefDir in refDirs)
	{
		refFiles <- c(refFiles, list.files(myRefDir, pattern="callReference.tsv", full.names=TRUE, recursive=FALSE))
	}
	for (myRef in refFiles)
	{
		message(myRef)
		outDir <- dirname(myRef)
		outFile <- file.path(outDir, "index.html")
		if (file.exists(outFile))
		{
			file.remove(outFile)
		}
		myDF <- readAsGenericDataframe(myRef)
		myDF$Disease <- sapply(myDF$MutationFile, function(theEntry)
		{
			strsplit(theEntry, ".", fixed=TRUE)[[1]][1]
		})
		myDF$Caller <- sapply(myDF$MutationFile, function(theEntry)
		{
			paste(tail(strsplit(theEntry, ".", fixed=TRUE)[[1]], -1), collapse=".")
		})
		myDF$Source <- baseName
		myDF$Batches <- lapply(myDF$Batches, function(theEntry)
		{
			strsplit(gsub(pattern=")", replacement="", x=gsub(pattern="(", replacement="", x=theEntry, fixed=TRUE), fixed=TRUE), ", ", fixed=TRUE)[[1]]
		})
		cat("<html>", "\n", file=outFile, append=TRUE, sep="")
		cat("\t", "<head>", "\n", file=outFile, append=TRUE, sep="")
		cat("\t\t", "<style>", "\n", file=outFile, append=TRUE, sep="")
		#cat("\t\t", "img { transition: -webkit-transform 0.25s ease; transition: transform 0.25s ease; }", "\n", file=outFile, append=TRUE, sep="")
		#cat("\t\t", "img:active { -webkit-transform: scale(5); transform: scale(5); }", "\n", file=outFile, append=TRUE, sep="")
		cat("\t\t", ".hideme { display: none; }", "\n", file=outFile, append=TRUE, sep="")
		cat("\t\t", ".box { height: 500px; width: 500px; object-fit: contain; }", "\n", file=outFile, append=TRUE, sep="")
		cat("\t\t", ".box.big { height: auto; width: auto; }", "\n", file=outFile, append=TRUE, sep="")
		#cat("\t\t", ".rotbox { height: 500px; width: 500px; transform: rotate(90deg); object-fit: contain; }", "\n", file=outFile, append=TRUE, sep="")
		#cat("\t\t", ".rotbox.big { height: auto; width: auto; transform: rotate(90deg); }", "\n", file=outFile, append=TRUE, sep="")
		cat("\t\t", ".rotbox { height: 500px; width: 500px; object-fit: contain; }", "\n", file=outFile, append=TRUE, sep="")
		cat("\t\t", ".rotbox.big { height: auto; width: auto; }", "\n", file=outFile, append=TRUE, sep="")
		cat("\t\t", "</style>", "\n", file=outFile, append=TRUE, sep="")
		cat("\t", "</head>", "\n", file=outFile, append=TRUE, sep="")
		cat("\t", "<body>", "\n", file=outFile, append=TRUE, sep="")
		# iterate over rows of DF
		for (index in 1:nrow(myDF))
		{
			row = myDF[index, ];
			possibleFiles <- c(
				# MutBatch Output
				# MutDots image:		FullMutCounts_<BatchType>_<Mutation_Type>_<Disease>_MutDots_Diagram.PNG
				file.path(outDir, paste("FullMutCounts", row$BatchType, row$MutationType, row$Disease, "MutDots_Diagram.PNG", sep="_")),
				# Narrow Boxplots:	NarrowBoxplot_<BatchType>_<Mutation_Type>_<Disease>_Log10_<Caller>_<Mutation_Type>_Diagram.PNG
				file.path(outDir, paste( paste("NarrowBoxplot", row$BatchType, row$MutationType, row$Disease, "Log10",
																			 strsplit(gsub(pattern="AggregationandMasking", x=row$Caller, replacement="", fixed=TRUE), ".", fixed=TRUE)[[1]][1],
																			 row$MutationType, sep="_"),
																 "_Diagram.PNG", sep="")),
				# Narrow Boxplots:	NarrowBoxplot_<BatchType>_<Mutation_Type>_<Disease>_Log10_<Caller>_Diagram.PNG
				file.path(outDir, paste( paste("NarrowBoxplot", row$BatchType, row$MutationType, row$Disease, "Log10",
																			 strsplit(gsub(pattern="AggregationandMasking", x=row$Caller, replacement="", fixed=TRUE), ".", fixed=TRUE)[[1]][1],
																			 sep="_"),
																 "_Diagram.PNG", sep="")),
				# Narrow Boxplots:	NarrowBoxplot_<BatchType>_<Mutation_Type>_<Disease>_ZScore_<Caller>_<Mutation_Type>_Diagram.PNG
				file.path(outDir, paste( paste("NarrowBoxplot", row$BatchType, row$MutationType, row$Disease, "ZScore",
																			 strsplit(gsub(pattern="AggregationandMasking", x=row$Caller, replacement="", fixed=TRUE), ".", fixed=TRUE)[[1]][1],
																			 row$MutationType, sep="_"),
																 "_Diagram.PNG", sep="")),
				# Narrow Boxplots:	NarrowBoxplot_<BatchType>_<Mutation_Type>_<Disease>_ZScore_<Caller>_Diagram.PNG
				file.path(outDir, paste( paste("NarrowBoxplot", row$BatchType, row$MutationType, row$Disease, "ZScore",
																			 strsplit(gsub(pattern="AggregationandMasking", x=row$Caller, replacement="", fixed=TRUE), ".", fixed=TRUE)[[1]][1],
																			 sep="_"),
																 "_Diagram.PNG", sep="")),
				# Wide Boxplot:			WideBoxplot_<BatchType>_<Mutation_Type>_<Disease>_Log10_Diagram.PNG
				file.path(outDir, paste( paste("WideBoxplot", row$BatchType, row$MutationType, row$Disease, "Log10", sep="_"),
																 "_Diagram.PNG", sep="")),
				# Wide Boxplot:			WideBoxplot_<BatchType>_<Mutation_Type>_<Disease>_ZScore_Diagram.PNG
				file.path(outDir, paste( paste("WideBoxplot", row$BatchType, row$MutationType, row$Disease, "ZScore", sep="_"),
																 "_Diagram.PNG", sep="")),
				# MBatch Output
				# BoxPlot Diagram:	<MutationFile>_<MutationType>/BoxPlot/Group-MEAN/BoxPlot_Group-MEAN_Diagram-<BatchType>_Diagram.PNG
				file.path(outDir, paste(row$MutationFile, row$MutationType, sep="_"), "BoxPlot", "Group-MEAN",
									paste( "BoxPlot_Group-MEAN_Diagram-", row$BatchType, "_Diagram.PNG", sep="")),
				# BoxPlot Legend:		<MutationFile>_<MutationType>/BoxPlot/Group-MEAN/BoxPlot_Group-MEAN_Legend-<BatchType>_Diagram.PNG
				file.path(outDir, paste(row$MutationFile, row$MutationType, sep="_"), "BoxPlot", "Group-MEAN",
									paste( "BoxPlot_Group-MEAN_Legend-", row$BatchType, "_Diagram.PNG", sep="")),
				# HierClust	Diagram:	<MutationFile>_<MutationType>/HierarchicalClustering/HierarchicalClustering_Diagram_Diagram.PNG
				file.path(outDir, paste(row$MutationFile, row$MutationType, sep="_"), "HierarchicalClustering",
									"HierarchicalClustering_Diagram_Diagram.PNG"),
				# HierClust	Legend:		<MutationFile>_<MutationType>/HierarchicalClustering/HierarchicalClustering_Legend-<BatchType>_Diagram.PNG
				file.path(outDir, paste(row$MutationFile, row$MutationType, sep="_"), "HierarchicalClustering",
									paste( "HierarchicalClustering_Legend-", row$BatchType, "_Diagram.PNG", sep="")),
				# SuperClust	Diagram:	<MutationFile>_<MutationType>/SupervisedClustering/Batches/SupervisedClust_Diagram-<BatchType>_Diagram.PNG
				file.path(outDir, paste(row$MutationFile, row$MutationType, sep="_"), "SupervisedClustering", "Batches",
									paste( "SupervisedClust_Diagram-", row$BatchType, "_Diagram.PNG", sep="")),
				# SuperClust	Legend:		<MutationFile>_<MutationType>/SupervisedClustering/Batches/SupervisedClust_Legend-<BatchType>_Diagram.PNG
				file.path(outDir, paste(row$MutationFile, row$MutationType, sep="_"), "SupervisedClustering", "Batches",
									paste( "SupervisedClust_Legend-", row$BatchType, "_Diagram.PNG", sep="")),
				# PCAPlus Diagram:	<MutationFile>_<MutationType>/PCA/<BatchType>/ManyToMany/PCA-Plus/ALL_Comp1_Comp2_Diagram_Diagram.PNG
				file.path(outDir, paste(row$MutationFile, row$MutationType, sep="_"), "ManyToMany", "PCA-Plus", "ALL_Comp1_Comp2_Diagram_Diagram.PNG"),
				# PCAPlus Legend:		<MutationFile>_<MutationType>/PCA/<BatchType>/ManyToMany/PCA-Plus/ALL_Comp1_Comp2_Legend-ALL_Diagram.PNG
				file.path(outDir, paste(row$MutationFile, row$MutationType, sep="_"), "ManyToMany", "PCA-Plus", "ALL_Comp1_Comp2_Legend-ALL_Diagram.PNG")
			)
			cat("\t\t", "<br><hr>", "\n", file=outFile, append=TRUE, sep="")
			# table, with data on left, and buttons to display images at the right
			cat("\t\t", "<table>", "\n", file=outFile, append=TRUE, sep="")
			cat("\t\t\t", "<tr>", "\n", file=outFile, append=TRUE, sep="")
			cat("\t\t\t\t", "<td>", "\n", file=outFile, append=TRUE, sep="")
			cat("\t\t\t\t\t", "<strong>", row$BatchType, " ", row$MutationType, "</strong><br>", "\n", file=outFile, append=TRUE, sep="")
			cat("\t\t\t\t\t", "<strong>", row$Disease, " ", row$Caller, "</strong><br>", "\n", file=outFile, append=TRUE, sep="")
			cat("\t\t\t\t\t", "<strong>", row$MutationFile, "</strong><br>", "\n", file=outFile, append=TRUE, sep="")
			cat("\t\t\t\t\t", "<strong>", paste(row$Batches[[1]], collapse=", ", sep=""), "</strong><br>", "\n", file=outFile, append=TRUE, sep="")
			cat("\t\t\t\t\t", "<textarea rows=\"4\" cols=\"40\"></textarea><br>", "\n", file=outFile, append=TRUE, sep="")
			cat("\t\t\t\t", "</td>", "\n", file=outFile, append=TRUE, sep="")
			cat("\t\t\t\t", "<td>", "\n", file=outFile, append=TRUE, sep="")
			for (myFile in possibleFiles)
			{
				if (file.exists(myFile))
				{
					cat("\t\t\t\t\t", "<button onclick=\"{ var ele= document.getElementById(\'",
							paste(index,
										gsub(pattern="/", x=gsub(pattern=dirname(outDir), x=myFile, replacement="", fixed=TRUE), replacement="", fixed=TRUE),
										sep="_"),
							"\').classList.toggle(\'hideme\'); return false; }\">",
							basename(myFile), "</button><br>",
							"\n", file=outFile, append=TRUE, sep="")
				}
			}
			cat("\t\t\t\t", "</td>", "\n", file=outFile, append=TRUE, sep="")
			cat("\t\t\t", "</tr>", "\n", file=outFile, append=TRUE, sep="")
			cat("\t\t", "</table>", "\n", file=outFile, append=TRUE, sep="")
			for (myFile in possibleFiles)
			{
				if (file.exists(myFile))
				{
					# id to img entries to allow displaying on and off
					if ((grepl(pattern =  "NarrowBoxplot", x = basename(myFile), fixed = TRUE))&&(!(grepl(pattern =  "Hierarchical", x = basename(myFile), fixed = TRUE))))
					{
						cat("\t\t\t", "<img id=\"", paste(index,
																							 gsub(pattern="/", x=gsub(pattern=dirname(outDir), x=myFile, replacement="", fixed=TRUE), replacement="", fixed=TRUE),
																							 sep="_"),
								"\" class=\"rotbox\" onclick=\"this.classList.toggle('big'); return false;\" alt=\"", basename(myFile), "\" style=\"object-fit: contain;\" src=\"", gsub(paste(outDir, "/", sep=""), myFile, replacement="", fixed=TRUE), "\" >", "\n",
								file=outFile, append=TRUE, sep="")
					}
					else
					{
						cat("\t\t\t", "<img id=\"", paste(index,
																							gsub(pattern="/", x=gsub(pattern=dirname(outDir), x=myFile, replacement="", fixed=TRUE), replacement="", fixed=TRUE),
																							sep="_"),
								"\" class=\"box\" onclick=\"this.classList.toggle('big'); return false;\" alt=\"", basename(myFile), "\" style=\"object-fit: contain;\" src=\"", gsub(paste(outDir, "/", sep=""), myFile, replacement="", fixed=TRUE), "\" >", "\n", file=outFile, append=TRUE, sep="")
					}
				}
			}
		}
		cat("\t", "</body>", "\n", file=outFile, append=TRUE, sep="")
		cat("</html>", "\n", file=outFile, append=TRUE, sep="")
	}
	NULL
}

########################################
########################################
########################################
########################################

collectMutationBatchEffects <- function(theBaseDirs)
{
	baseNames <- basename(theBaseDirs)
	# a data.frame with (from callReference.tsv) MutationType, BatchType, MutationFile, Batches,
	# and (calculated) Disease, Source, Caller
	# Batches will be converted to a vector, but since it is inside a dataframe, this will be a list with a vector in it
	# Source is basename(baseDir)
	# Disease is taken from first half of MutationFile
	# Caller is taken from second half MutationFile
	# TCGA-COAD.MuTect2VariantAggregationandMasking_HG38
	# brca.illuminaga.genome_wustl_edu.Level_2_HG19
	#
	effectsDF <- NULL
	for(baseDir in theBaseDirs)
	{
		refDirs <- list.dirs(baseDir, full.names=TRUE, recursive=FALSE)
		refFiles <- c()
		for (myRefDir in refDirs)
		{
			refFiles <- c(refFiles, list.files(myRefDir, pattern="callReference.tsv", full.names=TRUE, recursive=FALSE))
		}
		for (myRef in refFiles)
		{
			message(myRef)
			myDF <- readAsGenericDataframe(myRef)
			myDF$Disease <- sapply(myDF$MutationFile, function(theEntry)
			{
				toupper(gsub(pattern="TCGA-", replacement="", x=strsplit(theEntry, ".", fixed=TRUE)[[1]][1], fixed=TRUE))
			})
			myDF$Caller <- sapply(myDF$MutationFile, function(theEntry)
			{
				paste(tail(strsplit(theEntry, ".", fixed=TRUE)[[1]], -1), collapse=".")
			})
			myDF$Source <- basename(baseDir)
			myDF$Batches <- lapply(myDF$Batches, function(theEntry)
			{
				strsplit(gsub(pattern=")", replacement="", x=gsub(pattern="(", replacement="", x=theEntry, fixed=TRUE), fixed=TRUE), ", ", fixed=TRUE)[[1]]
			})
			if (is.null(effectsDF))
			{
				effectsDF <- myDF
			}
			else
			{
				#effectsDF <- merge(effectsDF, myDF, by=c("MutationType", "BatchType", "MutationFile", "Batches", "Disease", "Source"))
				effectsDF <- rbind(effectsDF, myDF)
			}
		}
	}
	# collate by MutationType, BatchType, Disease
	# for each disease
	collatedDF <- NULL
	# Disease, BatchType, MutationType, MatchingBatches, UnmatchedBatches
	for (myDisease in unique(sort(effectsDF$Disease)))
	{
		disDF <- effectsDF[effectsDF$Disease==myDisease,]
		for (myBatchType in unique(sort(disDF$BatchType)))
		{
			mutTypDF <- disDF[disDF$BatchType==myBatchType,]
			for (myMutType in unique(sort(mutTypDF$MutationType)))
			{
				message(myDisease, " ", myBatchType, " ", myMutType)
				myData <- mutTypDF[mutTypDF$MutationType==myMutType,]
				myDF <- data.frame(Disease=myDisease,
													 BatchType=myBatchType,
													 MutationType=myMutType,
													 SingleBatches=getSingleBatches(myData),
													 MatchingBatches=getMatchingBatches(myData),
													 UnmatchedBatches=getUnmatchedBatches(myData))
				if (is.null(myDF))
				{
					collatedDF <- myDF
				}
				else
				{
					collatedDF <- rbind(collatedDF, myDF)
				}
			}
		}
	}
	collatedDF
}

getSingleBatches <- function(theData)
{
	result <- ""
	if (length(theData$Batches)==1)
	{
		intersection <- theData$Batches[[1]]
		result <- paste(paste(intersection, collapse=", "),
										paste(unique(sort(theData$Caller)), collapse=", "),
										sep=" : ")
	}
	result
}

getMatchingBatches <- function(theData)
{
	result <- ""
	intersection <- c()
	if (length(theData$Batches)>1)
	{
		intersection <- theData$Batches[[1]]
		for (index in 2:length(theData$Batches))
		{
			intersection <- intersect(intersection, theData$Batches[[index]])
		}
	}
	if (length(intersection)>0)
	{
		result <- paste(paste(intersection, collapse=", "),
										paste(unique(sort(theData$Caller)), collapse=", "),
										sep=" : ")
	}
	result
}

getUnmatchedBatches <- function(theData)
{

	intersection <- theData$Batches[[1]]
	if (length(theData$Batches)>1)
	{
		for (index in 2:length(theData$Batches))
		{
			intersection <- intersect(intersection, theData$Batches[[index]])
		}
	}
	########
	unmatched <- lapply(1:length(theData$Batches), function(theIndex, theBatches, theCallers)
	{
		result <- c()
		if (length(setdiff(theBatches[[theIndex]], intersection))>0)
		{
			result <- paste(paste(setdiff(theBatches[[theIndex]], intersection), collapse=", "),
											theCallers[theIndex],
											sep=" : ")
		}
		result
	}, theData$Batches, theData$Caller)
	paste(as.vector(unlist(unmatched)), collapse=" \n")
}



########################################
########################################
########################################
########################################
