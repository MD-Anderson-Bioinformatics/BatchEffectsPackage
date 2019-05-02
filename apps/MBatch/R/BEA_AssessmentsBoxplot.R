#MBatch Copyright ? 2011, 2012, 2013, 2014, 2015, 2016, 2017 University of Texas MD Anderson Cancer Center
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

boxplotJinit <- function(theJavaParameters=c("-Xms8000m", "-Djava.awt.headless=true"))
{
	myClass1 <- system.file("BoxplotJava", "jcommon-1.0.17.jar", package="MBatch")
	myClass2 <- system.file("BoxplotJava", "jfreechart-1.0.14.jar", package="MBatch")
	myClass3 <- system.file("BoxplotJava", "commons-lang3-3.3.2.jar", package="MBatch")
	myClass4 <- system.file("BoxplotJava", "commons-math3-3.3.jar", package="MBatch")
	myClass5 <- system.file("BoxplotJava", "BoxplotJava.jar", package="MBatch")
	myClass6 <- system.file("BoxplotJava", "LegendJava.jar", package="MBatch")
	myJavaJars <- file.path(myClass1, myClass2, myClass3, myClass4, myClass5, myClass6, fsep=.Platform$path.sep)
	#logDebug("boxplotJinit - Calling .jinit ", myJavaJars)
	.jinit(classpath=myJavaJars, force.init = TRUE, parameters=theJavaParameters)
}

###########################################################################################
###########################################################################################

createBatchEffectsOutput_BoxPlot_AllSampleRLE<-function(theMatrixGeneData, theDataframeBatchData,
																												theTitle, theOutputPath, thePngFlag,
																												theJavaParameters=c("-Xms8000m", "-Djava.awt.headless=true"))
{
	title <- paste(theTitle, "AllSample-RLE", sep=" ")
	checkCreateDir(theOutputPath)
	boxplotJinit(theJavaParameters)
	logDebug("dim(theMatrixGeneData)", dim(theMatrixGeneData))
	logDebug("length(colnames(theMatrixGeneData))", length(colnames(theMatrixGeneData)))
	logDebug("length(rownames(theMatrixGeneData))", length(rownames(theMatrixGeneData)))
	logDebug("dim(theDataframeBatchData)", dim(theDataframeBatchData))
	logDebug("length(names(theDataframeBatchData))", length(names(theDataframeBatchData)))
	success <- .jcall("edu/mda/bioinfo/boxplotjava/BoxplotJava",
										returnSig = "Z",
										method="allSampleRLE",
										.jarray(as.vector(as.numeric(theMatrixGeneData)), contents.class="[D"),
										.jarray(as.vector(as.character(colnames(theMatrixGeneData)))),
										.jarray(as.vector(as.character(rownames(theMatrixGeneData)))),
										.jarray(as.vector(as.character(unlist(theDataframeBatchData[2:length(names(theDataframeBatchData))])))),
										.jarray(as.vector(as.character(names(theDataframeBatchData)[2:length(names(theDataframeBatchData))]))),
										.jnew("java/lang/String",as.character(theOutputPath)),
										.jnew("java/lang/String",as.character(title)))
	if (FALSE==success)
	{
		logError("createBatchEffectsOutput_BoxPlot_AllSampleRLE Java indicated error")
	}
	logDebug("after allSampleRLE call")
	success
}

createBatchEffectsOutput_BoxPlot_AllSampleData<-function(theMatrixGeneData, theDataframeBatchData,
																							 theTitle, theOutputPath, thePngFlag,
																							 theJavaParameters=c("-Xms8000m", "-Djava.awt.headless=true"))
{
	title <- paste(theTitle, "AllSample-Data", sep=" ")
	checkCreateDir(theOutputPath)
	boxplotJinit(theJavaParameters)
# 	logDebug("dim(theMatrixGeneData)", dim(theMatrixGeneData))
# 	logDebug("length(colnames(theMatrixGeneData))", length(colnames(theMatrixGeneData)))
# 	logDebug("length(rownames(theMatrixGeneData))", length(rownames(theMatrixGeneData)))
# 	logDebug("dim(theDataframeBatchData)", dim(theDataframeBatchData))
# 	logDebug("length(names(theDataframeBatchData))", length(names(theDataframeBatchData)))
# 	logDebug("theMatrixGeneData[1:10,1]", paste(theMatrixGeneData[1:10,1], sep=", ", collapse="; "))
# 	logDebug("theMatrixGeneData[1:10,2]", paste(theMatrixGeneData[1:10,2], sep=", ", collapse="; "))
# 	logDebug("theDataframeBatchData$Sample[1:5]", paste(theDataframeBatchData$Sample[1:5], sep=", ", collapse="; "))
# 	logDebug("theDataframeBatchData$BatchId[1:5]", paste(theDataframeBatchData$BatchId[1:5], sep=", ", collapse="; "))
# 	logDebug("theDataframeBatchData$PlateId[1:5]", paste(theDataframeBatchData$PlateId[1:5], sep=", ", collapse="; "))
# 	logDebug("theDataframeBatchData$ShipDate[1:5]", paste(theDataframeBatchData$ShipDate[1:5], sep=", ", collapse="; "))
# 	logDebug("theDataframeBatchData$TSS[1:5]", paste(theDataframeBatchData$TSS[1:5], sep=", ", collapse="; "))
#

	success <- .jcall("edu/mda/bioinfo/boxplotjava/BoxplotJava",
										returnSig = "Z",
										method="allSampleData",
										.jarray(as.vector(as.numeric(theMatrixGeneData)), contents.class="[D"),
										.jarray(as.vector(as.character(colnames(theMatrixGeneData)))),
										.jarray(as.vector(as.character(rownames(theMatrixGeneData)))),
										.jarray(as.vector(as.character(unlist(theDataframeBatchData[2:length(names(theDataframeBatchData))])))),
										.jarray(as.vector(as.character(names(theDataframeBatchData)[2:length(names(theDataframeBatchData))]))),
										.jnew("java/lang/String",as.character(theOutputPath)),
										.jnew("java/lang/String",as.character(title)))
	if (FALSE==success)
	{
		logError("createBatchEffectsOutput_BoxPlot_AllSampleData Java indicated error")
	}
	logDebug("after allSampleData call")
	success
}

createBatchEffectsOutput_BoxPlot_Group<-function(theMatrixGeneData, theDataframeBatchData,
																								 theTitle, theOutputPath,
																								 theListOfGroupBoxFunction, theListOfGroupBoxLabels,
																								 thePngFlag, theJavaParameters=c("-Xms8000m", "-Djava.awt.headless=true"))
{
	title <- paste(theTitle, "Group", sep=" ")
	stopifnotWithLogging("Group box list and group box labels should be the same length", length(theListOfGroupBoxFunction)==length(theListOfGroupBoxLabels))
	checkCreateDir(theOutputPath)
	boxplotJinit(theJavaParameters)
	logDebug("dim(theMatrixGeneData)", dim(theMatrixGeneData))
	logDebug("length(colnames(theMatrixGeneData))", length(colnames(theMatrixGeneData)))
	logDebug("length(rownames(theMatrixGeneData))", length(rownames(theMatrixGeneData)))
	logDebug("dim(theDataframeBatchData)", dim(theDataframeBatchData))
	logDebug("length(names(theDataframeBatchData))", length(names(theDataframeBatchData)))
	success <- .jcall("edu/mda/bioinfo/boxplotjava/BoxplotJava",
										returnSig = "Z",
										method="groupFunction",
										.jarray(as.vector(as.numeric(theMatrixGeneData)), contents.class="[D"),
										.jarray(as.vector(as.character(colnames(theMatrixGeneData)))),
										.jarray(as.vector(as.character(rownames(theMatrixGeneData)))),
										.jarray(as.vector(as.character(unlist(theDataframeBatchData[2:length(names(theDataframeBatchData))])))),
										.jarray(as.vector(as.character(names(theDataframeBatchData)[2:length(names(theDataframeBatchData))]))),
										.jnew("java/lang/String",as.character(theOutputPath)),
										as.integer(2),
										.jnew("java/lang/String",as.character(title)))
	if (FALSE==success)
	{
		logError("createBatchEffectsOutput_BoxPlot_Group Java indicated error")
	}
	logDebug("after groupFunction call")
	success
}

###########################################################################################
###########################################################################################
