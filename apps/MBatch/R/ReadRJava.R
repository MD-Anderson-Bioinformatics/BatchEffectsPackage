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

################################################################################

## JavaFile loadData(String theFile, boolean theCols, boolean theRows, boolean theData)
readAsMatrix <- function(theFile, thePar="-Xmx2000m")
{
	myJavaJars <- getJarsFromDir(dirname(system.file("ReadRJava", "ReadRJava.jar", package="MBatch")))
	logDebug("readAsMatrix - thePar ", thePar)
	logDebug("readAsMatrix - theFile ", theFile)
	logDebug("readAsMatrix - Calling .jinit ", myJavaJars)
	.jinit(classpath=myJavaJars, force.init = TRUE, parameters = updateJavaParameters(thePar))
	logDebug("readAsMatrix - .jinit complete")
	logDebug("readAsMatrix before java")
	objJavaFile <- .jcall("edu/mda/bcb/readrjava/ReadRJava", returnSig="Ledu/mda/bcb/readrjava/JavaFile;",
												method="loadDoubleData",
												.jnew("java/lang/String",theFile),
												TRUE, TRUE, TRUE)
	logDebug("readAsMatrix after java")
	myData <- .jcall(objJavaFile, "[D", "getmDoubleData")
	myCols <- .jcall(objJavaFile, "[S", "getmColumns")
	myRows <- .jcall(objJavaFile, "[S", "getmRows")
	logDebug("readAsMatrix - length(myData) ", length(myData))
	logDebug("readAsMatrix - length(myCols) ", length(myCols))
	logDebug("readAsMatrix - length(myRows) ", length(myRows))
	return(matrixWithIssues(myData, nrow=length(myRows), ncol=length(myCols), byrow=TRUE,
													dimnames=list(myRows, myCols)))
}

## JavaFile loadData(String theFile)
readAsDataFrame <- function(theFile, thePar="-Xmx2000m", theUnknownString="Unknown")
{
	myJavaJars <- getJarsFromDir(dirname(system.file("ReadRJava", "ReadRJava.jar", package="MBatch")))
	logDebug("readAsDataFrame - thePar ", thePar)
	logDebug("readAsDataFrame - theFile ", theFile)
	logDebug("readAsDataFrame - Calling .jinit ", myJavaJars)
	.jinit(classpath=myJavaJars, force.init = TRUE, parameters = updateJavaParameters(thePar))
	logDebug("readAsDataFrame - .jinit complete")
	logDebug("readAsDataFrame before java")
	objJavaFile <- .jcall("edu/mda/bcb/readrjava/ReadRJava", returnSig="Ledu/mda/bcb/readrjava/JavaFile;",
												method="loadStringData",
												.jnew("java/lang/String",theFile))
	logDebug("readAsDataFrame after java")
	myData <- .jcall(objJavaFile, "[S", "getmStringData")
	myCols <- .jcall(objJavaFile, "[S", "getmColumns")
	myRows <- .jcall(objJavaFile, "[S", "getmRows")
	# drop first column as being empty since it is the row label column
	logDebug("readAsDataFrame - length(myData) ", length(myData))
	logDebug("readAsDataFrame - length(myCols) ", length(myCols))
	logDebug("readAsDataFrame - length(myRows) ", length(myRows))
	logDebug("readAsDataFrame - myCols ", paste(myCols, collapse=", "))
	logDebug("readAsDataFrame - myRows ", paste(myRows, collapse=", "))
	df <- data.frame(matrixWithIssues(myData, ncol=length(myCols), byrow=TRUE,
	                                  dimnames=list(myRows, myCols)), stringsAsFactors=FALSE)
	colnames(df) <- myCols
	df[df==""] <- theUnknownString
	return(df)
}

## boolean writeDoubleData(String theFile, String [] theCols, String [] theRows, double [] theData)
writeAsMatrix <- function(theFile, theMatrix, thePar="-Xmx2000m")
{
  myJavaJars <- getJarsFromDir(dirname(system.file("ReadRJava", "ReadRJava.jar", package="MBatch")))
	logDebug("writeAsMatrix - thePar ", thePar)
	logDebug("writeAsMatrix - theFile ", theFile)
	myCols <- as.vector(colnames(theMatrix))
	myRows <- rownames(theMatrix)
	myData <- as.numeric(as.vector(t(theMatrix)))
	if (!is.null(myRows))
	{
		myRows <- as.vector(myRows)
	}
	logDebug("writeAsMatrix - length(myData) ", length(myData))
	logDebug("writeAsMatrix - length(myCols) ", length(myCols))
	logDebug("writeAsMatrix - length(myRows) ", length(myRows))
	logDebug("writeAsMatrix - Calling .jinit ", myJavaJars)
	.jinit(classpath=myJavaJars, force.init = TRUE, parameters = updateJavaParameters(thePar))
	logDebug("writeAsMatrix - .jinit complete")
	logDebug("writeAsMatrix before java")
	if (is.null(myRows))
	{
		success <- .jcall("edu/mda/bcb/readrjava/ReadRJava", returnSig="Z",
											method="writeDoubleData_Column",
											.jnew("java/lang/String",theFile),
											.jarray(myCols),
											.jarray(myData))
	}
	else
	{
		success <- .jcall("edu/mda/bcb/readrjava/ReadRJava", returnSig="Z",
											method="writeDoubleData_All",
											.jnew("java/lang/String",theFile),
											.jarray(myCols),
											.jarray(myRows),
											.jarray(myData))
	}
	logDebug("writeAsMatrix after java")
	logDebug("writeAsMatrix success=", success)
	return(success)
}

## boolean writeStringData(String theFile, String [] theCols, String [] theRows, String [] theData)
writeAsDataframe <- function(theFile, theDataframe, thePar="-Xmx2000m", theIncludeRowNamesFlag=FALSE)
{
  myJavaJars <- getJarsFromDir(dirname(system.file("ReadRJava", "ReadRJava.jar", package="MBatch")))
	logDebug("writeAsDataframe - thePar ", thePar)
	logDebug("writeAsDataframe - theFile ", theFile)
	myCols <- as.vector(colnames(theDataframe))
	myRows <- NULL
	myData <- as.vector(t(theDataframe))
	if(TRUE==theIncludeRowNamesFlag)
	{
		myRows <- as.vector(rownames(theDataframe))
	}
	logDebug("writeAsDataframe - length(myData) ", length(myData))
	logDebug("writeAsDataframe - length(myCols) ", length(myCols))
	logDebug("writeAsDataframe - length(myRows) ", length(myRows))
	logDebug("writeAsDataframe - Calling .jinit ", myJavaJars)
	.jinit(classpath=myJavaJars, force.init = TRUE, parameters = updateJavaParameters(thePar))
	logDebug("writeAsDataframe - .jinit complete")
	logDebug("writeAsDataframe before java")
	if (is.null(myRows))
	{
		success <- .jcall("edu/mda/bcb/readrjava/ReadRJava", returnSig="Z",
										method="writeStringData_Column",
										.jnew("java/lang/String",theFile),
										.jcastToArray(myCols),
										.jcastToArray(myData))
	}
	else
	{
		success <- .jcall("edu/mda/bcb/readrjava/ReadRJava", returnSig="Z",
											method="writeStringData_All",
											.jnew("java/lang/String",theFile),
											.jcastToArray(myCols),
											.jcastToArray(myRows),
											.jcastToArray(myData))
	}
	logDebug("writeAsDataframe after java")
	logDebug("writeAsDataframe success=", success)
	return(success)
}

