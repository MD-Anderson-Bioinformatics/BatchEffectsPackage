#MBatch Copyright ? 2011, 2012, 2013, 2014, 2015, 2016, 2017 University of Texas MD Anderson Cancer Center
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.


################################################################################

createDirPlusFilename<-function(theDir, ...)
{
	file.path(theDir, paste(..., sep="", collapse=""))
}

createDirPath<-function(...)
{
	file.path(...)
}

checkDirForCreation <- function(thePath)
{
	if (thePath!=dirname(thePath))
	{
		checkDirForCreation(dirname(thePath))
		if (FALSE==file.exists(thePath))
		{
			dir.create(thePath, recursive=FALSE)
		}
	}
}

checkCreateDir<-function(...)
{
	myDir <-createDirPath(...)
	logDebug("checkCreateDir: ", myDir)
	checkDirForCreation(myDir)
	return(myDir)
}

################################################################################

readAsGenericMatrix_Samples<-function(theFile)
{
	myFile <- file(theFile, "r")
	on.exit(close(myFile))
	myHeaderString <- readLines(con=myFile, n=1)
	myHeaderList <- unlist(strsplit(myHeaderString, "\t"))
	myHeaderList
}

readAsGenericMatrix <- function(theFile)
{
	samples <- readAsGenericMatrix_Samples(theFile)
	samples <- samples[2:length(samples)]
	whatList <- lapply(0:length(samples), function(x)
	{
		if (0==x)
		{
			return(character())
		}
		else
		{
			return(double())
		}
	})
	myScan <- scan(file=theFile, skip=1, what=whatList, quote="",
								 sep="\t", na.strings="", flush=TRUE, fill=FALSE, multi.line=FALSE)
	genes <- as.vector(unlist(myScan[1]))
	data <- as.vector(unlist(myScan[2:length(myScan)]))
	temp<-matrixWithIssues(data,
												 ncol=length(samples),
												 nrow=length(genes),
												 dimnames=list(make.unique(genes), make.unique(samples)))
	temp <- temp[,sort(colnames(temp))]
	temp <- temp[sort(rownames(temp)),]
	temp
}

writeAsGenericMatrix <- function(theFile, theMatrix)
{
	write.table(theMatrix, file=theFile, quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)
	return(TRUE)
}

writeAsGenericMatrixNoRows <- function(theFile, theMatrix)
{
	write.table(theMatrix, file=theFile, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
	return(TRUE)
}

readAsGenericDataframe <- function(theFile, theNaString="NA")
{
	return( read.csv(theFile, header=TRUE, sep="\t", as.is=TRUE, check.names=FALSE, stringsAsFactors=FALSE, colClasses="character", na.strings=theNaString) )
}

writeAsGenericDataframe <- function(theFile, theDataframe)
{
	write.table(theDataframe, file=theFile, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
	return(TRUE)
}

################################################################################
