# MBatch Copyright (c) 2011-2022 University of Texas MD Anderson Cancer Center
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

cleanAndShortPath <- function(theDir)
{
  theDir <- gsub("[^[:alnum:]]", "", theDir)
  if (nchar(theDir)>32)
  {
    milliseconds <- substr(format(as.numeric(Sys.time())*100000, digits=15), 11, 15)
    theDir <- paste(substr(theDir, 0, 32), milliseconds, sep="", collapse="")
  }
  theDir
}

cleanFilePath <- function(theDir, theNewDir)
{
  #cleanDir <- iconv(theNewDir, from = "UTF-8", to = "ASCII", sub = "byte")
  #cleanDir <- gsub("><", "_", cleanDir, fixed=TRUE)
  #cleanDir <- gsub("<", "_", cleanDir, fixed=TRUE)
  #cleanDir <- gsub(">", "_", cleanDir, fixed=TRUE)
  #file.path(theDir, cleanDir)
  tmpPath <- file.path(theDir, theNewDir)
  # clean double // from the path
  tmpPath <- gsub("//", "/", tmpPath)
  tmpPath
}

createDirPlusFilename<-function(theDir, ...)
{
  cleanFilePath(theDir, paste(..., sep="", collapse=""))
}

addVersionsIfNeeded <- function(thePath, theDataVersion, theTestVersion)
{
  if (!is.null(thePath))
  {
    if ((!is.null(theDataVersion)) && (""!=theDataVersion))
    {
      thePath <- file.path(thePath, theDataVersion)
    }
    if ((!is.null(theTestVersion)) && (""!=theTestVersion))
    {
      thePath <- file.path(thePath, theTestVersion)
    }
  }
  thePath
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
  return(thePath)
}

checkCreateDir<-function(theBaseDir, theNewDir)
{
	myDir <-cleanFilePath(theBaseDir, theNewDir)
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
								 sep="\t", na.strings="", flush=TRUE, fill=FALSE, multi.line=FALSE,
								 allowEscapes=TRUE)
	genes <- as.vector(unlist(myScan[1]))
	data <- as.vector(unlist(myScan[2:length(myScan)]))
	temp<-matrixWithIssues(data,
												 ncol=length(samples),
												 nrow=length(genes),
												 dimnames=list(make.unique(genes), make.unique(samples)),
												 )
	temp <- temp[,sort(colnames(temp)),drop=FALSE]
	temp <- temp[sort(rownames(temp)),,drop=FALSE]
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

readAsGenericDataframe <- function(theFile, theNaString=NULL, theUnknownString="Unknown")
{
	df <- read.csv(theFile, header=TRUE, sep="\t", as.is=TRUE, check.names=FALSE,
	               stringsAsFactors=FALSE, colClasses="character", na.strings=theNaString,
	               allowEscapes=TRUE, quote="" )
	# check for colons in the batch type names (incompatible with MANOVA tests)
	coln  <- colnames(df)
	coln <- unlist(lapply(coln, function(theName)
	{
	  if (grepl(":", theName, fixed = TRUE))
	  {
	    logWarn("Replacing colon in batch type name ", theName, " with tilde")
	    theName <- gsub(":", "~", theName, fixed = TRUE)
	  }
	  theName
	}))
	colnames(df) <- coln
	df[df==""] <- theUnknownString
  return(df)
}

writeAsGenericDataframe <- function(theFile, theDataframe)
{
	write.table(theDataframe, file=theFile, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
	return(TRUE)
}

################################################################################
