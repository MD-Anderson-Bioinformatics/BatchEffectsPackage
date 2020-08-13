# MBatch Copyright (c) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020 University of Texas MD Anderson Cancer Center
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
# MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>

### Log a message for the developer 
logDebug<-function(...)
{
	logOutput(theLevelName="DEBUG", ... )
}

### Log a timing message
logTiming<-function(theAlgorithm, theData, theSystemTime)
{
	logOutput(theLevelName="TIMING", "\t", theSystemTime[1], "\t", theSystemTime[3], "\t", theAlgorithm, "\t", theData )
}

### Log item at in a plain text info message in English that could be easily understood by a user.  
logInfo<-function(...)
{
	logOutput(theLevelName="INFO", ... )
}

### Log a warning message.  i.e. If an number should be within a certain range to be 
### effective and it is outside the range, display a warning message.
logWarn<-function(...)
{
	logOutput(theLevelName="WARN", ... )
}

### Log message telling the percentage of the program that has completed. 
logPercent<-function(...)
{
	logOutput(theLevelName="PERCENT", ... )
}

### Log an error.
logError<-function(...)
{
	logOutput(theLevelName="ERROR", ... )
}

getLogDir <- function()
{
	logger<-get("logger", inherits=TRUE)
	if (is.null(logger)==TRUE)
	{
		return("")
	}
	else
	{
		dirname(logger@mFile)
	}
}

### Write to the console and/or the default log file or the console and/or a specified log file
logOutput<-function(theLevelName="INFO", ..., theLogFile=NULL)
{
	if (FALSE==exists("logger", inherits=TRUE))
	{
		assign("logger", NULL, inherits=TRUE)
	}
	logger<-get("logger", inherits=TRUE)
	if (is.null(logger)==TRUE)
	{
		logger<-new("Logging")
		assign("logger", logger, inherits=TRUE)
	}
	if (length(which(logger@mLevelNames == theLevelName))>0)
	{
		if (length(which(logger@mLevelNamesToLog == theLevelName))>0)
		{
			### Write a time stamp, log.seperator, then all the other entries sent in with log.seperator between them
			### and add a newline at the end
			logFileName <- NULL
			if (!is.null(theLogFile)) 
			{
				logFileName = theLogFile
			}
			else if (logger@mFile!="")
			{
				logFileName <- logger@mFile
			}
			
			
			if(is.null(logFileName) || logFileName=="" || logger@mConsole == TRUE) 
			{
				####cat(paste(format(Sys.time(), "%Y %m %d %H:%M:%OS3"), theLevelName, Sys.info()['nodename'], paste(lapply(..., function(x) (return(paste(x, collapse=", ")))), collapse="; "), sep=logger@mSeparator), "\n" )
				cat(paste(format(Sys.time(), "%Y %m %d %H:%M:%OS3"), theLevelName, Sys.info()['nodename'], paste(..., collapse=", "), sep=logger@mSeparator), "\n" )
			}
			if (!is.null(logFileName) && logFileName!="")
			{
				###cat(paste(format(Sys.time(), "%Y %m %d %H:%M:%OS3"),theLevelName, Sys.info()['nodename'], paste(lapply(..., function(x) (return(paste(x, collapse=", ")))), collapse="; "), sep=logger@mSeparator), "\n", file=logFileName, append=TRUE)
				cat(paste(format(Sys.time(), "%Y %m %d %H:%M:%OS3"),theLevelName, Sys.info()['nodename'], paste(..., collapse=", "), sep=logger@mSeparator), "\n", file=logFileName, append=TRUE)
			}
		}
	}
}

setLogging<-function(theLogger)
{
	assign("logger", theLogger, inherits=TRUE);
}

stopWithLogging<-function(msg )
{
	###last.function.name <- as.list(sys.call(which=-1)[1])
 	###full.message = paste(last.function.name, "():", msg, sep="")
	###logError(full.message)
 	###stop(full.message, call.=FALSE)
	logError(msg)
	stop(msg, call.=FALSE)
}

stopifnotWithLogging<-function(msg="", ... )
{
	if (sum(...)!=length(c(...))) 
	{
		stopWithLogging(msg)	
	}
}

warnifnotWithLogging<-function(msg="", ... )
{
	if (sum(...)!=length(c(...))) 
	{
		logWarn(msg)	
	}
}

###writeTraceBack<-function(theErrorFile, theTraceback)
###{
###	sink(file=theErrorFile, append=TRUE, split=TRUE)
###	traceback(theTraceback)
###	sink()
###}

handleIssuesFunction<-function(theError, theExtraErrorFile=NULL)
{
	dump.frames()
	n <- length(last.dump)
	calls <- names(last.dump)
	errorString <- paste(paste(sub("error", "issue", theError, ignore.case=TRUE, fixed=TRUE)), paste(paste("  ", 1L:n, ": ", calls, sep = "", collapse= "\n"), collapse= "\n"), "\n\n", sep = "\n")
	logError(errorString)
	if (!is.null(theExtraErrorFile))
	{
		checkCreateDir(dirname(theExtraErrorFile))
		cat(errorString, file=theExtraErrorFile, append=TRUE)
		cat("\n", file=theExtraErrorFile, append=TRUE)
	}
}
