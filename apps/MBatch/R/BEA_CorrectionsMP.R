#MBatch Copyright ? 2011, 2012, 2013, 2014, 2015, 2016, 2017 University of Texas MD Anderson Cancer Center
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

BeaMP<- function(subMatrixGeneData, subDataframeBatchData, by='Batch', overall=TRUE, theIssuesFile=NULL)
{
	logDebug("starting BeaMP")
	foo <- NULL
	tryCatch(
	foo <- MP(subMatrixGeneData, convertDataFrameToSi(subDataframeBatchData), by, overall)
	,error=function(e) {handleIssuesFunction(e, theIssuesFile) })
	logDebug("finishing BeaMP")
	return(foo)
}

MP<-function(dat, si, by='Batch', overall=FALSE, ...)
{
	logDebug("starting MP")
	if(missing(si) & overall==FALSE)
	{
		stop('sample information is needed for batch-wise median polish corrections')
	}
	if(overall==TRUE)
	{
		logDebug("MP overall")
		MPdat<-medpolish(dat, eps=0.0001, trace.iter=FALSE, na.rm=TRUE, ...)
		MPdat<-MPdat$residuals+MPdat$overall
		final<-MPdat
	}
	else
	{
		logDebug("MP batch")
		stopifnotWithLogging("Data sample names should match and be in same order as those for batch data", all(colnames(dat)==rownames(si)))
		stopifnotWithLogging("All requested batch types should be in batch data", all(by %in% colnames(si)))
		batch<-table(si[,by])
		###logDebug("MP batch 1")
		MPBdat<-dat
		###logDebug("MP batch 2")
		MPBcol<-rep(NA, ncol(dat))
		###logDebug("MP batch 3")
		for(i in names(batch))
		{
			###logDebug("MP batch 4")
			temp<-dat[, si[,by]==i]
			###logDebug("MP batch 5")
			temp.MP<-medpolish(temp,eps=0.0001, trace.iter=FALSE,na.rm=TRUE, ...)
			###logDebug("MP batch 6")
			MPBdat[, si[,by]==i]<-temp.MP$residuals
			###logDebug("MP batch 7")
		}
		###logDebug("MP batch 8")
		all.MP<-medpolish(dat, eps=0.0001, trace.iter=FALSE,na.rm=TRUE)
		###logDebug("MP batch 9")
		final<-MPBdat+all.MP$overall
		###logDebug("MP batch 10")
	}
	return(final)
}