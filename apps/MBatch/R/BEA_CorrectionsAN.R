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

BeaAN<- function(subMatrixGeneData, subDataframeBatchData, by='Batch', var.adj=TRUE, theIssuesFile)
{
	logDebug("starting BeaAN")
	foo <- NULL
	tryCatch(
		foo <- AN(subMatrixGeneData, convertDataFrameToSi(subDataframeBatchData), by, var.adj)
	,
		error=function(e)
		{
			handleIssuesFunction(e, theIssuesFile)
		})
	logDebug("finishing BeaAN")
	return(foo)
}

#function to calculate
adj.NA <- function(y,X,dat)
{
	#logDebug("adj.NA des")
	des=X[!is.na(y),]
	#logDebug("adj.NA y1")
	y1=y[!is.na(y)]
	#logDebug("adj.NA B.hat")
	B.hat <- solve(t(des)%*%des)%*%t(des)%*%y1
	#logDebug("adj.NA B.bar")
	B.bar<-mean(B.hat)
	#logDebug("adj.NA resid")
	resid<-y1-des %*% B.hat
	pool.var<-sum(resid^2/ncol(dat))
	#get the variance of each batch of each gene as the weight
	ns.batch<-apply(des, 2, sum)
	average<-sweep(x=des,2, ns.batch, "/")
	batchvar<-t(resid^2) %*%average
	#get the pooled weight for each batch of sample
	batchvar.pool<-as.vector(batchvar %*% t(des))
	batchvar.pool[batchvar.pool==0]<-mean(batchvar.pool)
	# TODO: this is where singular check goes
	#logDebug("adj.NA - this is where singular check goes")
	if(0==det(diag(batchvar.pool)))
	{
		#logWarn("System is singular--no correction made. Skip further processing")
		y <- NULL
	}
	else if(rcond((t(des)%*%solve(diag(batchvar.pool))%*%des))<.Machine$double.eps)
	{
			#logWarn("System is singular--no correction made. Skip further processing")
			y <- NULL
	}
	else
	{
		#logDebug("adj.NA - solve diag")
		weight<-solve(diag(batchvar.pool))
		#logDebug("adj.NA - solve weight")
		foo <- solve(t(des)%*%weight%*%des)
		bar <- foo%*%t(des)
		qwerty <- bar%*%weight
		B.hat.w <- qwerty%*%y1
		#logDebug("adj.NA - colMeans")
		B.bar.w<-colMeans(B.hat.w)
		resid.w<-y1-des %*% B.hat.w
		#adjust for variance of the residuals
		varadj<-sqrt(pool.var/t(batchvar %*% t(des)))
		resid.w<-resid.w*varadj
		y[!is.na(y)]<-B.bar.w+resid.w
		##logDebug("adj.NA post Y =", paste(y,collapse=","))
	}
	return(y)
}

AN<-function(dat, si, by='Batch', var.adj=TRUE)
{
	ANdat <- NA
	logDebug("AN names")
	stopifnotWithLogging("samples should be same and in same order in data and batch information", all(colnames(dat)==rownames(si)))
	logDebug("AN all")
	stopifnotWithLogging("All requested batch types should be in batch information", all(by %in% colnames(si)))
	logDebug("AN cbin")
	batchInfo<-cbind(rownames(si), si[,by])
	#dat is a matrix with gene in row sample in column, batchInfo is a vector indicates the batch, var.adj=FALSE assumes the variance in each batch are the same)
	#a function to build design matrix of batch
	logDebug("AN function")
	build.X <- function(vec, des=NULL, start=1)
	{
		logDebug("AN build.X")
		vec<-as.factor(vec)
		tmp <- matrix(0,length(vec),nlevels(vec)-start+1)
		for (i in 1:ncol(tmp))
		{
			tmp[,i] <- vec==levels(vec)[i+start-1]
		}
		cbind(des,tmp)
	}
	######check number of batch NZ 12/6/11
	logDebug("AN check number of batch")
	n.batch<-nlevels(factor(si[,by]))
	if(n.batch==1)
	{
		logWarn('There is only ONE batch, no adjustment made.')
		return(NULL)
	}
	### Check for missing values
	logDebug("AN Check for missing values")
	NAs = any(is.na(dat))
	##Check for genes with whole batch missing or no variation. added by NZ 10/31/2011
	logDebug("AN Check for genes with whole batch missing or no variation.")
	gene.var<-apply(dat, 1, stats::var, na.rm=T)
	VAR0<-any(gene.var==0)
	if(NAs | VAR0)
	{
		missbatch<-apply(!is.na(dat), 1, function(x) tapply(x, INDEX=factor(si[,by]), FUN=sum) )
		remove.gene.ind<-which(colSums(missbatch<1)!=0 | gene.var==0)
		if(length(remove.gene.ind)>0)
		{
			remove.gene<-dat[remove.gene.ind,,drop=F]
			Original.dat<-dat
			dat<-Original.dat[-remove.gene.ind,]
		}
	}
	#get the design matrix for batch
	logDebug("AN design")
	design<-build.X(si[,by])
	#get the LS estimates
	logDebug("AN NAs")
	### Check for missing values
	NAs = any(is.na(dat))
	if (!NAs & var.adj)
	{
		logDebug("AN B.hat")
		B.hat <- solve(t(design)%*%design)%*%t(design)%*%t(as.matrix(dat))
		logDebug("AN B.bar")
		B.bar<-apply(B.hat, 2, mean)
		logDebug("AN resid")
		resid<-dat-t(design %*% B.hat)
		#get the pooled variance for each gene
		logDebug("AN pool.var")
		pool.var<-resid^2 %*% rep(1/ncol(dat),ncol(dat))
		#get the variance of each batch of each gene as the weight
		logDebug("AN ns.batch")
		ns.batch<-apply(design, 2, sum)
		logDebug("AN average")
		average<-sweep(x=design,2, ns.batch, "/")
		logDebug("AN batchvar")
		batchvar<-resid^2 %*%average
		#get the pooled weight for each batch of sample
		logDebug("AN batchvar.pool")
		batchvar.pool<-colMeans(batchvar %*% t(design), na.rm=TRUE)
		logDebug("AN batchvar.pool")
		batchvar.pool[batchvar.pool==0]<-mean(batchvar.pool)
		logDebug("AN weight")
		if(0==det(diag(batchvar.pool)))
		{
			logWarn("det(diag(batchvar.pool)) --> System is singular--no correction made. Skip further processing")
			return(NULL)
		}
		weight<-solve(diag(batchvar.pool))
		logDebug("AN B.hat.w")
		B.hat.w <- solve(t(design)%*%weight%*%design)%*%t(design)%*%weight%*%t(as.matrix(dat))
		logDebug("AN B.bar.w")
		B.bar.w<-colMeans(B.hat.w)
		logDebug("AN resid.w")
		resid.w<-dat-t(design %*% B.hat.w)
		#adjust for variance of the residuals
		logDebug("AN varadj")
		varadj<-sqrt(matrix(pool.var, ncol=ncol(resid), nrow=nrow(batchvar))/(batchvar %*% t(design)))
		logDebug("AN resid.w")
		resid.w<-resid.w*varadj
		logDebug("AN ANdat")
		ANdat<-matrix(B.bar.w, nrow=nrow(dat), ncol=ncol(dat))+resid.w
		logDebug("AN after ANDat")
	}
	else if(NAs & var.adj)
	{
		logInfo("NAs & var.adj")
		logInfo("is.matrix dat", is.matrix(dat))
		ANdat <- apply(dat, 1, adj.NA, design, dat)
		logInfo("transpose")
		logInfo("is.matrix ANdat", is.matrix(ANdat))
		if (is.matrix(ANdat))
		{
			ANdat <- t(ANdat)
			logInfo("check nulls")
			nullTrue <- apply(ANdat, c(1,2), function(x) { is.null(unlist(x)) } )
			logInfo("sum nulls")
			if (sum(nullTrue)>0)
			{
				logWarn("System is singular--no correction made. Skip further processing")
				ANdat <- NULL
			}
		}
		else
		{
			logWarn("ANdat is not a matrix--no correction made. Skip further processing")
			ANdat <- NULL
		}
	}
	else
	{
		B.hat <- solve(t(design)%*%design)%*%t(design)%*%t(as.matrix(dat))
		B.bar<-apply(B.hat, 2, mean)
		resid<-dat-t(design %*% B.hat)
		ANdat<-matrix(B.bar, nrow=nrow(dat), ncol=ncol(dat))+resid
	}
	###add back the removed genes with missing data in whole batch NZ 10/31/2011
	if (is.null(ANdat))
	{
		logWarn("System is singular. Skip further processing")
		Original.dat <- NULL
	}
	else if(exists('Original.dat'))
	{
		Original.dat[-remove.gene.ind,]<-ANdat
	}
	else
	{
		Original.dat<-ANdat
	}
	return(Original.dat)
}
