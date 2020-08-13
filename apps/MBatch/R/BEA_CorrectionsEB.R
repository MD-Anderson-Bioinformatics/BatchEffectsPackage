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

# need library call, even through parallel is part of base.
library(parallel)

BeaEB <- function(subMatrixGeneData, subDataframeBatchData, par.prior=TRUE, by='Batch', covariates=NULL,
								 prior.plots=FALSE, thePlotsFile="plots.png", theIssuesFile, theNumberOfThreads)
{
	logDebug("starting BeaEB")
	foo <- NULL
	tryCatch(
	foo <- EB(dat=subMatrixGeneData,
						si=convertDataFrameToSi(subDataframeBatchData),
						par.prior=par.prior,
						by=by,
						covariates=covariates,
						prior.plots=prior.plots,
						thePlotsFile=thePlotsFile,
						theNumberOfThreads=theNumberOfThreads)
	,error=function(e) {handleIssuesFunction(e, theIssuesFile) })
	logDebug("finishing BeaEB")
	return(foo)
}

doPriorPlots <- function(thePlotsFile, gamma.hat, gamma.bar, t2, delta.hat, a.prior, b.prior)
{
	logDebug(paste("Print prior plots at ", thePlotsFile ))
	CairoPNG(filename=thePlotsFile, width = 1600, height = 1000, pointsize=24, bg = "transparent")
	on.exit(dev.off(), add = TRUE)
	par(mfrow=c(2,2))
	tmp <- density(gamma.hat[1,], na.rm =TRUE)
	plot(tmp,  type='l', main="Density Plot")
	xx <- seq(min(tmp$x), max(tmp$x), length=100)
	lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
	qqnorm(gamma.hat[1,])
	qqline(gamma.hat[1,], col=2)

	tmp <- density(delta.hat[1,], na.rm =TRUE)
	invgam <- 1/rgamma(ncol(delta.hat),a.prior[1],b.prior[1])
	if(length(invgam)==(sum(is.na(invgam))))
	{
		plot.new()
		mtext("INVGAM is all NaN or NA")
	}
	else
	{
		tmp1 <- density(invgam, na.rm =TRUE)
		plot(tmp,  typ='l', main="Density Plot", ylim=c(0,max(tmp$y,tmp1$y)))
		lines(tmp1, col=2)
		qqplot(delta.hat[1,], invgam, xlab="Sample Quantiles", ylab='Theoretical Quantiles')
		lines(c(0,max(invgam)),c(0,max(invgam)),col=2)
		title('Q-Q Plot')
	}
	logDebug("finished prior plots")
}

EB <- function(dat, si, par.prior=TRUE, by='Batch', covariates=NULL, prior.plots=FALSE,
							thePlotsFile="plots.png", theNumberOfThreads)
{
	###############################
	### End of loading functions
	##############################
	logDebug("EB start")
  stopifnotWithLogging("Data sample ids must match batch sample ids and be in sorted order", all(colnames(dat)==rownames(si)))
  stopifnotWithLogging("requested batches and covariates should be in batch information", all(c(by, covariates) %in% colnames(si)))
	if(any(apply(dat,2,mode)!='numeric'))
	{
		return('ERROR: Array expression columns contain non-numeric values! (Check your data for non-numeric values)')
	}
	saminfo <- data.frame(sample=rownames(si), PBatch=si[,by], si[,covariates, drop=FALSE], check.names=FALSE)
	if(any(covariates != 'all'))
	{
		saminfo <- saminfo[,c('sample', 'PBatch' ,covariates)]
	}
	design <- design.mat(saminfo)
	batches <- list.batch(saminfo)
	n.batch <- length(batches)
	######check number of batch NZ 10/31/11
	logDebug("EB check number of batches")
	if(n.batch==1)
	{
		logWarn('There is only ONE batch, no adjustment made.')
		return(NULL)
	}
	n.batches <- sapply(batches, length)
	n.array <- sum(n.batches)
	### Check for missing values
	logDebug("EB Check for missing values")
	NAs = any(is.na(dat))
	if(NAs)
	{
		##logDebug(c('Found',sum(is.na(dat)),'Missing Data Values'),sep=' ')
	}

	##Check for genes with whole batch missing or no variation. added by NZ 10/31/2011
	logDebug("Check for genes with whole batch missing or no variation")
	gene.var<-apply(dat, 1, var, na.rm=TRUE)
	VAR0<-any(gene.var==0)
	if(NAs | VAR0)
	{
		missbatch<-apply(!is.na(dat), 1, function(x) tapply(x, INDEX=factor(saminfo[,'PBatch']), FUN=sum) )
		remove.gene.ind<-which(colSums(missbatch<1)!=0 | gene.var==0)
		if(length(remove.gene.ind)>0)
		{
			remove.gene<-dat[remove.gene.ind,,drop=FALSE]
			Original.dat<-dat
			dat<-Original.dat[-remove.gene.ind,]
		}
	}
	# if 0==nrow(s.data), then we can't do this -- this is checked for later
	if (0==nrow(dat))
	{
		logWarn('After checking for missing batches or no variation, no data is left, no adjustment made.')
		return(NULL)
	}
	###Standardize Data across genes
	logDebug('Standardizing Data across genes')
	B.hat <- NULL
	if (!NAs)
	{
		#logDebug('sd inside !NAs')
		B.hat <- solve(t(design)%*%design)%*%t(design)%*%t(as.matrix(dat))
	}
	else
	{
		#logDebug('sd inside !NAs else')
		### DEBUG this isn't generating an array/matrix/dataframe, but a list
		B.hat <- apply(dat,1,Beta.NA,design)
	}
	###Standarization Model
	logDebug('Standarization Model')
	grand.mean <- t(n.batches/n.array)%*%B.hat[1:n.batch,]
	if (!NAs)
	{
		#logDebug('sm inside !NAs')
		var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)
	}
	else
	{
		#logDebug('sm inside !NAs else')
		var.pooled <- apply(dat-t(design%*%B.hat),1,var,na.rm=TRUE)
	}
	logDebug('stand.mean')
	stand.mean <- t(grand.mean)%*%t(rep(1,n.array))
	if(!is.null(design))
	{
		#logDebug('!is.null(design)')
		tmp <- design;tmp[,c(1:n.batch)] <- 0;stand.mean <- stand.mean+t(tmp%*%B.hat)
	}
	#logDebug('s.data')
	s.data <- (dat-stand.mean)/(sqrt(var.pooled)%*%t(rep(1,n.array)))
	##Get regression batch effect parameters
	logDebug("Fitting L/S model and finding priors")
	batch.design <- design[,1:n.batch]
	if (!NAs)
	{
		logDebug("no NAs")
		# if 0==nrow(s.data), then we can't do this -- this is checked earlier
		gamma.hat <- solve(t(batch.design)%*%batch.design)%*%t(batch.design)%*%t(as.matrix(s.data))
	}
	else
	{
		logDebug("with NAs")
		gamma.hat=apply(s.data,1,Beta.NA,batch.design)
	}
	delta.hat <- NULL
	for (i in batches)
	{
		####logDebug(paste(i, collapse=", "), " gives ", paste(dim(s.data[,i]), collapse=", "))
		#To make sure no NA generated for delta.hat NZ. 12/1/11
		delta.hat.i<-apply(s.data[,i, drop=FALSE], 1, var,na.rm=TRUE)
		if(!all(is.na(delta.hat.i)))
		{
			delta.hat.i[is.na(delta.hat.i)]<-mean(delta.hat.i, na.rm=T)
		}
		delta.hat <- rbind(delta.hat, delta.hat.i)
	}
	#set priors for batches with only 1 sample. NZ 10/31/11
	delta.hat[is.na(delta.hat[,1]),]<-apply(delta.hat,2, mean, na.rm=TRUE)
	###Find Priors
	logDebug("Find priors")
	gamma.bar <- apply(gamma.hat, 1, mean, na.rm=TRUE)
	t2 <- apply(gamma.hat, 1, var, na.rm=TRUE)
	a.prior <- apply(delta.hat, 1, aprior)
	b.prior <- apply(delta.hat, 1, bprior)
	###Plot empirical and parametric priors
	logDebug("Plot empirical and parametric priors")
	if (prior.plots & par.prior)
	{
		tryCatch(
			doPriorPlots(thePlotsFile, gamma.hat, gamma.bar, t2, delta.hat, a.prior, b.prior),
			error=function(e)
			{
				logWarn("Unable to do prior plots",  sub("error", "issue", e, ignore.case=TRUE, fixed=TRUE) )
			})
	}
	###Find EB batch adjustments
	logDebug("Find EB batch adjustments")
	gamma.star <- delta.star <- NULL
	if(par.prior)
	{
		logDebug("Finding parametric adjustments")
		for (i in 1:n.batch)
		{
			temp <- it.sol(s.data[,batches[[i]], drop=FALSE],gamma.hat[i,],delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
			gamma.star <- rbind(gamma.star,temp[1,])
			delta.star <- rbind(delta.star,temp[2,])
		}
	}
	else
	{
		logDebug("Finding nonparametric adjustments")
		for (i in 1:n.batch)
		{
			temp <- int.eprior(as.matrix(s.data[,batches[[i]]]), gamma.hat[i,], delta.hat[i,], theNumberOfThreads)
			gamma.star <- rbind(gamma.star,temp[1,])
			delta.star <- rbind(delta.star,temp[2,])
		}

	}
	### Normalize the Data ###
	logDebug("Adjusting the Data")
	bayesdata <- s.data
	j <- 1
	for (i in batches)
	{
		bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
		j <- j+1
	}
	bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean
	###add back the removed genes with missing data in whole batch NZ 10/31/2011
	logDebug("add back the removed genes with missing data in whole batch")
	if(exists('Original.dat'))
	{
		Original.dat[-remove.gene.ind,]<-bayesdata
	}
	else
	{
		Original.dat<-bayesdata
	}
	logDebug("EB done")
	return(Original.dat)
}


##############################
###This file is a modified version of Combat function published by Johnson etal in Cheng Li's lab
###load functions
###############

### Next two functions make the design matrix (X) from the sample info file
build.design <- function(vec, des=NULL, start=2)
{
	tmp <- matrix(0,length(vec),nlevels(vec)-start+1)
	for (i in 1:ncol(tmp))
	{
		tmp[,i] <- vec==levels(vec)[i+start-1]
	}
	cbind(des,tmp)
}

design.mat <- function(saminfo)
{
	tmp <- which(colnames(saminfo) == 'PBatch')
	tmp1 <- as.factor(saminfo[,tmp])
	#logDebug("Found",nlevels(tmp1),'batches')
	design <- build.design(tmp1,start=1)
	ncov <- ncol(as.matrix(saminfo[,-c(1:2,tmp)]))
	#logDebug("Found",ncov,'covariate(s)')
	if(ncov>0)
	{
		for (j in 1:ncov)
		{
			tmp1 <- as.factor(as.matrix(saminfo[,-c(1:2,tmp)])[,j])
			design <- build.design(tmp1,des=design)
		}
	}
	design
}

### Makes a list with elements pointing to which array belongs to which batch
list.batch <- function(saminfo)
{
	tmp1 <- as.factor(saminfo[,which(colnames(saminfo) == 'PBatch')])
	batches <- NULL
	for (i in 1:nlevels(tmp1))
	{
		batches <- append(batches, list((1:length(tmp1))[tmp1==levels(tmp1)[i]]))
	}
	batches
}

### Trims the data of extra columns, note your array names cannot be named 'X' or start with 'X.'
trim.dat <- function(dat)
{
	tmp <- strsplit(colnames(dat),'\\.')
	tr <- NULL
	for (i in 1:length(tmp))
	{
		tr <- c(tr,tmp[[i]][1]!='X')
	}
	tr
}

### Following four find empirical hyper-prior values
aprior <- function(gamma.hat)
{
	m=mean(gamma.hat, na.rm=TRUE); s2=var(gamma.hat, na.rm=TRUE); (2*s2+m^2)/s2
}

bprior <- function(gamma.hat)
{
	m=mean(gamma.hat, na.rm=TRUE); s2=var(gamma.hat, na.rm=TRUE); (m*s2+m^3)/s2
}

postmean <- function(g.hat,g.bar,n,d.star,t2)
{
	(t2*n*g.hat+d.star*g.bar)/(t2*n+d.star)
}

postvar <- function(sum2,n,a,b)
{
	(.5*sum2+b)/(n/2+a-1)
}


### Pass in entire data set, the design matrix for the entire data, the batch means, the batch variances, priors (m, t2, a, b), columns of the data  matrix for the batch. Uses the EM to find the parametric batch adjustments

it.sol <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001)
{
n <- apply(!is.na(sdat),1,sum)
	g.old <- g.hat
	d.old <- d.hat
	change <- 1
	count <- 0
	#logDebug("it.sol - starting loop")
	#logDebug("it.sol - starting loop change=", change, " and conv=", conv)
	### DEBUG it appears that g.hat and d.hat are vectors. Which means g.net and d.new end up being vectors, and change ends up beings a vector, and change>conv ends up being a vector instead of a logical
	while(change>conv)
	{
		g.new <- postmean(g.hat,g.bar,n,d.old,t2)
		sum2 <- apply((sdat-g.new%*%t(rep(1,ncol(sdat))))^2, 1, sum, na.rm=TRUE)
		d.new <- postvar(sum2,n,a,b)
		change <- max(abs(g.new-g.old)/g.old,abs(d.new-d.old)/d.old, na.rm=TRUE)
		g.old <- g.new
		d.old <- d.new
		count <- count+1
			}
	#logDebug("it.sol - This batch took", count, "iterations until convergence")
	adjust <- rbind(g.new, d.new)
	rownames(adjust) <- c("g.star","d.star")
	adjust
}

###likelihood function used below
L <- function(x,g.hat,d.hat)
{
	prod(dnorm(x,g.hat,sqrt(d.hat)))
}

parallelIntEprior<-function(i, g.hat, d.hat, sdat)
{
	g <- g.hat[-i]
	d <- d.hat[-i]
	x <- sdat[i,!is.na(sdat[i,])]
	n <- length(x)
	j <- numeric(n)+1
	dat <- matrix(as.numeric(x),length(g),n,byrow=TRUE)
	resid2 <- (dat-g)^2
	sum2 <- resid2%*%j
	LH <- 1/(2*pi*d)^(n/2)*exp(-sum2/(2*d))
	LH[LH=="NaN"]=0
	g_star <- sum(g*LH)/sum(LH)
	d_star <- sum(d*LH)/sum(LH)
	return( c(g_star, d_star))
}

### Monte Carlo integration function to find the nonparametric adjustments
int.eprior <- function(sdat,g.hat,d.hat, theNumberOfThreads)
{
	g.star <- d.star <- NULL
	r <- nrow(sdat)
	myCluster <- makePSOCKcluster(theNumberOfThreads)
	#, outfile=file.path(getLogDir(), "int_eprior.log"))
	on.exit(stopCluster(myCluster))
	results <- parLapply(cl=myCluster, 1:r, parallelIntEprior, g.hat, d.hat, sdat)
	for(mypair in results)
	{
		g.star <- c(g.star, mypair[[1]])
		d.star <- c(d.star, mypair[[2]])
	}
	adjust <- rbind(g.star,d.star)
	rownames(adjust) <- c("g.star","d.star")
	adjust
}

### fits the L/S model in the presence of missing data values
Beta.NA = function(y,X)
{
	#logDebug("Beta.NA - start")
	des=X[!is.na(y),]
	y1=y[!is.na(y)]
	## debug - for some data this give us system is exactly singular
	B <- solve(t(des)%*%des)%*%t(des)%*%y1
	#logDebug("Beta.NA - end (before return)")
	return(B)
}
