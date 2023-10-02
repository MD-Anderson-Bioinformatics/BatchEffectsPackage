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


openAndWriteIssuesLogFileSuperClust<-function(theOutputDir)
{
  myFile <- file(cleanFilePath(theOutputDir, "error.log"), "w+")
  on.exit(close(myFile))
  cat("Unable to Generate Supervised Clustering Heatmap\n", file=myFile, append=TRUE)
}

###########################################################################################
###########################################################################################

createBatchEffectsOutput_SupervisedClustering_batches<-function(theMatrixGeneData, theDataframeBatchData,
                                                                theDataVersion, theTestVersion,
																																theTitle, theOutputDir)
{
  lastPath <- ""
	tryCatch({
	  for(batchTypeIndex in c(2:length(theDataframeBatchData)))
  	{
  		### compile data and information for display
  		batchTypeName <- names(theDataframeBatchData)[batchTypeIndex]
  		myOutputPath <- cleanFilePath(theOutputDir, batchTypeName)
  		myOutputPath <- addVersionsIfNeeded(myOutputPath, theDataVersion, theTestVersion)
  		lastPath <- myOutputPath
  		checkDirForCreation(myOutputPath)
  		# minor slow down, but needs to be here so path and error.log work
  		checkIfTestError()
  		batchIdsForSamples <- as.character(as.vector(unlist(theDataframeBatchData[batchTypeIndex])))
  		logInfo("createBatchEffectsOutput_SupervisedClustering_batches - batchTypeName = ", batchTypeName)
  		color <- beaRainbow(length(unique(batchIdsForSamples)), v=0.7)
  		title <- paste(theTitle, "SupervisedClust", batchTypeName, sep=" ")
  		diagramFilename = createDirPlusFilename(myOutputPath, "SupervisedClust_Diagram-", batchTypeName, ".png")
  		legendFilename = createDirPlusFilename(myOutputPath, "SupervisedClust_Legend-", batchTypeName, ".png")
  		makeBiasClust(theMatrixGeneData, theDataframeBatchData, batchTypeName, color, title, diagramFilename)
  		supervisedClustLegend(color, batchIdsForSamples, sort(unique(batchIdsForSamples)), batchTypeName, legendFilename)
	  }
	},
	error=function(e)
	{
	  logWarn("2 Unable to generate SupervisedClustering_batches")
	  logDebug(lastPath)
	  unlink(lastPath, recursive = TRUE, force = TRUE)
	  dir.create(lastPath, showWarnings=FALSE, recursive=TRUE)
	  openAndWriteIssuesLogFileSuperClust(lastPath)
	})
}

createBatchEffectsOutput_SupervisedClustering_pairs<-function(theMatrixGeneData, theDataframeBatchData, thePairList,
																															theDataVersion, theTestVersion,
																															theTitle, theOutputDir)
{
  lastPath <- ""
	for(indexA in seq(1, length(thePairList), 2))
	{
		indexB <- indexA + 1
		### compile data and information for display
		batchTypeNameA <- thePairList[indexA]
		batchIdsForSamplesA <- as.character(as.vector(unlist(theDataframeBatchData[batchTypeNameA])))
		batchTypeNameB <- thePairList[indexB]
		batchIdsForSamplesB <- as.character(as.vector(unlist(theDataframeBatchData[batchTypeNameB])))
		combinedName = paste(batchTypeNameA, "-", batchTypeNameB, sep="")
		logInfo("createBatchEffectsOutput_SupervisedClustering_pairs - batchTypeNameA = ", batchTypeNameA)
		logInfo("createBatchEffectsOutput_SupervisedClustering_pairs - batchTypeNameB = ", batchTypeNameB)
		myOutputPath <- cleanFilePath(theOutputDir, combinedName)
		myOutputPath <- addVersionsIfNeeded(myOutputPath, theDataVersion, theTestVersion)
		checkDirForCreation(myOutputPath)
		lastPath <- myOutputPath
		tryCatch({
		  checkIfTestError()
  		color <- beaRainbow(max(length(unique(batchIdsForSamplesA)),length(unique(batchIdsForSamplesB))), v=0.7)
  		title <- paste(theTitle, combinedName, sep=" ")
  		diagramFilename = createDirPlusFilename(myOutputPath, "SupervisedClust_Diagram.png")
  		legendFilenameA = createDirPlusFilename(myOutputPath, "SupervisedClust_Legend", "-", batchTypeNameA,".png")
  		legendFilenameB = createDirPlusFilename(myOutputPath, "SupervisedClust_Legend", "-", batchTypeNameB,".png")
  		logInfo("createBatchEffectsOutput_SupervisedClustering_pairs - pre bias clust")
  		success <- makeBiasClust(theMatrixGeneData, theDataframeBatchData, batchTypeNameA, color, title,
  									diagramFilename, theMore=batchTypeNameB)
  		logInfo("createBatchEffectsOutput_SupervisedClustering_pairs - post bias clust")
  		supervisedClustLegend(color, batchIdsForSamplesA, sort(unique(batchIdsForSamplesA)), batchTypeNameA, legendFilenameA)
  		logInfo("createBatchEffectsOutput_SupervisedClustering_pairs - after first legend")
  		supervisedClustLegend(color, batchIdsForSamplesB, sort(unique(batchIdsForSamplesB)), batchTypeNameB, legendFilenameB)
  		logInfo("createBatchEffectsOutput_SupervisedClustering_pairs - after second legend")
		},
		error=function(e)
		{
		  logWarn("2 Unable to generate SupervisedClustering_pairs")
		  unlink(lastPath, recursive = TRUE, force = TRUE)
		  dir.create(lastPath, showWarnings=FALSE, recursive=TRUE)
		  openAndWriteIssuesLogFileSuperClust(lastPath)
		})
	}
}

###########################################################################################
###########################################################################################

#This function is for biasedcluster algorithms
biasedClusterFunction<-function(Dat, cat)
{
	stopifnotWithLogging("Samples should be same length in data and cat", ncol(Dat)==length(cat))
	cat<-factor(cat, exclude=NULL)
	Dat <- Dat[!is.na(rowSums(Dat)),]
	Dat <- Dat[!is.infinite(rowSums(Dat)),]
	if (0==nrow(Dat))
	{
		return(NULL)
	}
	dis<-as.matrix(distanceMatrix(Dat, metric='pearson'))
	max.dis<-max(dis)
	n.fact<-as.numeric(factor(cat, exclude=NULL))
	nl<-max(n.fact)
	dis[outer(n.fact, n.fact, '!=')]<-max.dis
	new.dis<-as.dist(dis)
	logInfo("makeBiasClust new.dis size - ", paste(dim(as.matrix(new.dis)), sep=",", collapse="-"))
	new.dis<-as.matrix(new.dis)
	logInfo("makeBiasClust orig - ", paste(dim(new.dis), sep=",", collapse="-"))
	new.dis <- new.dis[!is.na(rowSums(new.dis)),]
	logInfo("makeBiasClust is.na - ", paste(dim(new.dis), sep=",", collapse="-"))
	# NOT NEEDED as is.na also checks for NaN
	#new.dis <- new.dis[!is.nan(rowSums(new.dis)),]
	#logInfo("makeBiasClust is.nan - ", paste(dim(new.dis), sep=",", collapse="-"))
	new.dis <- new.dis[!is.infinite(rowSums(new.dis)),]
	logInfo("makeBiasClust is.infinite - ", paste(dim(new.dis), sep=",", collapse="-"))
	if (nrow(new.dis)>0)
	{
		new.dis <- as.dist(new.dis)
		hc<-hclust(new.dis, method='average')
		maxn<-sort(hc$height, decreasing=TRUE)[1:nl]
		hc$height[hc$height>min(maxn)]<-median(maxn)
		return(hc)
	}
	else
	{
		return(NULL)
	}
}

udataCalc <- function(dat, Mad)
{
	U.data <- NULL
	if(nrow(dat)>1000)
	{
		logInfo("makeBiasClust - quantile dat dim = ", paste(dim(dat),collapse=","))
		U.data<-dat[!is.na(Mad) & Mad> quantile(Mad, na.rm=TRUE)[4],]
		logInfo("makeBiasClust - quantile U.data is.data.frame = ", paste(is.data.frame(U.data),collapse=","))
		logInfo("makeBiasClust - quantile U.data is.array = ", paste(is.array(U.data),collapse=","))
		logInfo("makeBiasClust - quantile U.data is.list = ", paste(is.list(U.data),collapse=","))
		logInfo("makeBiasClust - quantile U.data nrow = ", paste(nrow(U.data),collapse=","))
		logInfo("makeBiasClust - quantile U.data ncol = ", paste(ncol(U.data),collapse=","))
		logInfo("makeBiasClust - quantile U.data length = ", paste(length(U.data),collapse=","))
		logInfo("makeBiasClust - quantile U.data dim = ", paste(dim(U.data),collapse=","))
		logInfo("makeBiasClust - quantile U.data is.null = ", is.null(U.data))
		if (FALSE==is.array(U.data))
		{
			logWarn("makeBiasClust - Too many NA values generated during correction")
			U.data <- NULL
		}
	}
	else
	{
		logInfo("makeBiasClust - dat dim = ", dim(dat))
		U.data<-dat[!is.na(Mad),]
		logInfo("makeBiasClust - U.data dim = ", dim(U.data))
	}
	return(U.data)
}

makeBiasClust<-function(theMatrix, theBatches, theBatchType, theColors, theTitle,
                        theDiagramFilename, theMore=NULL)
{
	logInfo("makeBiasClust - starting")
	stopifnotWithLogging("Samples should be same length in data and batch info", ncol(theMatrix)==nrow(theBatches))
	stopifnotWithLogging("Batch types requested should be in list of batch types", theBatchType %in% colnames(theBatches))
	##
	## select samples with data and no NAs
	##
	keepSamples <- colSums(!is.na(theMatrix))!=0
	theMatrix <- theMatrix[,keepSamples]
	theBatches <- theBatches[keepSamples,]
	##
	## add second batch type (theMore) to batches to use
	##
	theMore <- intersect(theMore, colnames(theBatches))
	theBatches <- theBatches[,c(theBatchType, theMore), drop=FALSE]
	##
	## get "interesting" in term of batch effects data
	##
	madValues <- apply(theMatrix, 1, mad)
	interestingData <- udataCalc(theMatrix, madValues)
	##
	## get batch infor for selected batch types
	##
	selectedBatches <- data.frame(theBatches, stringsAsFactors=FALSE, check.names=FALSE)
	#prepare matrix for side color matrix
	###Do biased hierarchical clustering
	logInfo("makeBiasClust - biasedDend <- biasedClusterFunction")
	if (!is.null(interestingData))
	{
	  biasedDend <- NULL
	  tryCatch(
	    {
	      biasedDend <- biasedClusterFunction(interestingData, factor(selectedBatches[,theBatchType]))
	    },
	    error=function(e)
	    {
	      logWarn("2 Unable to calculate biasedClusterFunction--too many NAs, Infinities or NaNs in data")
	      biasedDend <- NULL
	    })
		featureDend <- NULL
		if (!is.null(biasedDend))
		{
		  featureDend <- hierClust_calc(t(interestingData))
		}
		if (!is.null(featureDend))
		{
  		if (!is.null(biasedDend))
  		{
  			############################################
  			preFilename <- theBatchType
  			if (length(theMore)>0)
  			{
  			  preFilename <- paste(theBatchType, theMore, sep="-")
  			}
  			############################################
  			writeHCDataTSVs(biasedDend, dirname(theDiagramFilename),
  			                paste(preFilename, "HCData.tsv", sep="_"),
  			                paste(preFilename, "HCOrder.tsv", sep="_"),
  			                paste(preFilename, "uDend.RData", sep="_"))
  			writeHCDataTSVs(featureDend, dirname(theDiagramFilename),
  			                paste(preFilename, "HCData_feature.tsv", sep="_"),
  			                paste(preFilename, "HCOrder_feature.tsv", sep="_"),
  			                paste(preFilename, "uDend_feature.RData", sep="_"))
  			writeSCmatrix(interestingData, dirname(theDiagramFilename),
  			              paste(preFilename, "SCMatrix.tsv", sep="_"),
  			              theTitle,
  			              paste(preFilename, "title.tsv", sep="_"))
  			############################################
  			CairoPNG(filename=theDiagramFilename, width = 1000, height = 1000, pointsize=24, bg = "transparent")
  			on.exit(dev.off(), add = TRUE)
  			heatmap(interestingData, Rowv = as.dendrogram(featureDend), Colv = as.dendrogram(biasedDend))
  		}
  		else
  		{
  		  openAndWriteIssuesLogFileSC(dirname(theDiagramFilename))
  		}
		}
		else
		{
		  openAndWriteIssuesLogFileSC(dirname(theDiagramFilename))
		}
	}
	else
	{
	  openAndWriteIssuesLogFileSC(dirname(theDiagramFilename))
	}
}

###########################################################################################
###########################################################################################

heatmap.plus.mod<-function(
	x,
	Rowv = NULL,
	Colv = if (symm) "Rowv" else NULL,
	distfun = dist,
	hclustfun = hclust,
	reorderfun = function(d, w) reorder(d, w),
	add.expr,
	symm = FALSE,
	revC = identical(Colv, "Rowv"),
	scale = c("row", "column", "none"),
	na.rm = TRUE,
	margins = c(5, 5),
	ColSideColors,
	RowSideColors,
	cexRow = 0.2 + 1/log10(nr),
	cexCol = 0.2 + 1/log10(nc),
	labRow = NULL,
	labCol = NULL,
	main = NULL,
	xlab = NULL,
	ylab = NULL,
	keep.dendro = FALSE,
	SideColorHeight=1,
	Colsep=NULL,
	sepColor='white',
	sepLWD=2,
	verbose = getOption("verbose") )
{
	scale <- NULL
	if (symm && missing(scale))
	{
		scale <- "none"
	}
	else
	{
		scale <- match.arg(scale)
	}
	if (length(di <- dim(x)) != 2 || !is.numeric(x))
	{
		stop("'x' must be a numeric matrix")
	}
	nr <- di[1]
	nc <- di[2]
	if (nr <= 1 || nc <= 1)
	{
		stop("'x' must have at least 2 rows and 2 columns")
	}
	if (!is.numeric(margins) || length(margins) != 2)
	{
		stop("'margins' must be a numeric vector of length 2")
	}
	doRdend <- !identical(Rowv, NA)
	doCdend <- !identical(Colv, NA)
	if (is.null(Rowv))
	{
		Rowv <- rowMeans(x, na.rm = na.rm)
	}
	if (is.null(Colv))
	{
		Colv <- colMeans(x, na.rm = na.rm)
	}
	if (doRdend)
	{
		if (inherits(Rowv, "dendrogram"))
		{
			ddr <- Rowv
		}
		else
		{
			hcr <- hclustfun(distfun(x))
			ddr <- as.dendrogram(hcr)
			if (!is.logical(Rowv) || Rowv)
			{
				ddr <- reorderfun(ddr, Rowv)
			}
		}
		if (nr != length(rowInd <- order.dendrogram(ddr)))
		{
			stop("row dendrogram ordering gave index of wrong length")
		}
	}
	else
	{
		rowInd <- 1:nr
	}
	if (doCdend)
	{
		if (inherits(Colv, "dendrogram"))
		{
			ddc <- Colv
		}
		else if (identical(Colv, "Rowv"))
		{
			if (nr != nc)
			{
				stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
			}
			ddc <- ddr
		}
		else
		{
			hcc <- hclustfun(distfun(if (symm) x else t(x)))
			ddc <- as.dendrogram(hcc)
			if (!is.logical(Colv) || Colv)
			{
				ddc <- reorderfun(ddc, Colv)
			}
		}
		if (nc != length(colInd <- order.dendrogram(ddc)))
		{
			stop("column dendrogram ordering gave index of wrong length")
		}
	}
	else
	{
		colInd <- 1:nc
	}
	x <- x[rowInd, colInd]
	labRow <- NULL
	if (is.null(labRow))
	{
		if (is.null(rownames(x)))
		{
			labRow <- (1:nr)[rowInd]
		}
		else
		{
			labRow <- rownames(x)
		}
	}
	else
	{
		labRow <- labRow[rowInd]
	}
	labCol <- NULL
	if (is.null(labCol))
	{
		if (is.null(colnames(x)))
		{
			labCol <- (1:nc)[colInd]
		}
		else
		{
			colnames(x)
		}
	}
	else
	{
		labCol <- labCol[colInd]
	}
	if (scale == "row")
	{
		x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
		sx <- apply(x, 1, sd, na.rm = na.rm)
		x <- sweep(x, 1, sx, "/")
	}
	else if (scale == "column")
	{
		x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
		sx <- apply(x, 2, sd, na.rm = na.rm)
		x <- sweep(x, 2, sx, "/")
	}
	lmat <- rbind(c(NA, 3), 2:1)
	lwid <- c(if (doRdend) 1 else 0.05, 4)
	lhei <- c((if (doCdend) 1 else 0.05) + (if (!is.null(main)) 0.2 else 0), 4)
	if (!missing(ColSideColors))
	{
		if (!is.matrix(ColSideColors))
		{
			stop("'ColSideColors' must be a matrix")
		}
		if (!is.character(ColSideColors) || dim(ColSideColors)[1] != nc)
		{
			stop("'ColSideColors' dim()[2] must be of length ncol(x)")
		}
		lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
		lhei <- c(lhei[1], SideColorHeight, lhei[2])
	}
	if (!missing(RowSideColors))
	{
		if (!is.matrix(RowSideColors))
		{
			stop("'RowSideColors' must be a matrix")
		}
		if (!is.character(RowSideColors) || dim(RowSideColors)[1] != nr)
		{
			stop("'RowSideColors' must be a character vector of length nrow(x)")
		}
		lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[, 2] + 1)
		lwid <- c(lwid[1], SideColorHeight, lwid[2])
	}
	lmat[is.na(lmat)] <- 0
	if (verbose)
	{
		cat("layout: widths = ", lwid, ", heights = ", lhei, "; lmat=\n")
		print(lmat)
	}
	op <- par(no.readonly = TRUE)
	on.exit(par(op))
	layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
	if (!missing(RowSideColors))
	{
		par(mar = c(margins[1], 0, 0, 0.5))
		rsc = RowSideColors[rowInd, ]
		rsc.colors = matrix()
		rsc.names = names(table(rsc))
		rsc.i = 1
		for (rsc.name in rsc.names)
		{
			rsc.colors[rsc.i] = rsc.name
			rsc[rsc == rsc.name] = rsc.i
			rsc.i = rsc.i + 1
		}
		rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
		image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
		if (length(colnames(RowSideColors)) > 0)
		{
			axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), colnames(RowSideColors), las = 2, tick = FALSE)
		}
	}
	if (!missing(ColSideColors))
	{
		par(mar = c(0.5, 0, 0, margins[2]))
		csc = ColSideColors[colInd, ,drop=FALSE]
		csc.colors = matrix()
		csc.names = names(table(csc))
		csc.i = 1
		for (csc.name in csc.names)
		{
			csc.colors[csc.i] = csc.name
			csc[csc == csc.name] = csc.i
			csc.i = csc.i + 1
		}
		csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
		image(1:nrow(csc), 1:ncol(csc), csc, col = as.vector(csc.colors), axes = FALSE)
		if (!is.null(Colsep))  abline(v=(Colsep+0.5), col=sepColor, lwd=sepLWD)
		if (length(colnames(ColSideColors)) > 0)
		{
			axis(4, 1:ncol(csc), colnames(ColSideColors), las = 2, tick = FALSE)
		}
	}
	par(mar = c(margins[1], 0, 0, margins[2]))
	if (!symm || scale != "none")
	{
		x <- t(x)
	}
	if (revC)
	{
		iy <- nr:1
		ddr <- rev(ddr)
		x <- x[, iy]
	}
	else
	{
		iy <- 1:nr
	}
	image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "")
	if (!is.null(Colsep))
	{
		abline(v=Colsep+0.5, col=sepColor, lwd=sepLWD)
	}
	axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, cex.axis = cexCol)
	if (!is.null(xlab))
	{
		mtext(xlab, side = 1, line = margins[1] - 1.25)
	}
	axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, cex.axis = cexRow)
	if (!is.null(ylab))
	{
		mtext(ylab, side = 4, line = margins[2] - 1.25)
	}
	if (!missing(add.expr))
	{
		eval(substitute(add.expr))
	}
	par(mar = c(margins[1], 0, 0, 0))
	if (doRdend)
	{
		plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
	}
	else
	{
		frame()
	}
	par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2]))
	if (doCdend)
	{
		plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
	}
	else if (!is.null(main))
	{
		frame()
	}
	if (!is.null(main))
	{
		title(breakIntoTitle(main), cex.main = 1.5 * op[["cex.main"]])
	}
	invisible(list(rowInd = rowInd,
								 colInd = colInd,
								 Rowv = if (keep.dendro && doRdend) ddr,
								 Colv = if (keep.dendro && doCdend) ddc))
}

###########################################################################################
###########################################################################################

scale.trunc<-function(x, di=c("all", "row","col"), scal=FALSE, tp=0.8, ...)
{
	x<-as.matrix(x)
	nc<-ncol(x)
	nr<-nrow(x)
	if(scal)
	{
		dir<-match.arg(di)
		if(dir=="row")
		{
			temp<-t(scale(t(x),...))
		}
		else if(dir=="col")
		{
			temp<-scale(x,...)
		}
		else if(dir=="all")
		{
			temp<-scale(as.vector(x),...)
			temp<-matrix(temp, nr,nc, dimnames=list(rownames(x), colnames(x)))
		}
	}
	else
	{
		temp<-x
	}
	# added na.rm=TRUE to prevent problems
	bound<-quantile(temp, c((1-tp)/2, 1-(1-tp)/2), na.rm=TRUE)
	temp[temp<bound[1]]<-bound[1]
	temp[temp>bound[2]]<-bound[2]
	return(temp)
}

###########################################################################################
###########################################################################################

match.cols<-function(si)
{
	si<-data.frame(si, check.names=FALSE)
	for(i in 1:ncol(si))
	{
		si[,i]<-factor(si[,i])
	}
	nl<-sapply(si, nlevels)
	ref.ind<-which(nl==max(nl))[1]
	Lev<-levels(si[,ref.ind])
	for(i in 1:ncol(si))
	{
		ori.lev<-levels(si[,i])
		temp.tab<-table(si[,ref.ind], si[,i])
		indx<-unlist(sapply(Lev, function(x) which(temp.tab[x,]==max(temp.tab[x,]))[1]))
		if(any(duplicated(indx)))
		{
			dup<-unique(indx[duplicated(indx)])
			for(ii in dup)
			{
				keep<-which(temp.tab[,ii]==max(temp.tab[,ii]))[1]
				indx[indx==ii]<-NA
				indx[keep]<-ii
			}
		}
		new.lev<-ori.lev[indx]
		nfill<-sum(is.na(new.lev))
		not.in<-ori.lev[!ori.lev %in% new.lev]
		if(nfill>0)
		{
			new.lev[is.na(new.lev)]<-paste('NAfill', 1:nfill, sep="")
		}
		new.lev<-c(new.lev, not.in)
		si[,i]<-factor(as.vector(si[,i]), levels=new.lev)
	}
	return(si)
}

###########################################################################################
###########################################################################################

supervisedClustLegend<-function(theBatchIdColors, theBatchIdsForSamples, theSortedListOfBatchIds, theTitle, theFilename)
{
	version <- paste("MBatch", packageDescription("MBatch")["Version"], sep=" ")
	sortedBatch <- getListOfBatchNamesWithCounts(theBatchIdsForSamples, theSortedListOfBatchIds)
	myColors <- theBatchIdColors
	if (length(theBatchIdColors)>length(sortedBatch))
	{
		myColors <- theBatchIdColors[1:length(sortedBatch)]
	}
	mbatchStandardLegend(theTitle, version, sortedBatch, as.vector(myColors), NULL, theFilename)
	#Cairo PNG(filename=theFilename, width = 800, height = 1200, pointsize=24, bg = "transparent")
	#on.exit(dev.off(), add = TRUE)
	#plot.new()
	#logDebug("writePcaAnalysis_Legend -- theSortedListOfBatchIds = ", paste(theSortedListOfBatchIds, collapse=","))
	###logDebug("writePcaAnalysis_Legend -- sortedBatch = ", paste(sortedBatch, collapse=","))
	#leg end("top", legend=sortedBatch, title=theTitle, title.col="black",
	#			 text.col=as.vector(theBatchIdColors), col=as.vector(theBatchIdColors),
	#			 lty=1, cex=0.70)
}

###########################################################################################
###########################################################################################

openAndWriteIssuesLogFileSC<-function(theOutputDir)
{
  myFile <- file(cleanFilePath(theOutputDir, "error.log"), "w+")
  on.exit(close(myFile))
  cat("Clustering not possible\n", file=myFile, append=TRUE)
}

writeSCmatrix<-function(theMatrix, theOutdir, theMatrixFile, theTitle, theTitleFile)
{
  writeAsGenericMatrix(cleanFilePath(theOutdir, theMatrixFile), theMatrix)
  myFile <- file(cleanFilePath(theOutdir, theTitleFile), "w+")
  on.exit(close(myFile))
  cat(theTitle, file=myFile)
}

###########################################################################################
###########################################################################################
