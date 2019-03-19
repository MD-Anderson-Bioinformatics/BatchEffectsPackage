#MBatch Copyright ? 2011, 2012, 2013, 2014, 2015, 2016, 2017 University of Texas MD Anderson Cancer Center
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(Vennerable, warn.conflicts=FALSE, verbose=FALSE)

createBatchEffectsOutput_MANOVA<-function(theMatrixGeneData, theDataframeBatchData,
																					theTitle, theOutputPath)
{
	#filename <- file.path(theOutputPath, "MANOVA.svg")
	filename <- file.path(theOutputPath, "MANOVA_diagram.png")

	mySuccess <- FALSE
	myError <- NULL
	tryCatch({
		#if VennMasterMANOVA() completes, it returns TRUE; otherwise any error occurs, it returns FALSE.
		mySuccess <- VennMasterMANOVA(theMatrixGeneData, theDataframeBatchData, colnames(theDataframeBatchData)[2:4], filename, theOutputPath)
	},error=function(myError)
	{
		logError("createBatchEffectsOutput_MANOVA(): generating MANOVA results, ERROR= ", myError)
	})

	if (FALSE==mySuccess)
	{
		#clean up directory
		file.remove(filename)
		CairoPNG(filename=filename, width = 1000, height = 1000, pointsize=14)
		on.exit(dev.off(), add = TRUE)
		boxplot(c(0,0),
						cex=0.2, axes=FALSE,
						main = paste("Unable to calculate--problem processing MANOVA ", theTitle, " \nsee log file.",sep=""))

	}
}

makeMANOVAFileName_PNG<-function(theDir, theDiagramOrLegend, theLegendType="")
{
	stopifnotWithLogging("theDiagramOrLegend should be Diagram|Legend", (("Diagram"==theDiagramOrLegend)||("Legend"==theDiagramOrLegend)))
	if ("Legend"!=theDiagramOrLegend)
	{
		if (""!=theLegendType)
		{
			stopWithLogging(paste("theLegendType is ", theLegendType, " which must be empty for theDiagramOrLegend=", theDiagramOrLegend, sep=""))
		}
	}
	if ("Legend"==theDiagramOrLegend)
	{
		createDirPlusFilename(theDir, "MANOVA_", theDiagramOrLegend, "-", theLegendType, ".png")
	}
	else
	{
		createDirPlusFilename(theDir, "MANOVA_", theDiagramOrLegend, ".png")
	}
}

VennMasterMANOVA<-function(dat, si, By, theFilename, theOutputPath)
{
	logDebug("VennMasterMANOVA")
	#vennpath <- system.file("venn", "venn.jar", package="MBatch")
	#Xmx <-"1024M"
	#prepare the data for the pipeline
	logDebug("VennMasterMANOVA - length(By)==3")
	stopifnotWithLogging("Need at least three entries by By", length(By)==3)
	logDebug("VennMasterMANOVA - all(colnames(dat)==si$Sample)")
	stopifnotWithLogging("Samples in data and batches should match and be in same order", all(colnames(dat)==si$Sample))
	logDebug("VennMasterMANOVA - all(By %in% colnames(si))")
	stopifnotWithLogging("By should exist in batch type list", all(By %in% colnames(si)))
	logDebug("VennMasterMANOVA - input.si")
	input.si<-si[,By]
	manova.result<-NULL

	logDebug("doTriManova - input.si")
	manova.result<-doTriManova(dat, input.si)
	#venn.list<-gen.venn.dat(manova.result)
	#simfile <- file.path(theOutputPath, ".simfile.tmp")
	#infile <- file.path(theOutputPath,".vennlistin.tmp")
	#svg <- theFilename
	#write.table(venn.list, infile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
	CairoPNG(filename=theFilename, width = 1000, height = 1000, pointsize=14, bg = "transparent")
	on.exit(dev.off(), add = TRUE)
	setNames <- colnames(si)[2:4]
	residual <- manova.result$partition[["Residual"]]
	b1 <- manova.result$partition[[setNames[1]]]
	b2 <- manova.result$partition[[setNames[2]]]
	b12 <- manova.result$partition[[paste(setNames[1], setNames[2], sep=".")]]
	b3 <- manova.result$partition[[setNames[3]]]
	b13 <- manova.result$partition[[paste(setNames[1], setNames[3], sep=".")]]
	b23 <- manova.result$partition[[paste(setNames[2], setNames[3], sep=".")]]
	b123 <- manova.result$partition[[paste(setNames[1], setNames[2], setNames[3], sep=".")]]
	VenTest <- Venn(SetNames=setNames, Weight=c(residual, b1, b2, b12, b3, b13, b23, b123))
	plot(VenTest, doWeights=FALSE)

	# plot and get midpoints
	##mp <- barplot(manova.result$partition)
	# use that to add actual values
	##text(mp, manova.result$partition, labels=formatC(manova.result$partition, format="f", digits=4), pos=1)
	#text(mp, manova.result$partition, labels=manova.result$partition, pos=1)
	## call vennmaster cmdline mode
	#system(paste("java -Xmx", Xmx, " -jar ", vennpath, " --list ", infile, " --svg ", svg, " --sim ", simfile, sep=""))
	#cost <- read.table(simfile, header=FALSE, skip=1, sep="\t")
	## delete tmp input files
	#file.remove(infile)
	#file.remove(simfile)
	#stop("VennMasterMANOVA(): any error happens in this function.")
	return(TRUE)
}

## Function doTriManova to perform PCA and ANOVA. Author: Nianxiang Zhang Date: 2011-06-22
doTriManova<-function(dat, si, usecor=FALSE, center=TRUE)
{
	cn=colnames(si)
	stopifnotWithLogging("Number of sample should match in data and batches", ncol(dat)==nrow(si))
	stopifnotWithLogging("Need three entries in batch types", length(cn)==3)
	pca<-SamplePCA(dat[!is.na(rowSums(dat)),], usecor=usecor, center=center)
	pc<-pca@scores
	Var.pct<-pca@variances/sum(pca@variances)


	all.anova<-apply(pc, 2, Anova3, si=si)
	adj.manova<-sweep(all.anova,2,Var.pct, '*')
	results<-rowSums(adj.manova)
	#results1<-rbind(sum.square=s, percent=pct)
	names(results)<-c(cn, paste(cn[-3], collapse='.'),paste(cn[-2], collapse='.'), paste(cn[-1], collapse='.'), paste(cn, collapse='.'), 'Residual')
	Names<-cn
	names(Names)<-c('V1', 'V2', 'V3')
	return(list(partition=results, name=Names))
}

## Function gen.venn.dat to prepare a matrix to be used by VennMaster Author: Nianxiang Zhang Date: 2011-08-22
gen.venn.dat<-function(manova.result, len=2000)
{
	vardat<-manova.result[[1]]

	Name<-manova.result[[2]]
	stopifnotWithLogging("Need numeric results of length 8", is.numeric(vardat) & length(vardat)==8)
	stopifnotWithLogging("Need three results from manova call", length(Name)==3)
	counts<-round(vardat*len, 0)
	out<-matrix(character(0), 0, 2)
	for(i in 1:7)
	{
		nam<-paste('V',i, sep='')
		mati<-switch( i ,
									'1'=cbind(rep(nam, counts[1]), rep(Name[1], counts[1])),
									'2'=cbind(rep(nam, counts[2]), rep(Name[2], counts[2])),
									'3'=cbind(rep(nam, counts[3]), rep(Name[3], counts[3])),
									'4'=cbind(rep(nam, 2*counts[i]), rep(c(Name[1],Name[2]), each=counts[4])),
									'5'=cbind(rep(nam, 2*counts[i]), rep(c(Name[1],Name[3]), each=counts[5])),
									'6'=cbind(rep(nam, 2*counts[i]), rep(c(Name[2],Name[3]), each=counts[6])),
									'7'=cbind(rep(nam, 3*counts[i]), rep(c(Name[1],Name[2],Name[3]), each=counts[7]))
		)
		out<-rbind(out, mati)
	}
	return(out)
}

Anova3<-function(vec, si)
{
	#with(COAD.rna.si, table(TSS, Batch))
	y<-vec
	x1<-factor(si[,1])
	x2<-factor(si[,2])
	x3<-factor(si[,3])
	x<-list(x1,x2,x3)
	warnLevel<-getOption("warn")
	on.exit(options(warn=warnLevel))
	options(warn=-1)
	sst<-stats::anova(stats::lm(y~1))[,2]
	temp<-stats::anova(lm(y~x1+x2+x3))
	z123<-temp['Residuals','Sum Sq']
	temp<-stats::anova(lm(y~x3+x2))
	z23<-temp['Residuals','Sum Sq']
	temp<-stats::anova(lm(y~x1+x2))
	z12<-temp['Residuals','Sum Sq']
	temp<-stats::anova(lm(y~x3+x1))
	z13<-temp['Residuals','Sum Sq']
	temp<-stats::anova(lm(y~x1))
	z1<-temp['Residuals','Sum Sq']
	temp<-stats::anova(lm(y~x2))
	z2<-temp['Residuals','Sum Sq']
	temp<-stats::anova(lm(y~x3))
	z3<-temp['Residuals','Sum Sq']
	options(warn=warnLevel)
	V1<-z23-z123
	V2<-z13-z123
	V3<-z12-z123
	V12<-z3-z123-V1-V2
	V12<-ifelse(V12<0,0,V12)
	V13<-z2-z123-V1-V3
	V13<-ifelse(V13<0,0,V13)
	V23<-z1-z123-V3-V2
	V23<-ifelse(V23<0,0,V23)
	V123<-sst-z123-V1-V2-V3-V12-V13-V23
	V123<-ifelse(V123<0,0,V123)
	Res<-sst-V1-V2-V3-V12-V13-V23-V123
	s<-c(V1, V2, V3, V12, V13, V23, V123, Res)
	pct<-s/sst
	return(pct)
}
