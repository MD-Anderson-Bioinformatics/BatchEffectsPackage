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

createBatchEffectsOutput_TRINOVA<-function(theMatrixGeneData, theDataframeBatchData,
																					theTitle, theOutputPath)
{
  textFile <- file.path(theOutputPath, "TRINOVA.tsv")
  pngFile <- file.path(theOutputPath, "TRINOVA.png")

	mySuccess <- FALSE
	myError <- NULL
	tryCatch({
		mySuccess <- calcTRINOVA(theMatrixGeneData, theDataframeBatchData,
		                        colnames(theDataframeBatchData)[2:4],
		                        textFile, pngFile, theTitle)
	},warn=function(myError)
	{
	  logError("createBatchEffectsOutput_TRINOVA(): generating TRINOVA results, ERROR= ", myError)
	  openAndWriteIssuesLogFile(theOutputPath)
	},error=function(myError)
	{
	  logError("createBatchEffectsOutput_TRINOVA(): generating TRINOVA results, ERROR= ", myError)
	  openAndWriteIssuesLogFile(theOutputPath)
	})

	if (FALSE==mySuccess)
	{
		#clean up directory
	  if (file.exists(textFile))
	  {
	    file.remove(textFile)
	  }
	  if (file.exists(pngFile))
	  {
	    file.remove(pngFile)
	  }
	}
	return(mySuccess)
}


openAndWriteIssuesLogFile<-function(theOutputDir)
{
  myFile <- file(file.path(theOutputDir, "error.log"), "w+")
  on.exit(close(myFile))
  cat("Not calculated\n", file=myFile, append=TRUE)
}

calcTRINOVA<-function(dat, si, By, theTextFile, thePngFile, theTitle)
{
  logDebug("calcTRINOVA")
  logDebug("theTextFile=", theTextFile)
  logDebug("thePngFile=", thePngFile)
	logDebug("calcTRINOVA - length(By)==", length(By))
	stopifnotWithLogging("Need at least three entries by By", length(By)==3)
	logDebug("calcTRINOVA - all(colnames(dat)==si$Sample)")
	stopifnotWithLogging("Samples in data and batches should match and be in same order", all(colnames(dat)==si$Sample))
	logDebug("calcTRINOVA - all(By %in% colnames(si))")
	stopifnotWithLogging("By should exist in batch type list", all(By %in% colnames(si)))
	logDebug("calcTRINOVA - input.si")
	input.si<-si[,By]
	trinova.result<-NULL

	logDebug("doTriNova - input.si")
	trinova.result<-doTriNova(dat, input.si)
	setNames <- colnames(si)[2:4]
	residual <- trinova.result$partition[["Residual"]]
	b1 <- trinova.result$partition[[setNames[1]]]
	b2 <- trinova.result$partition[[setNames[2]]]
	b3 <- trinova.result$partition[[setNames[3]]]
	b12 <- trinova.result$partition[[paste(setNames[1], setNames[2], sep=".")]]
	b13 <- trinova.result$partition[[paste(setNames[1], setNames[3], sep=".")]]
	b23 <- trinova.result$partition[[paste(setNames[2], setNames[3], sep=".")]]
	b123 <- trinova.result$partition[[paste(setNames[1], setNames[2], setNames[3], sep=".")]]

	file = writeTrinovaTSV(theTextFile, theTitle, setNames, b1, b2, b3, b12, b13, b23, b123, residual)
	writeTrinovaPNG(thePngFile, theTitle, setNames, b1, b2, b3, b12, b13, b23, b123, residual)
	return(file)
}

writeTrinovaPNG <- function(thePngFile, theTitle, theSetNames, b1, b2, b3, b12, b13, b23, b123, residual)
{
  CairoPNG(filename = thePngFile, width = 1000, height = 1000, pointsize=24)
  on.exit(dev.off(), add = TRUE)
  vennLike <- euler(
    c(
      "A" = round(residual, 3),
      "A&B" = round(b1, 3),
      "A&C" = round(b2, 3),
      "A&D" = round(b3, 3),
      "A&B&C" = round(b12, 3),
      "A&B&D" = round(b13, 3),
      "A&C&D" = round(b23, 3),
      "A&B&C&D" = round(b123, 3)
    )
  )
  myLabels = c("Residual", theSetNames[1], theSetNames[2], theSetNames[3])
  # not in docs, but works like ggplot/lattice, needs print - plot
  print(
    plot(vennLike,
         quantities = list(type = c("percent")),
         labels=myLabels, main=theTitle, adjust_labels = TRUE))
}

writeTrinovaTSV <- function(theTextFile, theTitle, theSetNames, b1, b2, b3, b12, b13, b23, b123, residual)
{
  df = data.frame(keys = c("Title" , "Residual",
                           theSetNames[1],
                           theSetNames[2],
                           theSetNames[3],
                           paste(theSetNames[1], theSetNames[2], sep=" - "),
                           paste(theSetNames[1], theSetNames[3], sep=" - "),
                           paste(theSetNames[2], theSetNames[3], sep=" - "),
                           paste(theSetNames[1], theSetNames[2], theSetNames[3], sep=" - ")),
                  values = c(theTitle, round(residual, 3),
                             round(b1, 3), round(b2, 3), round(b3, 3),
                             round(b12, 3), round(b13, 3), round(b23, 3),
                             round(b123, 3)))
  writeAsDataframe(theTextFile, df)
  theTextFile
}

## Function doTriNova to perform PCA and ANOVA. Author: Nianxiang Zhang Date: 2011-06-22
doTriNova<-function(dat, si, usecor=FALSE, center=TRUE)
{
  warnLevel<-getOption("warn")
  on.exit(options(warn=warnLevel))
  options(warn=-1) # disable warnings
  logDebug("doTriNova - disable warnings")
  cn=colnames(si)
	stopifnotWithLogging("Number of sample should match in data and batches", ncol(dat)==nrow(si))
	stopifnotWithLogging("Need three entries in batch types", length(cn)==3)
	pca<-SamplePCA(dat[!is.na(rowSums(dat)),], usecor=usecor, center=center)
	pc<-pca@scores
	Var.pct<-pca@variances/sum(pca@variances)


	all.anova<-apply(pc, 2, Anova3, si=si)
	adj.trinova<-sweep(all.anova,2,Var.pct, '*')
	results<-rowSums(adj.trinova)
	#results1<-rbind(sum.square=s, percent=pct)
	names(results)<-c(cn, paste(cn[-3], collapse='.'),paste(cn[-2], collapse='.'), paste(cn[-1], collapse='.'), paste(cn, collapse='.'), 'Residual')
	Names<-cn
	names(Names)<-c('V1', 'V2', 'V3')
	return(list(partition=results, name=Names))
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
