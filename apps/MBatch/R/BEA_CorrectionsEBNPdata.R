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

# A class to examine the statistical characteristics of and between two dataframes or matrixs
# the section and subsection comments are also used by a python program which can document the codes in pdf with latex directly

# define the class
# The input of this class are two data frames
# The first column of data frame should be entry name
# Samples arranged in the column
setClass("EBNPdata", representation(
	mData1="matrix",
	mData2="matrix",
	mat1Com="matrix",
	mat2Com="matrix",
	pairCorMat = "data.frame",
	randCorMat = "data.frame"
	))

setGeneric(name="makeCommonRows", def=function(Object, ...) {standardGeneric("makeCommonRows")})
setGeneric(name="makeCommonCols", def=function(Object, ...) {standardGeneric("makeCommonCols")})
#setGeneric(name="evaluateMatch", def=function(Object, ...) {standardGeneric("evaluateMatch")})
setGeneric(name="asSameOrder", def=function(Object, ...) {standardGeneric("asSameOrder")})
#setGeneric(name="scatPlotOnSelectCols", def=function(Object, ...) {standardGeneric("scatPlotOnSelectCols")})
setGeneric(name="corPairIn2Set", def=function(Object, ...) {standardGeneric("corPairIn2Set")})
setGeneric(name="corUnpairIn2Set", def=function(Object, ...) {standardGeneric("corUnpairIn2Set")})
setGeneric(name="plotPairUnpairCor", def=function(Object, ...) {standardGeneric("plotPairUnpairCor")})
setGeneric(name="getBiComOrder", def=function(Object, ...) {standardGeneric("getBiComOrder")})

# Initialize function of the class
setMethod("initialize", "EBNPdata", function(.Object, theData1, theData2, theColNameConvertFun, theColNameDuplicateFun, theRowNameConvertFun, theRowNameDuplicateFun)
{
	.Object@mData1 <- theData1
	.Object@mData2 <- theData2
	logDebug(paste("dim(mData1)=", paste(dim(.Object@mData1), collapse=" ", sep=" "), sep=""))
	logDebug(paste("dim(mData2)=", paste(dim(.Object@mData2), collapse=" ", sep=" "), sep=""))
	if (FALSE==is.null(theColNameConvertFun))
	{
		logDebug("converting column names")
		colnames(.Object@mData1) <- theColNameConvertFun(colnames(.Object@mData1))
		colnames(.Object@mData2) <- theColNameConvertFun(colnames(.Object@mData2))
		logDebug(paste("dim(mData1)=", paste(dim(.Object@mData1), collapse=" ", sep=" "), sep=""))
		logDebug(paste("dim(mData2)=", paste(dim(.Object@mData2), collapse=" ", sep=" "), sep=""))
	}
	if (FALSE==is.null(theColNameDuplicateFun))
	{
		logDebug("removing column duplicates")
		.Object@mData1 <- theColNameDuplicateFun(.Object@mData1)
		.Object@mData2 <- theColNameDuplicateFun(.Object@mData2)
		logDebug(paste("dim(mData1)=", paste(dim(.Object@mData1), collapse=" ", sep=" "), sep=""))
		logDebug(paste("dim(mData2)=", paste(dim(.Object@mData2), collapse=" ", sep=" "), sep=""))
	}
	if (FALSE==is.null(theRowNameConvertFun))
	{
		logDebug("converting row names")
		rownames(.Object@mData1) <- theRowNameConvertFun(rownames(.Object@mData1))
		rownames(.Object@mData2) <- theRowNameConvertFun(rownames(.Object@mData2))
		logDebug(paste("dim(mData1)=", paste(dim(.Object@mData1), collapse=" ", sep=" "), sep=""))
		logDebug(paste("dim(mData2)=", paste(dim(.Object@mData2), collapse=" ", sep=" "), sep=""))
	}
	if (FALSE==is.null(theRowNameDuplicateFun))
	{
		logDebug("removing row duplicates")
		.Object@mData1 <- theRowNameDuplicateFun(.Object@mData1)
		.Object@mData2 <- theRowNameDuplicateFun(.Object@mData2)
		logDebug(paste("dim(mData1)=", paste(dim(.Object@mData1), collapse=" ", sep=" "), sep=""))
		logDebug(paste("dim(mData2)=", paste(dim(.Object@mData2), collapse=" ", sep=" "), sep=""))
	}
	.Object
})

# Evaluate whether there is enough matched items in two sets
# Evaluate how many rows and columns of these two dataframes match
# This method is used to remind the users to format the entry names in column and row automatically at the beginning of the study
# When there is more than 20% of entry name of the smaller data frame that doesn't match to another data frame, this method will give the warning
#setMethod("evaluateMatch", "EBNPdata", function(Object, byColnames= TRUE, byCol) {
#  df1  <- Object@mData1
#  df2  <- Object@mData2
#  if(byColnames==TRUE) {
#    if( ncol(df2) > ncol(df1)) {
#      match.i <- which(colnames(df1) %in% colnames(df2))
#      matchRatio  <- length(match.i)/ncol(df1)
#      if(matchRatio < 0.8) {
#        warning("there is more than 20% columns of data1 which are not in data2")
#      }
#    } else {
#      match.i <- which(colnames(df2) %in% colnames(df1))
#      if(length(match.i)/ncol(df2) < 0.8) {
#        warning("there is more than 20% columns of data2 which are not in data1")
#      }
#    }
#  }
#
#  if(byCol=="row.names") {
#    rownamesdf1  <- row.names(df1)
#    rownamesdf2    <- row.names(df2)
#  } else {
#    rownamesdf1  <- df1[ ,byCol]
#    rownamesdf2  <- df2[ ,byCol]
#  }
#
#  if (nrow(df2) > nrow(df1))  {
#    match.i <- which(rownamesdf1 %in% rownamesdf2)
#    if(length(match.i)/nrow(df1) < 0.8) {
#      warning("there are more than 20% rows of data1 which are not in data2")
#    }
#  } else {
#    match.i <- which(rownamesdf2 %in% rownamesdf1)
#    if(length(match.i)/nrow(df2) < 0.8) {
#      warning("there are more than 20% rows of data2 which are not in data1")
#    }
#  }
#})

# this method arrange the data frame in same order by entry names
# The two data frame should have same items
setMethod("asSameOrder","EBNPdata", function(Object, by="both") {
  mat1  <- Object@mData1
  mat2  <- Object@mData2
	logDebug("asSameOrder before")
	printMatrix(mat1)
	printMatrix(mat2)
  # first check the entry name is completely match
  if( !(by %in% c("row","col","both"))) {
    message ("please select arrange order by row or col or both")
  }
	logDebug("asSameOrder after if 1")
	printMatrix(mat1)
	printMatrix(mat2)
	if( (by =="col") | (by=="both")) {
    m.i <- match(colnames(mat1),colnames(mat2))
    mat2 <- mat2[ ,m.i]
  }
	logDebug("asSameOrder after if 2")
	printMatrix(mat1)
	printMatrix(mat2)
	if( (by =="row") | (by=="both"))
	{
		logDebug("row.names(mat1)")
		printElements(head(row.names(mat1), n=100))
		logDebug("rownames(mat1)")
		printElements(head(rownames(mat1), n=100))
		logDebug("row.names(mat2)")
		printElements(head(row.names(mat2), n=100))
		logDebug("rownames(mat2)")
		printElements(head(rownames(mat2), n=100))
	  m.i <- match(row.names(mat1),row.names(mat2))
		logDebug("m.i")
		printElements(head(m.i, n=100))
	  mat2 <- mat2[m.i, ]
  }
	logDebug("asSameOrder after if 3")
	printMatrix(mat1)
	printMatrix(mat2)
	Object@mData2  <- mat2
  return(Object)
})

#setMethod("scatPlotOnSelectCols" , "EBNPdata", function(Object, items, ...){
#  mat1  <- Object@mData1
#  mat2  <- Object@mData2
#  numItems <- length(items)
#  colorList <- rainbow(numList)
#  xRange <- c(min(mat1[ ,items]), max(mat1[ ,items]))
#  yRange <- c(min(mat2[ ,items]), max(mat2[ ,items]))
#  for ( i in 1: numList ) {
#    pair <- merge(mat1[ ,items[i]], mat2[ ,items[i]], by.x="row.names", by.y="row.names")
#    if ( i == 1 ) {
#      # when there is only one column, ??
#      reg1 <- lm(pair[ ,2]~ pair[ ,1])
#      plot(pair[,1],pair[,2], col= colorList[i], xlim= xRange, ylim= yRange, ...)
#      abline(reg1, col=colorList[i])
#    } else {
#      par(new = TRUE)
#      reg1 <- lm(pair[ ,3]~ pair[ ,2])
#      plot(pair[,2],pair[,3], col= colorList[i], xlim= xRange, ylim= yRange, ann=FALSE)
#      abline(reg1, col=colorList[i])
#    }
#  }
#})


# Determine the correlation of paired samples in two datasets
# the IDs/names for compare will arranged on row
# the IDs/names should be in the first column
# require the IDs/names of the row are completely same, and require the columns are the same and in the same order
# The reason to use loop is because there may be one to multiple
setMethod("corPairIn2Set" , "EBNPdata", function(Object, method="pearson"){
  #mat1  <-  moveColumn2Rowname(Object@mData1)
  #mat2  <-  moveColumn2Rowname(Object@mData2)
  mat1 <- data.matrix(Object@mData1)
  mat2 <- data.matrix(Object@mData2)
  numSamples <- nrow(mat1)
  pairCorMat <- data.frame(IDs=character(),R=numeric(), pvalue=numeric(),meanOfmat1=numeric(),stdOfmat1=numeric(),meanOfmat2=numeric(),stdOfmat2=numeric() , stringsAsFactors=FALSE)

  # with loop, can handle the situation with mutilple replicates
  for ( i in 1: numSamples) {
    ID <- row.names(mat1)[i]
    match.i <- which(row.names(mat2) == ID)

    if(length(match.i)==0) next

    if(length(match.i)==1)
    {
      vec1 <- mat1[i, ]
      vec2 <- mat2[match.i, ]
      cor <- cor.test(vec1,vec2,method=method)
      pRmat <- data.frame(samples=ID, R= cor[[4]], pvalue= cor[[3]] , meanOfmat1=mean(vec1),stdOfmat1=sd(vec1), meanOfmat2=mean(vec2),stdOfmat2=sd(vec2), stringsAsFactors=FALSE)
    } else {
      data1 <- mat1[i, ]
      data2 <- mat2[match.i, ]

      corList <- apply(data2,1, function(x) cor.test(data1,x, method=method))
      pR <- lapply(corList, function(x) c(x[[4]],x[[3]]))
      pRmat <- as.data.frame(matrix(unlist(pR),ncol=2,byrow=TRUE))
      pRmat <- cbind(ID, pRmat, mean(data1),sd(data1), apply(data2,1,mean),apply(data2,1,sd))
      colnames(pRmat)<- c("IDs", "R","pvalue", "meanOfset1", "stdOfset1", "meanOfset2", "stdOfset2")
    }
    pairCorMat<- rbind(pairCorMat, pRmat)
  }

  if( method == "spearman" ) {
    colnames(pairCorMat)[2]="rho"
  }
  Object@pairCorMat  <- pairCorMat
  return(Object)
})

# Determine random correlation using non-pairing samples in two matrixs
setMethod("corUnpairIn2Set" , "EBNPdata", function(Object, method="pearson"){
  #mat1  <-  moveColumn2Rowname(Object@mData1)
  #mat2  <-  moveColumn2Rowname(Object@mData2)
  mat1 <- as.matrix(Object@mData1)
  mat2 <- as.matrix(Object@mData2)
  numSamples <- nrow(mat1)
  randCorMat <- data.frame(setmat1ID=character(),R=numeric(0),pvalue=numeric(0), setmat2ID=character(), stringsAsFactors=FALSE)

  for ( i in 1: numSamples) {
    ID <- row.names(mat1)[i]
    non.i <- which(rownames(mat2)!=ID)
    data1 <- mat1[i,]
    data2 <- mat2[non.i,]
    corList <- apply(data2,1, function(x) cor.test(data1,x, method=method))
    pR <- lapply(corList, function(x) c(x[[4]],x[[3]]))
    pRmat <- as.data.frame(matrix(unlist(pR),ncol=2,byrow=TRUE))
    pRmat <- cbind(ID, pRmat, row.names(mat2)[non.i])
    colnames(pRmat)<- c("setmat1ID", "R", "pvalue", "setmat2ID")
    randCorMat<- rbind(randCorMat, pRmat)
  }

  if( method == "spearman" ) {
    colnames(randCorMat)[2]="rho"
  }
  Object@randCorMat  <- randCorMat
  return(Object)
})

# Plot the histogram and density of the correlations between pair and unpair items in two sets
setMethod("plotPairUnpairCor" , "EBNPdata", function(Object, histogram=FALSE, legendPos="topleft",legend=NULL, xRange = NULL, breaks=NA, ...){
  vec1  <- Object@pairCorMat[,2]
  vec2  <- Object@randCorMat[ ,2]

  if ( histogram == TRUE ){
  	par(mfrow=c(2,1))
  	xRange <- c(min(c(vec1,vec2)), max(c(vec1,vec2)))
  	hist(vec1, ann=FALSE, col="skyblue2", xlim=xRange, breaks=breaks,...)
  	par(new=TRUE)
  	hist(vec2, xlim = xRange, axes = FALSE, col="green", breaks=breaks,...)
  	legend(legendPos, legend, lty=1, col=c("skyblue2","green"), bty='n')
  }

  x <- density(vec1)
  y <- density(vec2)
  if ( is.null(xRange)) {
    xRange <- c(min(c(vec1,vec2)), max(c(vec1,vec2)))
  }
  plot(x, type='l',xlim=xRange,ann=FALSE, col="deepskyblue", ...)
  par(new =TRUE)
  plot(y, type='l',col="firebrick2",xlim=xRange, axes = FALSE, ...)
  legend(legendPos, legend, lty=1, col=c("deepskyblue","firebrick2"), bty='n')
})

# Get common elements in two matrix
# get a new EBNPdata object with only common elements in two matrixs
setMethod("makeCommonRows" , "EBNPdata", function(Object)
{
	mat1  <- Object@mData1
	mat2  <- Object@mData2
	xCol <- rownames(mat1)
	yCol <- rownames(mat2)
	rowCom = intersect(xCol, yCol)
	com.i = match(rowCom, xCol)
	mat1com <- mat1[com.i, ]
	com.i <- match( rowCom, yCol)
	mat2com <- mat2[com.i, ]
	return(new("EBNPdata",mat1com, mat2com, NULL, NULL, NULL, NULL))
})

# Get common elements in two matrix
# get a new EBNPdata object with only common elements in two matrixs
setMethod("makeCommonCols" , "EBNPdata", function(Object)
{
	mat1  <- Object@mData1
	mat2  <- Object@mData2
	colnameCom  = intersect(colnames(mat1), colnames(mat2))
	mat1com = mat1[,colnameCom]
	mat2com = mat2[,colnameCom]
	return(new("EBNPdata",mat1com, mat2com, NULL, NULL, NULL, NULL))
})


# get common parts of two DF in both column and row and match the order of the common parts
setMethod( "getBiComOrder", "EBNPdata", function(Object)
{
	# Remove uncommon rows in the two datasets
	logDebug("makeCommonRows")
	dataCom <- makeCommonRows(Object)
	# Remove uncommons column in two dataframe
	logDebug("makeCommonCols")
	dataCom <- makeCommonCols(dataCom)
	logDebug("asSameOrder")
	dataCom <- asSameOrder(dataCom)
	Object@mat1Com  <- dataCom@mData1
	Object@mat2Com  <- dataCom@mData2
	Object
})
