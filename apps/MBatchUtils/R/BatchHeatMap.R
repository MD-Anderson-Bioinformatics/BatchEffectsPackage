# MBatchUtils Copyright (c) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021 University of Texas MD Anderson Cancer Center
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
# MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>

buildBatchHeatMap_Files <- function(theMatrixFile, theBatchFile, theTitle, theOutputFile, theSortByType,
                                    theRowType="scholar", theColType="bio.tcga.barcode.sample",
                                    theRowCluster=NULL, theColCluster=NULL,
                                    theShaidyMapGen, theNgchmWidgetJs,
                                    theShaidyMapGenJava="/usr/bin/java",
                                    theShaidyMapGenArgs="-Xmx16G")
{
  message("buildBatchHeatMap_Files")
  message(theMatrixFile)
  matData <- readAsGenericMatrix(theMatrixFile)
  message(theBatchFile)
  dfData <- readAsGenericDataframe(theBatchFile)
  buildBatchHeatMap_Structures(matData, dfData, theTitle, theOutputFile, theSortByType,
                               theRowType, theColType,
                               theRowCluster, theColCluster,
                               theShaidyMapGen, theNgchmWidgetJs,
                               theShaidyMapGenJava,
                               theShaidyMapGenArgs)
}

compressIntoFilename<-function(theString)
{
  ### listing whole list of characters out looks wrong, but is locale independent
  theString <- gsub("[^ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789/\\]", "", theString)
  theString <- gsub("\\", "_", theString, fixed=TRUE)
  theString <- gsub("/", "_", theString, fixed=TRUE)
  return(theString)
}

buildBatchHeatMap_Structures <- function(theMatrixData, theBatchData, theTitle, theOutputFile, theSortByType,
                                         theRowType="scholar", theColType="bio.tcga.barcode.sample",
                                         theRowCluster=NULL, theColCluster=NULL,
                                         theShaidyMapGen,
                                         theNgchmWidgetJs,
                                         theShaidyMapGenJava="/usr/bin/java",
                                         theShaidyMapGenArgs="-Xmx16G",
                                         theEncoding="en_US.UTF-8")
{
  # need to do this, since NGCHM uses title as filename, so slashed become directories.
  theTitle <- compressIntoFilename(theTitle)
  message("buildBatchHeatMap_Structures")
  print(theTitle)
  print(theOutputFile)
  dir.create(dirname(theOutputFile), showWarnings=FALSE, recursive=TRUE)
  print(theSortByType)
  message("theMatrixData size ", dim(theMatrixData)[1], " ", dim(theMatrixData)[2])
  message("theBatchData size ", dim(theBatchData)[1], " ", dim(theBatchData)[2])
  message("sort batches")
  sortOrder <- order(theBatchData[theSortByType])
  theBatchData <- theBatchData[sortOrder,]
  theMatrixData <- theMatrixData[,sortOrder]
  message("compute row clusters")
  rowClusters <- rownames(theMatrixData)
  if (!is.null(theRowCluster))
  {
    # c("pearson", "ward.D2")
    rowClusters <- as.dendrogram(hierClustForNgchmCalc(t(theMatrixData), theRowCluster[1], theRowCluster[2]))
  }
  message("compute col clusters")
  colClusters <- colnames(theMatrixData)
  if (!is.null(theColCluster))
  {
    # c("pearson", "ward.D2")
    colClusters <- as.dendrogram(hierClustForNgchmCalc(theMatrixData, theColCluster[1], theColCluster[2]))
  }
  message("make colormapBreaks")
  centeredMatrixData <- bidirectionalCentering(theMatrixData)
  colormapBreaks1 <- makeGenericColorMap(centeredMatrixData)
  colormapBreaks2 <- makeGenericColorMap(theMatrixData)
  print(colormapBreaks1)
  print(colormapBreaks2)
  message("make color map")
  # High = red
  # Medium (usually corresponding to the median) = light gray
  # Low = blue
  # Missing = white
  cmap1 <- chmNewColorMap(colormapBreaks1, colors = c('blue','lightgrey', 'red'), missing.color='white');
  cmap2 <- chmNewColorMap(colormapBreaks2, colors = c('blue','lightgrey', 'red'), missing.color='white');
  message("make data layer")
  layer1 <- chmNewDataLayer('bi-directional median centered', centeredMatrixData, colors=cmap1, summarizationMethod = "average")
  layer2 <- chmNewDataLayer('original', theMatrixData, colors=cmap2, summarizationMethod = "average")
  message("make CHM")
  chm <- chmNew(name=theTitle, layer1, layer2, rowOrder=rowClusters, colOrder=colClusters, rowAxisType=theRowType, colAxisType=theColType)
  message("add caption")
  chm <- chmAddProperty(chm, "chm.info.caption", paste("MBatch NGCHM: ", theTitle))
  message("add covariates")
  for(myCovariate in colnames(theBatchData))
  {
    message(myCovariate)
    if("Sample"!=myCovariate)
    {
      myData <- as.vector(unlist(theBatchData[myCovariate]))
      message(length(myData))
      names(myData) <- as.vector(unlist(theBatchData["Sample"]))
      message("covariate length ", length(myData))
      covar <- chmNewCovariate(myCovariate, myData, chmNewColorMap(unique(sort(myData))))
      barvar <- chmNewCovariateBar(covar, thickness=as.integer(20))
      chm <- chmAddCovariateBar(chm, 'column', barvar)
    }
  }
  message("Set Environment LC_CTYPE=en_US.UTF-8")
  Sys.setenv(LC_CTYPE=theEncoding)
  message("chmExportToHTML")
  htmlNgchm <- paste(theOutputFile, ".html", sep="")
  message(htmlNgchm)
  chmExportToHTML(chm, htmlNgchm, overwrite=TRUE, shaidyMapGen=theShaidyMapGen, shaidyMapGenJava=theShaidyMapGenJava,
                  shaidyMapGenArgs=theShaidyMapGenArgs, ngchmWidgetPath=theNgchmWidgetJs)
  message("chmExportToFile")
  result <- chmExportToFile(chm, theOutputFile, overwrite=TRUE, shaidyMapGen=theShaidyMapGen, shaidyMapGenJava=theShaidyMapGenJava, shaidyMapGenArgs=theShaidyMapGenArgs)
  result
}

# ====================================================================================
# ====================================================================================

#hierClustForNgchmCalc<-function(theMatrixGeneData, theDist="pearson", theClust="ward.D2")
hierClustForNgchmCalc<-function(theMatrixGeneData, theDist="pearson", theClust="ward")
{
  message("hierClustForNgchmCalc")
  collateOrigValue<-Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
  Sys.setlocale("LC_COLLATE","C")
  message("Changing LC_COLLATE to C for duration of run")
  ###Do hierarchical clustering
  # do not need to check for NaN, as NaN is an NA
  if(is.null(theMatrixGeneData))
  {
    message("Unable to calculate distanceMatrix--only one dimension")
    return(NULL)
  }
  else
  {
    subMatrix <- theMatrixGeneData[!is.na(rowSums(theMatrixGeneData)),]
    subMatrix <- subMatrix[!is.infinite(rowSums(subMatrix)),]
    if(nrow(subMatrix)>0)
    {
      message("calculating HC")
      d <- NULL
      tryCatch(
        d<-distanceMatrix(subMatrix, metric=theDist),
        error = function(e)
        {
          message("distanceMatrix failed")
          d <- NULL
        },
        warning = function(e)
        {
          message("distanceMatrix failed")
          d <- NULL
        })

      #message("checking distance matrix results")
      #subMatrix <- d[!is.na(rowSums(d)),]
      #subMatrix <- subMatrix[!is.infinite(rowSums(subMatrix)),]
      uDend <- NULL
      tryCatch(
        uDend<-hclust(d, method=theClust),
        error=function(e)
        {
          message("Unable to calculate hclust--too many NAs, Infinities or NaNs in data")
          uDend <- NULL
        })
      return(uDend)
    }
    else
    {
      message("Unable to calculate distanceMatrix--too many NAs, Infinities or NaNs in data")
      return(NULL)
    }
  }
}

makeGenericColorMap <- function(theMatrix)
{
  # Breakpoint 1 = 95th percentile
  # Breakpoint 2 = 50th percentile (median)
  # Breakpoint 3 = 5th percentile
  breakpoints <- quantile(theMatrix, c(.05, .50, .95), na.rm=TRUE)
  message("breakpoints=", breakpoints)
  if (sum(duplicated(breakpoints))>0)
  {
    myMin <- min(theMatrix, na.rm=TRUE)
    myMedian <- median(theMatrix, na.rm=TRUE)
    myMax <- max(theMatrix, na.rm=TRUE)
    divided <- (myMax+myMin)/5
    message("myMin=",myMin)
    message("myMedian=",myMedian)
    message("myMax=",myMax)
    message("divided=",divided)
    if (0!=divided)
    {
      breakpoints <- c(myMin+divided, myMedian, myMax-divided )
    }
    else
    {
      breakpoints <- c(myMedian-1, myMedian, myMedian+1 )
    }
  }
  breakpoints
}

# ====================================================================================
# ====================================================================================

updateVectorForLinkoutFeatures <- function(theNgchmFeatureMapFile, theVector)
{
  mapping <- readAsGenericDataframe(theNgchmFeatureMapFile)
  currentFeatures <- theVector
  newFeatures <- as.vector(unlist(mapping["linkout"]))
  names(newFeatures) <- as.vector(unlist(mapping["feature"]))
  ##
  indexFromNewFeatures <- match(currentFeatures, names(newFeatures))
  indexFromNewFeatures <- indexFromNewFeatures[!is.na(indexFromNewFeatures)]
  ##
  indexToFeatures <- match(names(newFeatures), currentFeatures)
  indexToFeatures <- indexToFeatures[!is.na(indexToFeatures)]
  ##
  newRowNames <- currentFeatures
  newRowNames[indexToFeatures] <- as.vector(unlist(newFeatures[indexFromNewFeatures]))
  ##
  theVector <- newRowNames
  theVector
}

updateDendrogramForLinkoutFeatures <- function(theNgchmFeatureMapFile, theUDend)
{
  mapping <- readAsGenericDataframe(theNgchmFeatureMapFile)
  currentFeatures <- theUDend$labels
  newFeatures <- as.vector(unlist(mapping["linkout"]))
  names(newFeatures) <- as.vector(unlist(mapping["feature"]))
  ##
  indexFromNewFeatures <- match(currentFeatures, names(newFeatures))
  indexFromNewFeatures <- indexFromNewFeatures[!is.na(indexFromNewFeatures)]
  ##
  indexToFeatures <- match(names(newFeatures), currentFeatures)
  indexToFeatures <- indexToFeatures[!is.na(indexToFeatures)]
  ##
  newRowNames <- theUDend$labels
  newRowNames[indexToFeatures] <- as.vector(unlist(newFeatures[indexFromNewFeatures]))
  ##
  theUDend$labels <- newRowNames
  theUDend
}

updateForLinkoutFeatures <- function(theNgchmFeatureMapFile, theMatrixData)
{
  mapping <- readAsGenericDataframe(theNgchmFeatureMapFile)
  currentFeatures <- rownames(theMatrixData)
  newFeatures <- as.vector(unlist(mapping["linkout"]))
  names(newFeatures) <- as.vector(unlist(mapping["feature"]))
  ##
  indexFromNewFeatures <- match(currentFeatures, names(newFeatures))
  indexFromNewFeatures <- indexFromNewFeatures[!is.na(indexFromNewFeatures)]
  ##
  indexToFeatures <- match(names(newFeatures), currentFeatures)
  indexToFeatures <- indexToFeatures[!is.na(indexToFeatures)]
  ##
  newRowNames <- rownames(theMatrixData)
  newRowNames[indexToFeatures] <- as.vector(unlist(newFeatures[indexFromNewFeatures]))
  ##
  rownames(theMatrixData) <- newRowNames
  theMatrixData
}

# ====================================================================================
# ====================================================================================

buildBatchHeatMapFromHC_Structures <- function(theMatrixData, theBatchData,
                                               theTitle, theOutputFile,
                                               theColDendRDataFile,
                                               theRowDendRDataFile,
                                               theNgchmFeatureMapFile,
                                               theRowType="scholar", theColType="bio.tcga.barcode.sample",
                                               theRowCluster=NULL,
                                               theShaidyMapGen,
                                               theNgchmWidgetJs,
                                               theShaidyMapGenJava="/usr/bin/java",
                                               theShaidyMapGenArgs="-Xmx16G",
                                               theEncoding="en_US.UTF-8")
{
  message(mbatchUtilVersion())
  # need to do this, since NGCHM uses title as filename, so slashed become directories.
  theTitle <- compressIntoFilename(theTitle)
  message("buildBatchHeatMapFromHC_Structures")
  print(theTitle)
  print(theOutputFile)
  dir.create(dirname(theOutputFile), showWarnings=FALSE, recursive=TRUE)
  message("theMatrixData size ", dim(theMatrixData)[1], " ", dim(theMatrixData)[2])
  message("theBatchData size ", dim(theBatchData)[1], " ", dim(theBatchData)[2])
  message("original row names")
  message(length(rownames(theMatrixData)))
  print(rownames(theMatrixData)[1:10])
  message("compute row clusters")
  #print(rev(rownames(theMatrixData)))
  rowClusters <- rownames(theMatrixData)
  if (is.null(theRowDendRDataFile))
  {
    if (!is.null(theRowCluster))
    {
      # c("pearson", "ward.D2")
      message("theRowCluster[1]=", theRowCluster[1])
      message("theRowCluster[2]=", theRowCluster[2])
      rowClusters <- hierClustForNgchmCalc(t(theMatrixData), theRowCluster[1], theRowCluster[2])
      if (is.null(rowClusters))
      {
        message("got NULL, use rownames")
        rowClusters <- rownames(theMatrixData)
      }
      else
      {
        message("cast as.dendrogram")
        if (!is.null(theNgchmFeatureMapFile))
        {
          message("update rowClusters for linkouts")
          rowClusters <- updateDendrogramForLinkoutFeatures(theNgchmFeatureMapFile, rowClusters)
        }
        rowClusters <- as.dendrogram(rowClusters)
      }
    }
  }
  else
  {
    # loaded from file, but needed for check
    uDend <- NULL
    load(theRowDendRDataFile)
    if (!is.null(uDend))
    {
      message("use pre-calc theRowDendRDataFile")
      if (!is.null(theNgchmFeatureMapFile))
      {
        message("update row dendrogram for linkouts")
        uDend <- updateDendrogramForLinkoutFeatures(theNgchmFeatureMapFile, uDend)
      }
      rowClusters <-as.dendrogram(uDend)
    }
  }
  message("compute col clusters")
  colClusters <- colnames(theMatrixData)
  if (!is.null(theColDendRDataFile))
  {
    # loaded from file, but needed for check
    uDend <- NULL
    load(theColDendRDataFile)
    if (!is.null(uDend))
    {
      message("use pre-calc HC clustering")
      colClusters <-as.dendrogram(uDend)
    }
  }
  message("duplicate rownames theMatrixData before linkout")
  message(sum(duplicated(rownames(theMatrixData))))
  message("duplicate colnames theMatrixData before linkout")
  message(sum(duplicated(colnames(theMatrixData))))
  if (!is.null(theNgchmFeatureMapFile))
  {
    message("update features for linkouts")
    message(theNgchmFeatureMapFile)
    theMatrixData <- updateForLinkoutFeatures(theNgchmFeatureMapFile, theMatrixData)
  }
  message("duplicate rownames theMatrixData after linkout")
  message(sum(duplicated(rownames(theMatrixData))))
  message("duplicate colnames theMatrixData after linkout")
  message(sum(duplicated(colnames(theMatrixData))))
  message("make colormapBreaks")
  centeredMatrixData <- bidirectionalCentering(theMatrixData)
  colormapBreaks1 <- makeGenericColorMap(centeredMatrixData)
  colormapBreaks2 <- makeGenericColorMap(theMatrixData)
  print(colormapBreaks1)
  print(colormapBreaks2)
  message("make color map")
  # High = red
  # Medium (usually corresponding to the median) = light gray
  # Low = blue
  # Missing = white
  cmap1 <- chmNewColorMap(colormapBreaks1, colors = c('blue','lightgrey', 'red'), missing.color='white');
  cmap2 <- chmNewColorMap(colormapBreaks2, colors = c('blue','lightgrey', 'red'), missing.color='white');
  message("make data layer")
  message("duplicate rownames")
  message(sum(duplicated(rownames(centeredMatrixData))))
  message("duplicate colnames")
  message(sum(duplicated(colnames(centeredMatrixData))))
  message("centered col names")
  message(length(colnames(centeredMatrixData)))
  print(colnames(centeredMatrixData)[1:10])
  message("centered row names")
  message(length(rownames(centeredMatrixData)))
  print(rownames(centeredMatrixData)[1:10])
  message("original col names")
  print(length(colnames(theMatrixData)))
  print(colnames(theMatrixData)[1:10])
  message("original row names")
  print(length(rownames(theMatrixData)))
  print(rownames(theMatrixData)[1:10])
  layer1 <- chmNewDataLayer('bi-directional median centered', centeredMatrixData, colors=cmap1, summarizationMethod = "average")
  layer2 <- chmNewDataLayer('original', theMatrixData, colors=cmap2, summarizationMethod = "average")
  message("make CHM")
  message("theRowType")
  message(theRowType)
  message("theColType")
  message(theColType)
  # try not using clusters to see if that breaks linkouts
  chm <- chmNew(name=theTitle, layer1, layer2, rowOrder=rowClusters, colOrder=colClusters, rowAxisType=theRowType, colAxisType=theColType)
  message("add caption")
  chm <- chmAddProperty(chm, "chm.info.caption", paste("TCGA heatmap: ", theTitle))
  message("add covariates")
  for(myCovariate in colnames(theBatchData))
  {
    message(myCovariate)
    if("Sample"!=myCovariate)
    {
      myData <- as.vector(unlist(theBatchData[myCovariate]))
      message(length(myData))
      names(myData) <- as.vector(unlist(theBatchData["Sample"]))
      message("covariate length ", length(myData))
      #print(myData)
      covar <- chmNewCovariate(myCovariate, myData, chmNewColorMap(unique(sort(myData))))
      barvar <- chmNewCovariateBar(covar, thickness=as.integer(20))
      chm <- chmAddCovariateBar(chm, 'column', barvar)
    }
  }
  message("Set Environment LC_CTYPE=", theEncoding)
  Sys.setenv(LC_CTYPE=theEncoding)
  message("chmExportToHTML")
  htmlNgchm <- paste(theOutputFile, ".html", sep="")
  message(htmlNgchm)
  chmExportToHTML(chm, htmlNgchm, overwrite=TRUE, shaidyMapGen=theShaidyMapGen, shaidyMapGenJava=theShaidyMapGenJava,
                  shaidyMapGenArgs=theShaidyMapGenArgs, ngchmWidgetPath=theNgchmWidgetJs)
  message("chmExportToFile")
  result <- chmExportToFile(chm, theOutputFile, overwrite=TRUE, shaidyMapGen=theShaidyMapGen, shaidyMapGenJava=theShaidyMapGenJava, shaidyMapGenArgs=theShaidyMapGenArgs)
  result
}

# ====================================================================================
# ====================================================================================

# args <- commandArgs(TRUE)
# print(args)
# THE_MATRIX <- args[1]
# THE_BATYPE <- args[2]
# THE_OUTDIR <- args[3]
# THE_GENJAR <- args[4]
# THE_JAVABIN <- args[5]
#
# processSingleStdFile(THE_MATRIX, THE_BATYPE, THE_OUTDIR, THE_GENJAR, THE_JAVABIN)

# ====================================================================================
# ====================================================================================
