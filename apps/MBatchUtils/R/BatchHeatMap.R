#MBatchUtils Copyright ? 2018 University of Texas MD Anderson Cancer Center
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
# test_ngchm_be.R

library(methods)
library(NGCHM)
library(MBatch)
library(ClassDiscovery)

buildBatchHeatMapForSTDdata <- function(theFile, theBatchType, theOutputDir,
                                 theShaidyMapGen,
                                 theShaidyMapGenJava="/usr/bin/java")
{
  message(theFile)
  dir.create(file.path(theOutputDir,
                       basename(dirname(dirname(dirname(dirname(dirname(theFile)))))),
                       basename(dirname(dirname(dirname(dirname(theFile))))),
                       basename(dirname(dirname(dirname(theFile)))),
                       basename(dirname(dirname(theFile))),
                       basename(dirname(theFile))),
             recursive=TRUE, showWarnings=TRUE)
  batchTypes <- tail(names(readAsDataFrame(file.path(dirname(theFile), "batches.tsv"))), -1)
  message("dvlpStandardizedDataToProcess batchType")
  message(theBatchType)
  myTitle <- paste( theBatchType,
                    basename(dirname(dirname(dirname(dirname(dirname(theFile)))))),
                    basename(dirname(dirname(dirname(dirname(theFile))))),
                    basename(dirname(dirname(dirname(theFile)))),
                    basename(dirname(dirname(theFile))),
                    basename(dirname(theFile)))
  buildBatchHeatMap_Files(theMatrixFile=theFile,
                          theBatchFile=file.path(dirname(theFile), "batches.tsv"),
                          theTitle=myTitle,
                          theOutputFile=file.path(theOutputDir,
                                                  basename(dirname(dirname(dirname(dirname(dirname(theFile)))))),
                                                  basename(dirname(dirname(dirname(dirname(theFile))))),
                                                  basename(dirname(dirname(dirname(theFile)))),
                                                  basename(dirname(dirname(theFile))),
                                                  basename(dirname(theFile)),
                                                  paste(theBatchType, "_matrix.ngchm", sep="")),
                          theSortByType=theBatchType,
                          theShaidyMapGen=theShaidyMapGen,
                          theShaidyMapGenJava=theShaidyMapGenJava,
                          theShaidyMapGenArgs="-Xmx16G")
}

buildBatchHeatMap_Files <- function(theMatrixFile, theBatchFile, theTitle, theOutputFile, theSortByType,
                                    theRowType="labels", theColType="bio.tcga.barcode.sample",
                                    theRowCluster=NULL, theColCluster=NULL,
                                    theShaidyMapGen,
                                    theShaidyMapGenJava="/usr/bin/java",
                                    theShaidyMapGenArgs="-Xmx16G")
{
  message("buildBatchHeatMap_Files")
  message(theMatrixFile)
  matData <- readAsGenericMatrix(theMatrixFile)
  message(theBatchFile)
  dfData <- readAsDataFrame(theBatchFile)
  buildBatchHeatMap_Structures(matData, dfData, theTitle, theOutputFile, theSortByType,
                               theRowType, theColType,
                               theRowCluster, theColCluster,
                               theShaidyMapGen,
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
                                         theRowType="labels", theColType="bio.tcga.barcode.sample",
                                         theRowCluster=NULL, theColCluster=NULL,
                                         theShaidyMapGen,
                                         theShaidyMapGenJava="/usr/bin/java",
                                         theShaidyMapGenArgs="-Xmx16G")
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
  message("make quartiles")
  quartiles <- makeGenericColorMap(theMatrixData)
  print(quartiles)
  message("make color map")
  cmap1 <- chmNewColorMap(quartiles, colors = c('black', 'blue','yellow', 'orange'), missing.color='red');
  message("make data layer")
  layer1 <- chmNewDataLayer('matrix data', theMatrixData, colors=cmap1, summarizationMethod = "average")
  message("make CHM")
  chm <- chmNew(name=theTitle, layer1, rowOrder=rowClusters, colOrder=colClusters, rowAxisType=theRowType, colAxisType=theColType)
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
      covar <- chmNewCovariate(myCovariate, myData)
      chm <- chmAddCovariateBar(chm, 'column', covar)
    }
  }
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
  quantiles <- quantile(theMatrix, c(.25, .50,  .75, .90), na.rm=TRUE)
  if (sum(duplicated(quantiles))>0)
  {
    myMin <- min(theMatrix)
    myMax <- max(theMatrix)
    divided <- (myMax+myMin)/5
    if (0!=divided)
    {
      quantiles <- c(myMin+divided, myMin+divided+divided, myMin+divided+divided+divided, myMin+divided+divided+divided+divided )
    }
    else
    {
      middle <- (myMin+myMax)/2
      quantiles <- c(myMin-1, middle+(myMin/2), middle+(myMax/2), myMax+1 )
    }
  }
  quantiles
}

# ====================================================================================
# ====================================================================================

buildBatchHeatMapFromHC_Structures <- function(theMatrixData, theBatchData,
                                               theTitle, theOutputFile,
                                               theUDendRDataFile,
                                         theRowType="labels", theColType="bio.tcga.barcode.sample",
                                         theRowCluster=NULL,
                                         theShaidyMapGen,
                                         theShaidyMapGenJava="/usr/bin/java",
                                         theShaidyMapGenArgs="-Xmx16G")
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
  message("compute row clusters")
  rowClusters <- rownames(theMatrixData)
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
      rowClusters <- as.dendrogram(rowClusters)
    }
  }
  message("compute col clusters")
  colClusters <- colnames(theMatrixData)
  if (!is.null(theUDendRDataFile))
  {
    # loaded from file, but needed for check
    uDend <- NULL
    load(theUDendRDataFile)
    if (!is.null(uDend))
    {
      message("use HC clustering")
      colClusters <-as.dendrogram(uDend)
    }
  }
  message("make quartiles")
  quartiles <- makeGenericColorMap(theMatrixData)
  print(quartiles)
  message("make color map")
  cmap1 <- chmNewColorMap(quartiles, colors = c('black', 'blue','yellow', 'orange'), missing.color='red');
  message("make data layer")
  layer1 <- chmNewDataLayer('matrix data', theMatrixData, colors=cmap1, summarizationMethod = "average")
  message("make CHM")
  chm <- chmNew(name=theTitle, layer1, rowOrder=rowClusters, colOrder=colClusters, rowAxisType=theRowType, colAxisType=theColType)
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
      covar <- chmNewCovariate(myCovariate, myData)
      chm <- chmAddCovariateBar(chm, 'column', covar)
    }
  }
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
