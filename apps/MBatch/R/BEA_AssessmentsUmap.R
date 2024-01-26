# MBatch Copyright (c) 2011-2024 University of Texas MD Anderson Cancer Center
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
# MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>

embed_img <- function(X, Y, k = 15, ...)
{
  args <- list(...)
  args$coords <- Y
  args$x <- X

  do.call(vizier::embed_plot, args)
}

writeErrorsForAllBatchTypesUMAP<-function(theMessage, theDataframeBatchData, theUmapOutputDir,
                                          theTitle, theDataVersion, theTestVersion)
{
  unlink(theUmapOutputDir, recursive = TRUE, force = TRUE)
  dir.create(theUmapOutputDir, showWarnings=FALSE, recursive=TRUE)
  for(batchTypeIndex in c(2:length(theDataframeBatchData)))
  {
    ### compile data and information for display
    batchTypeName <- names(theDataframeBatchData)[batchTypeIndex]
    logInfo("createBatchEffectsOutput_umap batchTypeName=", batchTypeName)
    umapBatchDir <- checkCreateDir(theUmapOutputDir, batchTypeName)
    umapBatchDir <- addVersionsIfNeeded(umapBatchDir, theDataVersion, theTestVersion)
    # do errors for each batch type
    dir.create(umapBatchDir, showWarnings=FALSE, recursive=TRUE)
    mesg=paste("No UMAP -", theMessage, " - ", theTitle, sep="")
    openAndWriteIssuesUmap(umapBatchDir, mesg)
  }
}

createBatchEffectsOutput_umap<-function(theMatrix, theDataframeBatchData, theTitle,
                                        theDataVersion, theTestVersion,
                                        theUmapOutputDir,
                                        theDoDSCFlag, theDSCPermutations,
                                        theDSCThreads, theDoDscPermsFileFlag, theSeed)
{
  logDebug("createBatchEffectsOutput_umap")
  collateOrigValue<-Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
  Sys.setlocale("LC_COLLATE","C")
  logDebug("Changing LC_COLLATE to C for duration of run")
  checkPackageSettings()
  stopifnotWithLogging("The number of columns of gene data must equal the number of rows of batch data.", ncol(theMatrix)==nrow(theDataframeBatchData))
  rdataFiles <- NULL
  tryCatch({
    pca_data <- doPrcompCall(theMatrix)
    if (is.null(pca_data))
    {
      stop("Unable to calculate--too many NAs, Infinities or NaNs in data")
    }
    else
    {
      neighbors <- 15
      if (neighbors > (dim(pca_data$x)[2]))
      {
        neighbors <- dim(pca_data$x)[2] -1
        if (neighbors < 1)
        {
          neighbors <- 1
        }
      }
      dscOutputDir <- NULL
      if (TRUE==theDoDSCFlag)
      {
        # make directory for later use
        dscOutputDir <- addVersionsIfNeeded(cleanFilePath(theUmapOutputDir, "DSC"), theDataVersion, theTestVersion)
        checkDirForCreation(dscOutputDir)
        file.create(file.path(dscOutputDir, "DSCOverview.tsv"))
        logDebug("make future dir", dscOutputDir)
        # write complist file
        # NOTE - MUST BE AT TOP LEVEL FOR DSC OVERVIEW FUNCTION TO WORK
        logDebug("write ALL__CompListDSC.RData", theUmapOutputDir)
        saveCompListDscData(cleanFilePath(theUmapOutputDir, "ALL__CompListDSC.RData"), c("PC1", "PC2", "PC3", "PC4"))
      }
      logInfo("createBatchEffectsOutput_umap neighbors=", neighbors)
      umap_data <- uwot::umap(pca_data$x, n_neighbors = neighbors)
      # do not recast to dataframe, it changes the column names in a way that breaks UTF-8 characters
      #myDF <- data.frame(theDataframeBatchData, stringsAsFactors=FALSE)
      samplesIds <- as.vector(unlist(theDataframeBatchData[1]))
      for(batchTypeIndex in c(2:length(theDataframeBatchData)))
      {
        ### compile data and information for display
        batchTypeName <- names(theDataframeBatchData)[batchTypeIndex]
        logInfo("createBatchEffectsOutput_umap batchTypeName=", batchTypeName)
        umapBatchDir <- checkCreateDir(theUmapOutputDir, batchTypeName)
        umapBatchDir <- addVersionsIfNeeded(umapBatchDir, theDataVersion, theTestVersion)
        checkDirForCreation(umapBatchDir)
        logInfo("createBatchEffectsOutput_umap umapBatchDir=", umapBatchDir)
        batchIdsForSamples <- theDataframeBatchData[batchTypeName]
        # needs to be factors for embed_img/embed_plot to work
        batchIdsForSamples <- as.vector(unlist(batchIdsForSamples))
        uniqueBatchIds <- sort(unique(sort(batchIdsForSamples)))
        batchIdsForSamples <- factor(batchIdsForSamples, levels=uniqueBatchIds)
        # ############################
        dscAllResults <- NULL
        if(TRUE==theDoDSCFlag)
        {
          logDebug("call pvalueDSC")
          pca_scores <- SamplePCA(t(umap_data))
          dscAllResults <- pvalueDSC(pca_scores, batchIdsForSamples, theDSCPermutations, 0, 0, theDSCThreads, theSeed)
          logDebug("openAndWriteDscAllFile 2")
          openAndWriteDscAllFile(pca_scores, dscAllResults, umapBatchDir, "ANY", "NA", "NA")
          if (TRUE==theDoDscPermsFileFlag)
          {
            openAndWriteDscPermsFile(dscAllResults, umapBatchDir, "ANY", "NA", "NA")
          }
          logDebug("after openAndWriteDscAllFile 2")
        }
        # ############################
        logDebug("UMAP diagrams")
        diagramFilename = createDirPlusFilename(umapBatchDir, "UMAP_Diagram.png")
        mytitle <- paste(theTitle, "/", "UMAP", batchTypeName, sep=" ")
        writeTitleFile(mytitle, diagramFilename)
        legendFilename = createDirPlusFilename(umapBatchDir, "UMAP_Legend.png")
        color <- beaRainbow(length(uniqueBatchIds), v=0.7)
        title <- paste(theTitle, "UMAP (PCA)", batchTypeName, sep=" ")
        umapDiagram(diagramFilename, batchIdsForSamples, umap_data, color, title, levels(batchIdsForSamples))
        umapLegend(legendFilename, title, color, uniqueBatchIds, as.vector(unlist(batchIdsForSamples)))
        meta <- data.frame(key=c("neighbors"), value=c(neighbors))
        rdataFiles <- c(rdataFiles, writeUmapDataTSVs(samplesIds, as.vector(unlist(batchIdsForSamples)),
                                                      umap_data, umapBatchDir, "UMAP_Data-", meta))
        for(i in seq(from=1, to=2, by=2))
        {
          componentA <- colnames(pca_data$x)[[i]]
          componentB <- colnames(pca_data$x)[[i+1]]
          ###
          dscTxtFile <- makePcaFileName_TXT(umapBatchDir, "ANY", componentA, componentB, "DSC")
          errorFile <- file.path(umapBatchDir, "error.log")
          ###
          results <- NULL
          if (TRUE==theDoDSCFlag)
          {
            logDebug("doInternalPca - pvalueDSC 2")
            ### TODO: make calls to doInternalPca return the list of DSC objects created here
            pca_scores <- SamplePCA(t(umap_data))
            # use i and i+1, since component variables are strings here
            compA <- i
            compB <- i+1
            results <- pvalueDSC(pca_scores, batchIdsForSamples, theDSCPermutations, compA, compB, theDSCThreads, theSeed)
            if(TRUE==theDoDscPermsFileFlag)
            {
              openAndWriteDscPermsFile(results, dirname(dscTxtFile), "ANY", compA, compB)
            }
            diagramLabelDsc <- dscAllResults@mDSC
            diagramLabelDscXY <- results@mDSC
            DSCLabels <- c(
              "First PCA Component:",
              "First Component FVE (%):",
              "Second PCA Component:",
              "Second Component FVE (%):",
              "",
              paste("Disp. Sep. Crit. (DSC) (", compA, ",", compB, "):", sep=""),
              paste("Disp. within groups (Dw) (", compA, ",", compB, "):", sep=""),
              paste("Disp. between groups (Db) (", compA, ",", compB, "):", sep=""),
              paste("DSC pvalue(", compA, ",", compB, "):", sep=""),
              "",
              "Disp. Sep. Crit. (DSC):",
              "Disp. within groups (Dw):",
              "Disp. between groups (Db):",
              "DSC pvalue:")
            DSCvalues <- c(
              componentA,
              calculateFVEForSingleComponent(pca_scores, compA)*100,
              componentB,
              calculateFVEForSingleComponent(pca_scores, compB)*100,
              "",
              results@mDSC,
              results@mDW,
              results@mDB,
              results@mPvalue,
              "",
              dscAllResults@mDSC,
              dscAllResults@mDW,
              dscAllResults@mDB,
              dscAllResults@mPvalue)
            ### combine DSC and DSC pvalue data
            extraLegend <- combineStringPairs(DSCLabels, DSCvalues, " ", FALSE, FALSE, FALSE)
            combineForFile <- combineStringPairs(DSCLabels, DSCvalues, "\t", FALSE, FALSE, FALSE)
            ### save value objects for later reuse in making DSC summary / passing in PCA-PVALUE-DSC / file is *__CompDSC.RData
            saveCompDscData(paste(dscTxtFile, "__CompDSC.RData", sep=""), results, compA, compB)
          }
          if(TRUE==theDoDSCFlag)
          {
            openAndWriteDSCFile(combineForFile, dscTxtFile)
          }
        }
      }
    }
  }, error = function(err)
  {
    logError("2 createBatchEffectsOutput_umap caught an error", err)
    rdataFiles <- NULL
    writeErrorsForAllBatchTypesUMAP(err, theDataframeBatchData, theUmapOutputDir,
                                    theTitle, theDataVersion, theTestVersion)
  })
  rdataFiles
}

umapLegend<-function(theFilename, theTitle, theColors, theUniqueBatches, theAllBatches)
{
  logDebug("umapLegend", theFilename)
  version <- paste("MBatch", packageDescription("MBatch")["Version"], sep=" ")
  stopifnotWithLogging("The number of sorted colors and number of sorted batch ids must match.",
                       length(theUniqueBatches)==length(theColors))
  sortedBatchIdsWithCount <- as.vector(unlist(sapply(theUniqueBatches, function (batchId)
  {
    return(paste(batchId," (",sum(theAllBatches==batchId), ")", sep=""))
  }, USE.NAMES=FALSE)))
  ##leg end("center", legend=sortedBatchIdsWithCount, fill=sortedColors, title=batchTypeName)
  mbatchStandardLegend(theTitle, version, sortedBatchIdsWithCount, theColors, NULL, theFilename)
  return(theFilename)
}

umapDiagram <- function(theFilename, theBatches, theUmap, theColors, theTitle, theLimits)
{
  logDebug("umapDiagram - start")
  CairoPNG(filename=theFilename, width = 1000, height = 1000, pointsize=24, bg = "white")
  on.exit(dev.off(), add = TRUE)
  embed_img(theBatches, theUmap, pc_axes=TRUE, equal_axes=TRUE, alpha_scale=0.5,
            title=theTitle, cex=1, colors=theColors, limits=theLimits)
}

doPrcompCall<-function(theMatrixGeneData, theRank=50)
{
  logDebug("doPrcompCall - start")
  checkIfTestError()
  pca <- NULL
  if (sum(!is.na(rowSums(theMatrixGeneData)))<2)
  {
    logWarn("All rows contain NA values")
  }
  else
  {
    ### calculate prcomp
    errorsContinued <- FALSE
    logOutput("keep only NA")
    myData <- theMatrixGeneData[!is.na(rowSums(theMatrixGeneData)),]
    logOutput("remove zero variance")
    myData <- myData[apply(myData, 1, var)!=0,]
    pca <- try( prcomp(t(myData), scale=TRUE, center=TRUE, rank=theRank) )
    if(is(pca, 'try-error'))
    {
      logDebug("prcomp first call threw a problem -- dropping outliers and retrying")
      myNoMissingMatrixGeneData<-theMatrixGeneData[!is.na(rowSums(theMatrixGeneData)),]
      myGeneVar<-apply(myNoMissingMatrixGeneData, 1, var)
      sortedGeneVar<-sort(myGeneVar, decreasing=TRUE)
      myOutlierNumber<-1
      logDebug("myOutlierNumber =", myOutlierNumber)
      while((is(pca,'try-error'))&&(myOutlierNumber<length(sortedGeneVar)))
      {
        logDebug("prcomp subsequent call generated a failed (not necessarily a problem yet) ", myOutlierNumber)
        subMatrix <- myNoMissingMatrixGeneData[myGeneVar<sortedGeneVar[myOutlierNumber],]
        pca<-try(prcomp(t(subMatrix), scale=TRUE, center=TRUE, rank=theRank))
        myOutlierNumber<-myOutlierNumber+1
      }
      if (myOutlierNumber>=length(sortedGeneVar))
      {
        logWarn("prcomp continued generating failures (data set cannot be used) ", myOutlierNumber)
        errorsContinued <- TRUE
      }
    }
  }
  if(is.null(pca))
  {
    logDebug("prcomp - pca is null")
  }
  return(pca)
}

umapPCA_calc<-function(theMatrix)
{
  logDebug("umapPCA_calc")
  collateOrigValue<-Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
  Sys.setlocale("LC_COLLATE","C")
  logDebug("Changing LC_COLLATE to C for duration of run")
  checkPackageSettings()
  ###Do hierarchical clustering
  # do not need to check for NaN, as NaN is an NA
  if (is.null(theMatrix))
  {
    logWarn("Unable to calculate prcomp--only one dimension")
    return(NULL)
  }
  else
  {
    subMatrix <- theMatrix[!is.na(rowSums(theMatrix)),]
    subMatrix <- subMatrix[!is.infinite(rowSums(subMatrix)),]
    if (nrow(subMatrix)>0)
    {
      logDebug("calculating distance")
      d<-distanceMatrix(subMatrix, metric="pearson")
      return(d)
    }
    else
    {
      logWarn("Unable to calculate distanceMatrix--too many NAs, Infinities or NaNs in data")
      return(NULL)
    }
  }
}

openAndWriteIssuesUmap<-function(theOutputDir, theMessage)
{
  myFile <- file(cleanFilePath(theOutputDir, "error.log"), "w+")
  on.exit(close(myFile))
  cat(theMessage, file=myFile, append=TRUE)
}

####################################################################
####################################################################

writeUmapDataTSVs<-function(theSamplesIds, theBatchIdsForSamples, theUmapData,
                            theOutputDir, theOutFilePrefix, theMeta)
{
  #rdataFile <- cleanFilePath(theOutputDir, paste(theOutFilePrefix, "all.RData", sep=""))
  # build dataframe from theUmapData with theSamplesIds in front
  writeUmapData <- cbind(Samples=theSamplesIds, theUmapData)
  colnames(writeUmapData) <- c("Samples", "V1", "V2")
  umapFile <- cleanFilePath(theOutputDir, paste(theOutFilePrefix, "umap.tsv", sep=""))
  write.table(writeUmapData, file=umapFile, row.names = FALSE, sep = "\t", quote = FALSE)
  #sampFile <- cleanFilePath(theOutputDir, paste(theOutFilePrefix, "samp.tsv", sep=""))
  #write.table(theSamplesIds, file=sampFile, row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
  writeBatchData <- cbind(Samples=theSamplesIds, Batches=theBatchIdsForSamples)
  colnames(writeBatchData) <- c("Samples", "Batches")
  batcFile <- cleanFilePath(theOutputDir, paste(theOutFilePrefix, "batc.tsv", sep=""))
  write.table(writeBatchData, file=batcFile, row.names = FALSE, sep = "\t", quote = FALSE)
  metaFile <- cleanFilePath(theOutputDir, paste(theOutFilePrefix, "meta.tsv", sep=""))
  write.table(theMeta, file=metaFile, row.names = FALSE, sep = "\t", quote = FALSE)
  #write.table(data, file = cleanFilePath(theHierClustOutputDir,theOutputHCOrderFileName), append = FALSE, quote = FALSE, sep = "\t", row.names=FALSE)
  # write udend RData file
  #save(theSamplesIds, theBatchIdsForSamples, theUmapData, file=rdataFile)
  umapFile
}

####################################################################
####################################################################
