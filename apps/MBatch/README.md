# MBatch R Package

This is for educational and research purposes only. 

Samples from large research projects are often processed and run in multiple batches at different times. Because the samples are processed in batches rather than all at once, the data can be vulnerable to systematic noise such as batch effects (unwanted variation between batches) and trend effects (unwanted variation over time), which can lead to misleading analysis results.

The MBatch R package is designed to help assess and correct for batch effects. It first allows the user to assess and quantify the presence of any batch effects via algorithms such as Hierarchical Clustering, Principal Component Analysis, and box plots. If significant batch effects are observed in the data, the user then has the option of selecting from a variety of correction algorithms, such as Empirical Bayes (aka Combat), ANOVA and Median Polish.

Additional information can be found at http://bioinformatics.mdanderson.org/main/TCGABatchEffects:Overview

The documentation directort contains several kinds of documentation for MBatch:

 * Files that start MBatch_01 are install documentations for Linux (Debian 9.1), Windows, and OS X.
 * Files that start MBatch_02 are additional details about the test files in the package.
 * Files that start MBatch_03 are detail the file formats used by MBatch and the associated "Standardized Data" files.
 * Files that start MBatch_04 are documentation of assessment algorithms/plots.
 * Files that start MBatch_05 are documentation of correction algorithms.

Downloads and details on Standardized Data are available at http://bioinformatics.mdanderson.org/TCGA/databrowser/

See main README.MD on install instructions.

# Seurat RDS files into Standardized Data Format

Here is example code that reads in an RDS file and converts it to usable Standardized Data format.

Notice in Build the Data Matrix, that you need to know which assay/slot to use. The three standard assay/slots are counts, data, and scale.data. (Scale.data is a technique that is supposed to normalize single cell data for quality – supposed to help with doublets, two or more cells processed as one, and empty results, where no cell was processed and you got random data. Anna can/will explain this technique better.)

Note that the batches dataframe has columns that are not useful as batches. 

```r
# READ THE SEURAT FILE
library(Seurat)
seuObj <- readRDS("~/Downloads/local.rds")

# BUILD THE DATA MATRIX
assayObj <- seuObj@assays[[assayToUse]]
# select slot with matrix info (usually: counts, data, scale.data)
assayToUse <- "RNA"
assayData <- slot(assayObj,"data")
assayMatrix <- as.matrix(assayData)

# BUILD THE BATCHES DATAFRAME
metaData <- seuObj@meta.data
# convert matrix to dataframe
metaDF <- as.data.frame(metaData)
# add label for first column
metaDF$Sample <- rownames(metaDF)
# move Row column to the front
metaDF <- metaDF[, c(ncol(metaDF), 1:(ncol(metaDF)-1))]
# remove rownames
rownames(metaDF) <- NULL

# SAVE IN STANDARDIZED DATA FORMAT
library(MBatch)
writeAsGenericMatrix("/scratch/tdcasasent/DEV/seurat_matrix.tsv", assayMatrix)
writeAsGenericDataframe("/scratch/tdcasasent/DEV/seurat_batches.tsv", metaDF)
```

See also the functions: convertSeuratObjtoStdDataObj, convertSeuratRDStoStdDataFiles, updateOrBuildSeuratObjectFromStdDataFiles, and updateOrBuildSeuratObjectFromStdDataObjs.

# Accessing Standardized Data from R

Here is example code for using the MBatch API for finding and downloading data.

```r
# URL of MBatch Standardized Data API endpoint
myUrl <- "https://bioinformatics.mdanderson.org/MQA"
# build base Python object for performing query
pyObj <- getDapiQuery(myUrl)
# update the query to look for TCGA LUSC, uncorrected data, of the "STAR - Counts" data type.
pyObj$selected_projects <- append(pyObj$selected_projects, "TCGA-LUSC")
pyObj$selected_jobtype <- append(pyObj$selected_jobtype, "Original")
pyObj$selected_data <- append(pyObj$selected_data, "STAR - Counts")
# update the query results to contain just the requested datasets
updateDapiQuery(pyObj)
# should be 3
length(pyObj$available_datasets)
# list is 1 indexed for Python object in R, also remove object from list
queryEntry <- pyObj$available_datasets[[1]]
retList <- getDownloadableForDataset(pyObj, queryEntry)
# should be 1
versionVector <- retList[[1]]
length(versionVector)
# should be 12
ngchmVector <- retList[[2]]
length(ngchmVector)
# download matrix files
# id and version are hard coded for testing, but available from the results
datasetId <- "8ff96a845261f8455b5b2698b5db776f~2022-12-12~2022_12_28_1300"
versionStr <- "DATA_2022-12-12"
downloadFile <- file.path(theOutputDir, "from_r_matrix_original.tsv")
downloadDataMatrix(pyObj, downloadFile, datasetId, versionStr, TRUE)
downloadFile <- file.path(theOutputDir, "from_r_matrix_pipeline.tsv")
downloadDataMatrix(pyObj, downloadFile, datasetId, versionStr, FALSE)
downloadFile <- file.path(theOutputDir, "from_r_batches_original.tsv")
downloadDataBatches(pyObj, downloadFile, datasetId, versionStr, TRUE)
downloadFile <- file.path(theOutputDir, "from_r_batches_pipeline.tsv")
downloadDataBatches(pyObj, downloadFile, datasetId, versionStr, FALSE)
# download batch files
# id is hard coded for testing, but available from the results
datasetId <- "8ff96a845261f8455b5b2698b5db776f~2022-12-12~2022_12_28_1300"
downloadFile <- file.path(theOutputDir, "batch_id_ngchm.ngchm")
zipFilePath <- "/analysis/NGCHM/DATA_2022-12-12/TEST_2022_12_28_1300/All_ngchm.ngchm.html"
downloadNgchmNgchm(pyObj, downloadFile, datasetId, zipFilePath)
downloadFile <- file.path(theOutputDir, "batch_id_ngchm.html")
downloadNgchmHtml(pyObj, downloadFile, datasetId, zipFilePath)
```



