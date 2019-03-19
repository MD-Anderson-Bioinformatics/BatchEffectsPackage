#buildBatchHeatMap_Files

library(MBatchUtils)

if (!is.null(getTestOutputDir()))
{
  #warnLevel<-getOption("warn")
  #on.exit(options(warn=warnLevel))
  # warnings are errors
  #options(warn=3)
  # if there is a warning, show the calls leading up to it
  #options(showWarnCalls=TRUE)
  # if there is an error, show the calls leading up to it
  #options(showErrorCalls=TRUE)
  ########################################################
  ########################################################
  jarDir=file.path(getTestInputDir(), "exe")
  baseTestDir=file.path(dirname(getTestInputDir()), "MATRIX_DATA")
  baseOutputDir=getTestOutputDir()
  javaExe=getJava()
  ####
  jarFile=file.path(jarDir, "ShaidyMapGen.jar")
  sourceMatrix=file.path(baseTestDir, "brca_agi4502_matrix_data.tsv")
  sourceBatches=file.path(baseTestDir, "brca_agi4502_batches.tsv")
  dir.create(file.path(baseOutputDir, "NGCHM"), showWarnings=FALSE)
  outputFile=file.path(baseOutputDir, "NGCHM", "buildBatchHeatMap_Files.ngchm")
  ########################################################
  ########################################################
  buildBatchHeatMap_Files(theMatrixFile=sourceMatrix,
                          theBatchFile=sourceBatches,
                          theTitle="Test for BRCA AGI4502 Data",
                          theOutputFile=outputFile,
                          theSortByType="BatchId",
                          theRowType="labels", theColType="bio.tcga.barcode.sample",
                          theRowCluster=NULL, theColCluster=NULL,
                          theShaidyMapGen=jarFile,
                          theShaidyMapGenJava=javaExe,
                          theShaidyMapGenArgs="-Xmx16G")
  TRUE
}
