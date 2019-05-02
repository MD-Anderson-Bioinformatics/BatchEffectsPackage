#mutationBatchAssess


library(MBatchUtils)
library(tools)


if (!is.null(getTestOutputDir()))
{
  warnLevel<-getOption("warn")
  on.exit(options(warn=warnLevel))
  # warnings are errors
  options(warn=3)
  # if there is a warning, show the calls leading up to it
  options(showWarnCalls=TRUE)
  # if there is an error, show the calls leading up to it
  options(showErrorCalls=TRUE)
  ########################################################
  ########################################################
  baseTestDir=getTestInputDir()
  baseOutputDir=getTestOutputDir()
  dccMinMut=file.path(baseTestDir, "DCCMin_MUT")
  gdcMinMut=file.path(baseTestDir, "GDCMin_MUT")
  dccMinOut=file.path(baseOutputDir, "DCCMin_OUT")
  gdcMinOut=file.path(baseOutputDir, "GDCMin_OUT")

  dccCompare=file.path(baseTestDir, "DCCMin_OUT")
  gdcCompare=file.path(baseTestDir, "GDCMin_OUT")
  ########################################################
  ########################################################
  unlink(dccMinOut, recursive=TRUE)
  dir.create(dccMinOut, showWarnings=FALSE, recursive=TRUE)
  unlink(gdcMinOut, recursive=TRUE)
  dir.create(gdcMinOut, showWarnings=FALSE, recursive=TRUE)
  ########################################################
  ########################################################
  mutationBatchAssess(theTypeCountDir=dccMinMut, theOutputDir=dccMinOut, theJavaArgs=c("-Xms16000m", "-Djava.awt.headless=true"), theThreads=5,
                      thePvalueCutoff=.00001, theZScoreCutoff=1.96, thePCAflag=TRUE,
                      theBatchTypes=c("BatchId"),
                      theMutationTypes=c("FrameShiftIns", "Total"))
  mutationBatchAssess(theTypeCountDir=gdcMinMut, theOutputDir=gdcMinOut, theJavaArgs=c("-Xms16000m", "-Djava.awt.headless=true"), theThreads=5,
                      thePvalueCutoff=.00001, theZScoreCutoff=1.96, thePCAflag=TRUE,
                      theBatchTypes=c("PlateId"),
                      theMutationTypes=c("FrameShiftIns", "Total"))
  ########################################################
  ########################################################
  htmlMutationBatchEffects(dccMinOut)
  htmlMutationBatchEffects(gdcMinOut)
  ########################################################
  ########################################################
  outIndex <- file.path(dccMinOut, "BatchId_FrameShiftIns", "index.html")
  cmpIndex <- file.path(gdcCompare, "BatchId_FrameShiftIns", "index.html")
  outCallRef <- file.path(dccMinOut, "BatchId_FrameShiftIns", "callReference.tsv")
  cmpCallRef <- file.path(gdcCompare, "BatchId_FrameShiftIns", "callReference.tsv")
  (md5sum(outCallRef) == md5sum(cmpCallRef))&&(md5sum(outIndex) == md5sum(cmpIndex))
} else {
  message("No test data. Skip test.")
  TRUE
}
