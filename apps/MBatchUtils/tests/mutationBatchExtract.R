#mutationBatchExtract


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
  baseTestDir=getTestInputDir()
  baseOutputDir=getTestOutputDir()
  dccMinMaf=file.path(baseTestDir, "DCCMin_MAF")
  gdcMinMaf=file.path(baseTestDir, "GDCMin_MAF")
  dccMinMut=file.path(baseOutputDir, "DCCMin_MUT")
  gdcMinMut=file.path(baseOutputDir, "GDCMin_MUT")
  dccMinMutCompare=file.path(baseTestDir, "DCCMin_MUT")
  gdcMinMutCompare=file.path(baseTestDir, "GDCMin_MUT")
  ########################################################
  ########################################################
  unlink(dccMinMut, recursive=TRUE)
  dir.create(dccMinMut, showWarnings=FALSE, recursive=TRUE)
  unlink(gdcMinMut, recursive=TRUE)
  dir.create(gdcMinMut, showWarnings=FALSE, recursive=TRUE)
  ########################################################
  ########################################################
  mutationBatchExtract(theMafDir=dccMinMaf, theGDCflag=FALSE, theTypeCountDir=dccMinMut)
  mutationBatchExtract(theMafDir=gdcMinMaf, theGDCflag=TRUE, theTypeCountDir=gdcMinMut)
  ########################################################
  ########################################################
  correctedMatrix <- readAsGenericMatrix(file.path(dccMinMut, "acc", "acc.illuminaga.hgsc_bcm_edu.Level_2_HG19_Total.tsv"))
  compareMatrix <- readAsGenericMatrix(file.path(dccMinMutCompare, "acc", "acc.illuminaga.hgsc_bcm_edu.Level_2_HG19_Total.tsv"))
  print("compared1")
  compared1 <- compareTwoMatrices(correctedMatrix, compareMatrix)
  print(compared1)
  correctedMatrix <- readAsGenericMatrix(file.path(gdcMinMut, "TCGA-ACC", "TCGA-ACC.MuSEVariantAggregationandMasking_HG38_Total.tsv"))
  compareMatrix <- readAsGenericMatrix(file.path(gdcMinMutCompare, "TCGA-ACC", "TCGA-ACC.MuSEVariantAggregationandMasking_HG38_Total.tsv"))
  print("compared2")
  print(file.path(gdcMinMut, "TCGA-ACC", "TCGA-ACC.MuSEVariantAggregationandMasking_HG38_Total.tsv"))
  print(file.path(gdcMinMutCompare, "TCGA-ACC", "TCGA-ACC.MuSEVariantAggregationandMasking_HG38_Total.tsv"))
  compared2 <- compareTwoMatrices(correctedMatrix, compareMatrix)
  print(compared2)
  compared1&&compared2
} else {
  message("No test data. Skip test.")
  TRUE
}
