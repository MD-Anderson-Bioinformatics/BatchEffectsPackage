#mbatchRunFromConfig


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
  # writes to input directory, so copy files to output
  outDir=file.path(getTestOutputDir(), "configout", "2018-07-11-1200")
  print(outDir)
  unlink(outDir, recursive=TRUE)
  dir.create(outDir, recursive=TRUE, showWarnings=FALSE)
  print(file.exists(outDir))
  print(file.path(getTestInputDir(), "config", "MBatchConfig.tsv"))
  print(file.path(outDir, "MBatchConfig.tsv"))
  file.copy(file.path(getTestInputDir(), "config", "MBatchConfig.tsv"), file.path(outDir, "MBatchConfig.tsv"))
  print(file.exists(file.path(outDir, "MBatchConfig.tsv")))
  file.copy(file.path(getTestInputDir(), "config", "matrix_data.tsv"), file.path(outDir, "matrix_data.tsv"))
  file.copy(file.path(getTestInputDir(), "config", "batches.tsv"), file.path(outDir, "batches.tsv"))
  configFile=file.path(outDir, "MBatchConfig.tsv")
  ########################################################
  ########################################################
  mbatchRunFromConfig(theConfigFile=configFile, theOutputDir=outDir, theNaStrings="NA")
  file.exists(file.path(outDir, "MBATCH_SUCCESS.txt"))
} else {
  message("No test data. Skip test.")
  TRUE
}
