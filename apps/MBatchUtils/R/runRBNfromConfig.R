
# once you are done, I'll show you how to add the R man page for this
runRBNfromConfig <- function(theConfigFile, theOutputDir)
{
  collateOrigValue<-Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
  Sys.setlocale("LC_COLLATE","C")
  message("Changing LC_COLLATE to C for duration of run")
  ####################################################################
  myConfig <- read.csv(theConfigFile, header=FALSE, sep="\t", as.is=TRUE, row.names=1 )
  ####################################################################
  RBN_UseFirstAsInvariant <- as.logical(convertNulls(convertNA(myConfig["RBN_UseFirstAsInvariantFlag",])))
  RBN_InvariantId <- convertNulls(myConfig["RBN_InvariantId",])
  RBN_VariantId <- convertNulls(myConfig["RBN_VariantId",])
  RBN_Matched <- as.logical(convertNulls(myConfig["RBN_Matched",]))
  RBN_InvariantRepsType <- convertNulls(myConfig["RBN_InvariantRepsType",])
  RBN_VariantRepsType <- convertNulls(myConfig["RBN_VariantRepsType",])
  RBN_InvariantReps <- convertStringArray(convertNulls(myConfig["RBN_InvariantRepsArray",]))
  RBN_VariantReps <- convertStringArray(convertNulls(myConfig["RBN_VariantRepsArray",]))
  message("RBN_UseFirstAsInvariant ", RBN_UseFirstAsInvariant)
  message("RBN_InvariantId ", RBN_InvariantId)
  message("RBN_VariantId ", RBN_VariantId)
  message("RBN_Matched ", RBN_Matched)
  message("RBN_InvariantRepsType ", RBN_InvariantRepsType)
  message("RBN_VariantRepsType ", RBN_VariantRepsType)
  message("RBN_InvariantReps ", RBN_InvariantReps)
  message("RBN_VariantReps ", RBN_VariantReps)
  ####################################################################
  sourceDir <-dirname(theConfigFile)
  logFile <- file.path(sourceDir, "mbatch.log")
  datFile <- file.path(sourceDir, "matrix_data.tsv")
  datFile2 <- file.path(sourceDir, "matrix_data2.tsv")
  #############################################################################
  # check directories
  message("theOutputDir=", theOutputDir)
  if(!dir.exists(theOutputDir))
  {
    message("create ", theOutputDir)
    dir.create(theOutputDir, showWarnings=FALSE, recursive=TRUE)
  }
  #############################################################################
  setLogging(new("Logging", theFile=logFile))
  ####################################################################
  # don't forget to call the RBN function
}
