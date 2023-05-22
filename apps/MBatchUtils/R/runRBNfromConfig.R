
# once you are done, I'll show you how to add the R man page for this
runRBNfromConfig <- function(theConfigFile, theOutputDir, theLogDir,
                             theMatrixFile, theMatrixFile2,
                             theDataVersion, theTestVersion)
{
  collateOrigValue<-Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE",collateOrigValue), add=TRUE)
  Sys.setlocale("LC_COLLATE","C")
  message("Changing LC_COLLATE to C for duration of run")
  ####################################################################
  myConfig <- read.csv(theConfigFile, header=FALSE, sep="\t", as.is=TRUE, row.names=1, quote="" )
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
  logFile <- cleanFilePath(theLogDir, "mbatch.log")
  setLogging(new("Logging", theFile=logFile))
  ####################################################################
  matrix1 <- readAsGenericMatrix(theMatrixFile)
  matrix2 <- readAsGenericMatrix(theMatrixFile2)
  RBN_Replicates(matrix1, matrix2,
                 theInvariantGroupId=RBN_InvariantId,
                 theVariantGroupId=RBN_VariantId,
                 theMatchedReplicatesFlag=RBN_Matched,
                 theCombineOnlyFlag=FALSE,
                 thePath=theOutputDir,
                 theDataVersion=theDataVersion,
                 theTestVersion=theTestVersion,
                 theWriteToFile=TRUE)
}
