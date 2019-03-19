library(MBatchUtils)

if (!is.null(getTestOutputDir()))
{
  ########################################################
  ########################################################
  # writes to input directory, so copy files to output
  sourceDir=file.path(getTestInputDir(), "configout")
  outDir=file.path(getTestOutputDir(), "archiveout", "2018-07-11-1200")
  unlink(outDir)
  dir.create(outDir, recursive=TRUE, showWarnings=FALSE)
  #file.copy(file.path(getTestInputDir(), "config", "MBatchConfig.tsv"), file.path(outDir, "MBatchConfig.tsv"))
  #buildSingleArchive
  buildSingleArchive(sourceDir, outDir, dirname(outDir),
                     theExternalIndexPath=file.path(outDir, "GDC_2018-07-11-1200.json"),
                     theDataSetName="GDC Name in Dropdown",
                     theDataSetLabel="Choose data results",
                     theDefaultLink =c("PCA"),
                     theLevelLabels = c("Test Run Options", "Algorithm", "Diagram Type", "Sub-Type"),
                     theTooltipFile=system.file("BEVIndex", "tooltip_gdc.tsv", package="MBatchUtils"))
  TRUE
} else {
  message("No test data. Skip test.")
  TRUE
}
