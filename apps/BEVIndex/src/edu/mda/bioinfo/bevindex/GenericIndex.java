/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.mda.bioinfo.bevindex;

import static edu.mda.bioinfo.bevindex.BEVIndex.mVersion;
import static edu.mda.bioinfo.bevindex.GDCRunIndex.addReleaseFiles;
import java.io.File;

public class GenericIndex
{
	////////////////////////////////////////////////////////////////////////////
	//// main
	////////////////////////////////////////////////////////////////////////////
	public static void main(String[] args)
	{
		System.out.println(mVersion);
		System.out.println("Start GenericIndex");
		try
		{
			////////////////////////////////////////////////////////////////////
			// release sub-directory which will be for this data release (releases are sorted in this order)
			// this directory will form the first entry in the dropdowns to select this data
			String DEST_DATA_DIR=args[0];
			// label to use in GUI for this release
			String RELEASE_LABEL=args[1];
			// processing sub-directory which CURRENTLY contains this data release
			// This must match base name for DEST_DATA_DIR
			// this directory will form the first entry in the dropdowns to select this data
			String SOURCE_DATA_DIR=args[2];
			// name of ZIP archive to create
			String ZIP_NAME=args[3];
			// name of files containing data (Strongly suggest using matrix_data.tsv)
			String DATA_FILES=args[4];
			// comma delimited list of level labels
			// labels start below the DEST_DATA_DIR (or SOURCE_DATA_DIR). The DEST/SOURCE_DATA_DIR will be labeled "Data Run"
			String LEVEL_LABELS=args[5];
			// comma delimited list of directory names listing default data to show
			// labels start below the DEST_DATA_DIR (or SOURCE_DATA_DIR). The DEST/SOURCE_DATA_DIR will be labeled "Data Run"
			String DEFAULT_PATH=args[6];
			// path and file for JSON outout index (releases are sorted in this order)
			// This should match base name for DEST_DATA_DIR
			String INDEX_JSON=args[7];
			// tooltip file
			String TOOLTIP=null;
			if (!"null".equalsIgnoreCase(args[8]))
			{
				TOOLTIP=args[8];
			}
			System.out.println("DEST_DATA_DIR=" + DEST_DATA_DIR);
			System.out.println("RELEASE_LABEL=" + RELEASE_LABEL);
			System.out.println("SOURCE_DATA_DIR=" + SOURCE_DATA_DIR);
			System.out.println("ZIP_NAME=" + ZIP_NAME);
			System.out.println("DATA_FILES=" + DATA_FILES);
			System.out.println("LEVEL_LABELS=" + LEVEL_LABELS);
			System.out.println("DEFAULT_PATH=" + DEFAULT_PATH);
			System.out.println("INDEX_JSON=" + INDEX_JSON);
			System.out.println("TOOLTIP=" + TOOLTIP);
			////////////////////////////////////////////////////////////////////
			System.out.println("////////////////////////////////////////////////////////////////////");
			System.out.println("BEVRunIndex.cleanOldZips");
			BEVRunIndex.cleanOldZips(SOURCE_DATA_DIR, ZIP_NAME);
			System.out.println("////////////////////////////////////////////////////////////////////");
			System.out.println("addReleaseFiles");
			addReleaseFiles(SOURCE_DATA_DIR, "RunInfo.tsv", "RunInfo.PNG", RELEASE_LABEL, DATA_FILES);
			String [] levelLabels = LEVEL_LABELS.split(",", -1);
			String [] defaultDisplay = DEFAULT_PATH.split(",", -1);
			// addIndex(String baseDir, String runDir, String runName, String tooltipTsv, String theArchiveMarker, String theArchiveName, 
			//          String [] levelLabels, String [] defaultDisplay) throws IOException
			System.out.println("////////////////////////////////////////////////////////////////////");
			System.out.println("BEVRunIndex.addIndex");
			BEVRunIndex.addIndex(SOURCE_DATA_DIR, SOURCE_DATA_DIR, RELEASE_LABEL, TOOLTIP, DATA_FILES, ZIP_NAME, levelLabels, defaultDisplay);
			System.out.println("////////////////////////////////////////////////////////////////////");
			System.out.println("BEVRunIndex.buildZipFiles");
			BEVRunIndex.buildZipFiles(SOURCE_DATA_DIR, SOURCE_DATA_DIR, DATA_FILES, ZIP_NAME);
			////////////////////////////////////////////////////////////////////
			System.out.println("////////////////////////////////////////////////////////////////////");
			System.out.println("BEVRunIndex.addOverallIndex");
			String indexFile = new File(INDEX_JSON).getAbsolutePath();
			// addOverallIndex(String baseDir, String runDir, String runName, String tooltipTsv,
			//   String indexFile, String finalDir, String theArchiveMarker, String theArchiveName,
			//   String [] levelLabels, String [] defaultDisplay
			BEVRunIndex.addOverallIndex(SOURCE_DATA_DIR, SOURCE_DATA_DIR, RELEASE_LABEL, TOOLTIP,
					indexFile, DEST_DATA_DIR, DATA_FILES, ZIP_NAME, 
					levelLabels, defaultDisplay, true);
			System.out.println("////////////////////////////////////////////////////////////////////");
		}
		catch (Exception exp)
		{
			exp.printStackTrace(System.err);
		}
		System.out.println("Done GenericIndex");
	}

}
