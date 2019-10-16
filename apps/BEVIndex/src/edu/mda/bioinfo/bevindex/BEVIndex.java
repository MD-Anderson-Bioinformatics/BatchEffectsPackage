/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.mda.bioinfo.bevindex;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import edu.mda.bioinfo.bevindex.display.DisplayRun;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;

/**
 *
 * @author linux
 */
public class BEVIndex
{

	static protected String mVersion = "BEVIndex 2019-09-04-1400";
	static public String [] mBEVlabels = { "Version", "Program", "Disease", "Workflow", "Data Type", "Algorithm", "Diagram Type", "Sub-Type" };
	static public String [] mBEVcurrentDefault = { "current", "TCGA", "KIRC", "methylation", "All-original", "PCA", "BatchId", "ManyToMany", "PCAValues" };
	static public String [] mBEVlegacyDefault = { "legacy", "TCGA", "KIRC", "methylation", "All-original", "PCA", "BatchId", "ManyToMany", "PCAValues" };
	static public String [] mSDBlabels = { "Version", "Status", "Program", "Disease", "Workflow" };
	static public String [] mSDBcurrentDefault = { "current", "standardized", "TCGA", "KIRC", "methylation" };
	static public String [] mSDBlegacyDefault = { "legacy", "standardized", "TCGA", "KIRC", "methylation" };
	
	// TODO: update BENIndex for R use
	public static void addIndexAndZipDataSet(String theDataRunDir, String theZipFilepath, String[] theDefault,
			String theRunDisplayName, String theRunLabel, String[] theLevelLabels, String theTooltipTsv) throws Exception
	{
		System.out.println(mVersion);
		System.out.println("theDataRunDir=" + theDataRunDir);
		System.out.println("theZipFilepath=" + theZipFilepath);
		System.out.println("theDefault=" + Arrays.toString(theDefault));
		System.out.println("theRunDisplayName=" + theRunDisplayName);
		System.out.println("theRunLabel=" + theRunLabel);
		System.out.println("theLevelLabels=" + Arrays.toString(theLevelLabels));
		System.out.println("theTooltipTsv=" + theTooltipTsv);

		BEVtooltips btt = new BEVtooltips();
		btt.readTsv(theTooltipTsv);
		String archiveMarker = "MBATCH_SUCCESS.txt";
		String archiveName = "ResultSet.zip";
		DisplayRun dr = new DisplayRun(archiveMarker, archiveName);
		dr.init(Paths.get(theDataRunDir), theRunDisplayName, theLevelLabels, theRunLabel,
				Arrays.copyOf(theDefault, theDefault.length), Paths.get(theZipFilepath).getParent(), btt);
		//
		GsonBuilder builder = new GsonBuilder();
		builder.setPrettyPrinting();
		Gson gson = builder.create();
		// Write file out (pre-ZIP), adding to existing ZIP corrupts file
		File indexFile = new File(new File(theZipFilepath).getParentFile(), "index.json");
		Path indexPath = indexFile.toPath();
		Files.write(indexPath, gson.toJson(dr).getBytes(), StandardOpenOption.WRITE, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.CREATE);
	}
}
