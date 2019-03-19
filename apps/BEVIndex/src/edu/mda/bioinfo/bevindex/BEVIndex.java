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
import java.net.URI;
import java.nio.file.FileSystem;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author linux
 */
public class BEVIndex
{

	static protected String mVersion = "BEVIndex 2019-03-07-1313";
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
		Map<String, String> env = new HashMap<>();
		env.put("create", "true");
		Path path = Paths.get(theZipFilepath);
		URI uri = URI.create("jar:" + path.toUri());
		System.out.println("uri=" + uri);
		try (FileSystem fs = FileSystems.newFileSystem(uri, env))
		{
			Path nf = fs.getPath("index.json");
			Files.write(nf, gson.toJson(dr).getBytes(), StandardOpenOption.WRITE, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.CREATE);
		}
	}

	public static void addExternalIndexToDataRun(String theDataRunDir, String[] theDataRuns, String[] theDefault,
			String[] theNameRuns, String[] theLabelRuns, String[] theFilesRuns, String[] theLevelLabels, String theTooltipTsv) throws Exception
	{
		System.out.println(mVersion);
		BEVtooltips btt = new BEVtooltips();
		btt.readTsv(theTooltipTsv);
		for (int index = 0; index < theDataRuns.length; index++)
		{
			String archiveMarker = "MBATCH_SUCCESS.txt";
			String archiveName = "ResultSet.zip";
			DisplayRun dr = new DisplayRun(archiveMarker, archiveName);
			//(Path theBaseLocation, String theName, boolean theGDCFlag, String theLabel, String [] theDefaultDiagram)
			dr.init(Paths.get(theDataRunDir), theNameRuns[index], theLevelLabels, theLabelRuns[index], theDefault, Paths.get(theDataRuns[index]), btt);
			//
			GsonBuilder builder = new GsonBuilder();
			builder.setPrettyPrinting();
			Gson gson = builder.create();
			File output = new File(theFilesRuns[index]);
			Files.write(output.toPath(), gson.toJson(dr).getBytes());
		}
	}
}
