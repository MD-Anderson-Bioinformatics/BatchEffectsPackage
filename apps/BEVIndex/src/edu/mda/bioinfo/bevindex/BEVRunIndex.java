/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.mda.bioinfo.bevindex;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import static edu.mda.bioinfo.bevindex.BEVIndex.mVersion;
import edu.mda.bioinfo.bevindex.display.DisplayRun;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.TreeSet;

/**
 *
 * @author linux
 */
public class BEVRunIndex
{

	////////////////////////////////////////////////////////////////////////////
	//// built overall run JSON index
	////////////////////////////////////////////////////////////////////////////
	public static void addOverallIndex(String baseDir, String runDir,
			String runName, String tooltipTsv, String indexFile, String finalDir, String theArchiveMarker, String theArchiveName,
			String [] levelLabels, String [] defaultDisplay, boolean useStdFlag) throws IOException, Exception
	{
		System.out.println("baseDir="+baseDir);
		System.out.println("runDir="+runDir);
		System.out.println("runName="+runName);
		System.out.println("tooltipTsv="+tooltipTsv);
		System.out.println("indexFile="+indexFile);
		System.out.println("finalDir="+finalDir);
		System.out.println("theArchiveMarker="+theArchiveMarker);
		System.out.println("theArchiveName="+theArchiveName);
		for (String level : levelLabels)
		{
			System.out.println("levelLabels="+level);
		}
		BEVtooltips btt = new BEVtooltips();
		btt.readTsv(tooltipTsv);
		System.out.println("build run object");
		DisplayRun dr = new DisplayRun(theArchiveMarker, theArchiveName);
		if (false==useStdFlag)
		{
			dr.initStart(Paths.get(runDir), runName, levelLabels, "Data Run", defaultDisplay, btt);
		}
		else
		{
			dr.init(Paths.get(runDir), runName, levelLabels, "Data Run", defaultDisplay, btt);
		}
		System.out.println("update paths");
		dr.updateZipPaths(baseDir, finalDir);
		//
		System.out.println("write to indexFile=" + indexFile);
		GsonBuilder builder = new GsonBuilder();
		builder.setPrettyPrinting();
		Gson gson = builder.create();
		File output = new File(indexFile);
		if (output.exists())
		{
			output.delete();
		}
		Files.write(output.toPath(), gson.toJson(dr).getBytes());
	}

	////////////////////////////////////////////////////////////////////////////
	//// add JSON index to existing ZIP files
	////////////////////////////////////////////////////////////////////////////
	public static void addIndexToDataSet(String theBaseDir, Path theMarkerFilepath, String[] theDefault,
			String theRunDisplayName, String theRunLabel, String[] theLevelLabels, String theTooltipTsv, String theArchiveMarker, String theArchiveName) throws Exception
	{
		System.out.println(mVersion);
		System.out.println("theBaseDir=" + theBaseDir);
		System.out.println("theMarkerFilepath=" + theMarkerFilepath);
		System.out.println("theDefault=" + Arrays.toString(theDefault));
		System.out.println("theRunDisplayName=" + theRunDisplayName);
		System.out.println("theRunLabel=" + theRunLabel);
		System.out.println("theLevelLabels=" + Arrays.toString(theLevelLabels));
		System.out.println("theTooltipTsv=" + theTooltipTsv);

		BEVtooltips btt = new BEVtooltips();
		btt.readTsv(theTooltipTsv);
		DisplayRun dr = new DisplayRun(theArchiveMarker, theArchiveName);
		dr.init(Paths.get(theBaseDir), theRunDisplayName, theLevelLabels, theRunLabel,
				Arrays.copyOf(theDefault, theDefault.length), theMarkerFilepath.getParent(), btt);
		dr.removeZipPaths();
		//
		GsonBuilder builder = new GsonBuilder();
		builder.setPrettyPrinting();
		Gson gson = builder.create();
		File newJson = new File(theMarkerFilepath.getParent().toFile(), "index.json");
		try (BufferedWriter bw = java.nio.file.Files.newBufferedWriter(Paths.get(newJson.getAbsolutePath()), Charset.availableCharsets().get("UTF-8")))
		{
			bw.write(gson.toJson(dr));
		}
	}

	public static void addIndex(String baseDir, String runDir, String runName, String tooltipTsv, String theArchiveMarker, String theArchiveName,
			String [] levelLabels, String [] defaultDisplay) throws IOException
	{
		System.out.println("Add index to " + runDir);
		FileFind ff = new FileFind();
		ff.find(Paths.get(runDir), theArchiveMarker);
		TreeSet<Path> toZipSet = ff.mMatches;
		List<Exception> errors = Collections.synchronizedList(new ArrayList<>());
		//toZipSet.stream()
		System.out.println("PARALLEL-THREADS = 5");
		System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", "5");
		toZipSet.parallelStream()
				.forEach((Path marker) ->
				{
					System.out.println("------------------------------------------");
					System.out.println("marker=" + marker);
					try
					{
						addIndexToDataSet(baseDir, marker, defaultDisplay, runName, "Data Run", levelLabels, tooltipTsv, theArchiveMarker, theArchiveName);
					}
					catch (Exception exp)
					{
						errors.add(new Exception("Error adding index to marker directory=" + marker, exp));
					}
				});
		for (Exception exp : errors)
		{
			exp.printStackTrace(System.err);
			exp.getCause().printStackTrace(System.err);
		}
		System.out.println("------------------------------------------");
	}

	////////////////////////////////////////////////////////////////////////////
	//// ZIP directories with successful MBatch runs
	////////////////////////////////////////////////////////////////////////////
	public static void buildZipFiles(String runDir, String baseDir, String theArchiveMarker, String theArchiveName) throws IOException
	{
		FileFind ff = new FileFind();
		ff.find(Paths.get(runDir), theArchiveMarker);
		TreeSet<Path> toZipSet = ff.mMatches;
		// note this solution changes the parallelism settings for the entire system
		// so if you have multiple parallel operations going on this could do weird things
		List<Exception> errors = Collections.synchronizedList(new ArrayList<>());
		//System.out.println("PARALLEL-THREADS = 5");
		//System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", "5");
		//toZipSet.parallelStream()
		toZipSet.stream()
				.forEach((Path archive) ->
				{
					Path archDir = archive.getParent();
					Path zipFile = new File(archDir.toFile(), theArchiveName).toPath();
					System.out.println("------------------------------------------");
					System.out.println("archDir=" + archDir);
					System.out.println("zipFile=" + zipFile);
					ZipDir zd = new ZipDir();
					try
					{
						zd.zipContents(archDir, zipFile);
					}
					catch (IOException exp)
					{
						errors.add(new Exception("Error zipping archDir=" + archDir + " to zipFile=" + zipFile, exp));
					}
				});
		for (Exception exp : errors)
		{
			exp.printStackTrace(System.err);
			exp.getCause().printStackTrace(System.err);
		}
		System.out.println("------------------------------------------");
	}

	////////////////////////////////////////////////////////////////////////////
	//// clean old zips from inclusion
	////////////////////////////////////////////////////////////////////////////
	public static void cleanOldZips(String runDir, String archiveName) throws IOException
	{
		System.out.println("Find and clear " + archiveName + " from " + runDir);
		FileFind ff = new FileFind();
		ff.find(Paths.get(runDir), archiveName);
		TreeSet<Path> oldZipSet = ff.mMatches;
		System.out.println("removing");
		oldZipSet.stream()
				.forEach((Path oldZip) ->
				{
					oldZip.toFile().delete();
				});
	}

	////////////////////////////////////////////////////////////////////////////
	//// main
	////////////////////////////////////////////////////////////////////////////
	public static void main(String[] args)
	{
		System.out.println(mVersion);
		try
		{
			////////////////////////////////////////////////////////////////////
			String baseDir = args[0];
			String runDir = args[1];
			String runName = args[2];
			String indexDir = args[3];
			String finalDir = args[4];
			String tooltipTsv = args[5];
			System.out.println("baseDir=" + baseDir);
			System.out.println("runDir=" + runDir);
			System.out.println("runName=" + runName);
			System.out.println("indexDir=" + indexDir);
			System.out.println("finalDir=" + finalDir);
			System.out.println("tooltipTsv=" + tooltipTsv);
			////////////////////////////////////////////////////////////////////
			String archiveMarker = "MBATCH_SUCCESS.txt";
			String archiveName = "ResultSet.zip";
			cleanOldZips(runDir, archiveName);
			String [] levelLabels = BEVIndex.mBEVlabels;
			String [] defaultDisplay = BEVIndex.mBEVcurrentDefault;
			if (runName.contains("DCC"))
			{
				defaultDisplay = BEVIndex.mBEVlegacyDefault;
			}
			addIndex(baseDir, runDir, runName, tooltipTsv, archiveMarker, archiveName, levelLabels, defaultDisplay);
			buildZipFiles(runDir, baseDir, archiveMarker, archiveName);
			////////////////////////////////////////////////////////////////////
			String indexFile = new File(indexDir, runName + ".json").getAbsolutePath();
			addOverallIndex(baseDir, runDir, runName, tooltipTsv,
					indexFile, finalDir, archiveMarker, archiveName, levelLabels, defaultDisplay, false);
		}
		catch (Exception exp)
		{
			exp.printStackTrace(System.err);
		}
	}
}
