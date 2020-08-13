// Copyright (c) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020 University of Texas MD Anderson Cancer Center
//
// This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
// MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>

package edu.mda.bioinfo.bevindex.utils;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import static edu.mda.bioinfo.bevindex.utils.BEVIndex_old.mVersion;
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
			String [] levelLabels, String [] defaultDisplay, boolean useStdFlag,
			String theNoticeText) throws IOException, Exception
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
		DisplayRun dr = new DisplayRun(theArchiveMarker, theArchiveName, theNoticeText);
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
			String theRunDisplayName, String theRunLabel, String[] theLevelLabels, String theTooltipTsv, 
			String theArchiveMarker, String theArchiveName, String theArchiveNotice) throws Exception
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
		DisplayRun dr = new DisplayRun(theArchiveMarker, theArchiveName, theArchiveNotice);
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
			String [] levelLabels, String [] defaultDisplay, String theArchiveNotice) throws IOException
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
						addIndexToDataSet(baseDir, marker, defaultDisplay, runName, "Data Run", levelLabels, tooltipTsv, theArchiveMarker, theArchiveName, theArchiveNotice);
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
			String noticeFile = args[6];
			System.out.println("baseDir=" + baseDir);
			System.out.println("runDir=" + runDir);
			System.out.println("runName=" + runName);
			System.out.println("indexDir=" + indexDir);
			System.out.println("finalDir=" + finalDir);
			System.out.println("tooltipTsv=" + tooltipTsv);
			System.out.println("noticeFile=" + noticeFile);
			if ("null".equalsIgnoreCase(noticeFile))
			{
				noticeFile = null;
			}
			System.out.println("noticeFile2=" + noticeFile);
			////////////////////////////////////////////////////////////////////
			String noticeString = null;
			if (null!=noticeFile)
			{
				noticeString = new String(Files.readAllBytes(Paths.get(noticeFile)));
			}
			System.out.println("noticeString=" + noticeString);
			////////////////////////////////////////////////////////////////////
			String archiveMarker = "MBATCH_SUCCESS.txt";
			String archiveName = "ResultSet.zip";
			cleanOldZips(runDir, archiveName);
			String [] levelLabels = BEVIndex_old.mBEVlabels;
			String [] defaultDisplay = BEVIndex_old.mBEVcurrentDefault;
			if (runName.contains("DCC"))
			{
				defaultDisplay = BEVIndex_old.mBEVlegacyDefault;
			}
			addIndex(baseDir, runDir, runName, tooltipTsv, archiveMarker, archiveName, levelLabels, defaultDisplay, noticeString);
			buildZipFiles(runDir, baseDir, archiveMarker, archiveName);
			////////////////////////////////////////////////////////////////////
			String indexFile = new File(indexDir, runName + ".json").getAbsolutePath();
			addOverallIndex(baseDir, runDir, runName, tooltipTsv,
					indexFile, finalDir, archiveMarker, archiveName, levelLabels, defaultDisplay, false, noticeString);
		}
		catch (Exception exp)
		{
			exp.printStackTrace(System.err);
		}
	}
}
