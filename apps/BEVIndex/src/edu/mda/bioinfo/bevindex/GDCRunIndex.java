/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.mda.bioinfo.bevindex;

import static edu.mda.bioinfo.bevindex.BEVIndex.mVersion;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeSet;
import javax.imageio.ImageIO;

/**
 *
 * @author linux
 */
public class GDCRunIndex
{

	////////////////////////////////////////////////////////////////////////////
	//// add release name to existing ZIP files
	////////////////////////////////////////////////////////////////////////////
	private static LinkedHashMap<String, String> getReleaseData(String runName, Path archiveLocation) throws IOException
	{
		LinkedHashMap<String, String> results = new LinkedHashMap<>();
		results.put("Run", runName);
		try (BufferedReader br = new BufferedReader(new FileReader(archiveLocation.toFile())))
		{
			int samples = br.readLine().split("\t", -1).length - 1;
			int features = 0;
			while ((br.readLine()) != null)
			{
				features += 1;
			}
			results.put("Samples", Integer.toString(samples));
			results.put("Features", Integer.toString(features));
		}
		return results;
	}
	
	// addReleaseFiles(runDir, "RunInfo.tsv", "RunInfo.PNG", archiveMarker);
	public static void addReleaseFiles(String runDir, String runFile, String imageFile, String runName, String archiveMarker) throws IOException
	{
		System.out.println("************************* addReleaseToZipFiles");
		FileFind ff = new FileFind();
		ff.find(Paths.get(runDir), archiveMarker);
		TreeSet<Path> locationSet = ff.mMatches;
		List<Exception> errors = Collections.synchronizedList(new ArrayList<>());
		System.out.println("PARALLEL-THREADS = 5");
		System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", "5");
		locationSet.parallelStream()
				.forEach((Path archiveLoc) ->
				{
					System.out.println("------------------------------------------");
					System.out.println("archiveLoc=" + archiveLoc);
					try
					{
						LinkedHashMap<String, String> release = getReleaseData(runName, archiveLoc);
						// write to file
						File tsvFile = new File(archiveLoc.getParent().toFile(), runFile);
						File pngFile = new File(archiveLoc.getParent().toFile(), imageFile);
						try(BufferedWriter bw = java.nio.file.Files.newBufferedWriter(Paths.get(tsvFile.getAbsolutePath()), Charset.availableCharsets().get("UTF-8")))
						{
							for (Entry<String, String> lineEntry : release.entrySet())
							{
								bw.write((lineEntry.getKey() + "\t" + lineEntry.getValue() + "\n"));
							}
						}
						//
						BufferedImage bufferedImage = new BufferedImage(400, 200, BufferedImage.TYPE_INT_RGB);
						Graphics g = bufferedImage.getGraphics();
						g.setColor(Color.WHITE);
						g.fillRect(0, 0, bufferedImage.getWidth(), bufferedImage.getHeight());
						g.setColor(Color.BLACK);
						int row = 10;
						for (Entry<String, String> lineEntry : release.entrySet())
						{
							g.drawString((lineEntry.getKey() + ": " + lineEntry.getValue()), 10, row);
							row += 11;
						}
						ImageIO.write(bufferedImage, "png", pngFile);
					}
					catch (Exception exp)
					{
						errors.add(new Exception("Error release file to archiveLoc=" + archiveLoc, exp));
						System.out.println("Error release file to archiveLoc=" + archiveLoc + " error: " + exp.getMessage());
						exp.printStackTrace(System.err);
						exp.getCause().printStackTrace(System.err);
					}
					System.out.flush();
				});
		for (Exception exp : errors)
		{
			exp.printStackTrace(System.err);
			exp.getCause().printStackTrace(System.err);
		}
		System.out.println("------------------------------------------");
	}

	////////////////////////////////////////////////////////////////////////////
	//// main
	////////////////////////////////////////////////////////////////////////////
	public static void main(String[] args)
	{
		System.out.println(mVersion);
		System.out.println("Start GDCRunIndex");
		try
		{
			////////////////////////////////////////////////////////////////////
			String baseDir = args[0];
			String runDir = args[1];
			String runName = args[2];
			String indexDir = args[3];
			String finalDir = args[4];
			//////////////////////////////////////////////////////////////////
			String archiveMarker = "matrix_data.tsv";
			String archiveName = "standardized.zip";
			BEVRunIndex.cleanOldZips(runDir, archiveName);
			addReleaseFiles(runDir, "RunInfo.tsv", "RunInfo.PNG", runName, archiveMarker);
			String [] levelLabels = BEVIndex.mSDBlabels;
			String [] defaultDisplay = BEVIndex.mSDBcurrentDefault;
			if (runName.contains("DCC"))
			{
				defaultDisplay = BEVIndex.mSDBlegacyDefault;
			}
			BEVRunIndex.addIndex(baseDir, runDir, runName, null, archiveMarker, archiveName, levelLabels, defaultDisplay);
			BEVRunIndex.buildZipFiles(runDir, baseDir, archiveMarker, archiveName);
			////////////////////////////////////////////////////////////////////
			String indexFile = new File(indexDir, runName + ".json").getAbsolutePath();
			BEVRunIndex.addOverallIndex(baseDir, runDir, runName, null,
					indexFile, finalDir, archiveMarker, archiveName, levelLabels, defaultDisplay, false);
		}
		catch (Exception exp)
		{
			exp.printStackTrace(System.err);
		}
		System.out.println("Done GDCRunIndex");
	}
}
