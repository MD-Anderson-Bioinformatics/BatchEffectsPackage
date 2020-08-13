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

import static edu.mda.bioinfo.bevindex.utils.BEVIndex_old.mVersion;
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
		if ("matrix_data.tsv".equals(archiveLocation.toFile().getName()))
		{
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
			// complete path to release directory (timestamp)
			String baseDir = args[0];
			// complete path to current or legacy directory under release
			String runDir = args[1];
			// run name of form 2018_11_20_1200_GDC (current)
			String runName = args[2];
			// directory into which to write index file
			String indexDir = args[3];
			// filname directory (usually SDB)
			String finalDir = args[4];
			System.out.println("baseDir=" + baseDir);
			System.out.println("runDir=" + runDir);
			System.out.println("runName=" + runName);
			System.out.println("indexDir=" + indexDir);
			System.out.println("finalDir=" + finalDir);
			//////////////////////////////////////////////////////////////////
			String archiveMarker = "matrix_data.tsv";
			String archiveName = "standardized.zip";
			BEVRunIndex.cleanOldZips(runDir, archiveName);
			addReleaseFiles(runDir, "RunInfo.tsv", "RunInfo.PNG", runName, archiveMarker);
			String [] levelLabels = BEVIndex_old.mSDBlabels;
			String [] defaultDisplay = BEVIndex_old.mSDBcurrentDefault;
			if (runName.contains("DCC"))
			{
				defaultDisplay = BEVIndex_old.mSDBlegacyDefault;
			}
			BEVRunIndex.addIndex(baseDir, runDir, runName, null, archiveMarker, archiveName, levelLabels, defaultDisplay, null);
			BEVRunIndex.buildZipFiles(runDir, baseDir, archiveMarker, archiveName);
			////////////////////////////////////////////////////////////////////
			String indexFile = new File(indexDir, runName + ".json").getAbsolutePath();
			BEVRunIndex.addOverallIndex(baseDir, runDir, runName, null,
					indexFile, finalDir, archiveMarker, archiveName, levelLabels, defaultDisplay, false, null);
		}
		catch (Exception exp)
		{
			exp.printStackTrace(System.err);
		}
		System.out.println("Done GDCRunIndex");
	}
}
