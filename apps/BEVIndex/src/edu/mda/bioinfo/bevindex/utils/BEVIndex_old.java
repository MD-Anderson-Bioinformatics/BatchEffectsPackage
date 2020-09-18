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
import edu.mda.bioinfo.bevindex.display.DisplayRun;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;

/**
 *
 * @author Tod-Casasent
 */
public class BEVIndex_old
{

	static protected String mVersion = "BEVIndex 2020-09-11-1000";
	static public String [] mBEVlabels = { "Version", "Program", "Disease", "Workflow", "Data Type", "Algorithm", "Diagram Type", "Sub-Type" };
	static public String [] mBEVcurrentDefault = { "current", "TCGA", "KIRC", "methylation", "All-original", "PCA", "BatchId", "ManyToMany", "PCAValues" };
	static public String [] mBEVlegacyDefault = { "legacy", "TCGA", "KIRC", "methylation", "All-original", "PCA", "BatchId", "ManyToMany", "PCAValues" };
	static public String [] mSDBlabels = { "Version", "Status", "Program", "Disease", "Workflow" };
	static public String [] mSDBcurrentDefault = { "current", "standardized", "TCGA", "KIRC", "methylation" };
	static public String [] mSDBlegacyDefault = { "legacy", "standardized", "TCGA", "KIRC", "methylation" };
	
	public static void addIndexAndZipDataSet(String theDataRunDir, String theZipFilepath, String[] theDefault,
			String theRunDisplayName, String theRunLabel, String[] theLevelLabels, String theTooltipTsv,
			String theNoticeText) throws Exception
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
		DisplayRun dr = new DisplayRun(archiveMarker, archiveName, theNoticeText);
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
