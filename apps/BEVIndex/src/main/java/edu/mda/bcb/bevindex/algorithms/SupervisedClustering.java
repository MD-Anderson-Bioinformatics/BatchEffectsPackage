// Copyright (c) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021 University of Texas MD Anderson Cancer Center
//
// This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
// MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>

package edu.mda.bcb.bevindex.algorithms;

import edu.mda.bcb.bevindex.gson.DropdownEntry;
import edu.mda.bcb.bevindex.gson.MBatchIndexData;
import java.io.File;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author Tod-Casasent
 */
public class SupervisedClustering
{
	public File mAlgDir = null;
	public MBatchIndexData mMID = null;
	
	public SupervisedClustering(File theAlgDir, MBatchIndexData theMID)
	{
		mAlgDir = theAlgDir;
		mMID = theMID;
	}
	
	public void checkAndAddAlgorithm()
	{
		if (mAlgDir.exists())
		{
			boolean added = false;
			DropdownEntry de = new DropdownEntry();
			de.entry_label = "Supervised Clustering";
			de.dropdown_label = "Batch Type";
			// iterate over batch type diagrams
			TreeSet<String> btSet = new TreeSet<>();
			String filePattern = "(.*?)" + "_ngchm\\.ngchm\\.html";
			Pattern pattern = Pattern.compile(filePattern);
			File batchDir = new File(mAlgDir, "Batches");
			File [] btFiles = batchDir.listFiles();
			for (File btFile : btFiles)
			{
				// BoxPlot_AllSample-Data_Diagram-BatchId.png
				Matcher matcher = pattern.matcher(btFile.getName());
				if (matcher.find())
				{
					btSet.add(matcher.group(1));
				}
			}
			for (String bt : btSet)
			{
				DropdownEntry btDe = new DropdownEntry();
				btDe.entry_label = bt;
				btDe.diagram_type = "ngchm";
				btDe.ngchm = "/SupervisedClustering/Batches/" + bt + "_ngchm.ngchm.html";
				btDe.diagram_image = "/SupervisedClustering/Batches/SupervisedClust_Diagram-" + bt + ".png";
				btDe.legend_image = "/SupervisedClustering/Batches/SupervisedClust_Legend-" + bt + ".png";
				de.addDropdownEntry(btDe);
				added = true;
			}
			if (added)
			{
				mMID.addDropdownEntry(de);
			}
		}
	}
}
