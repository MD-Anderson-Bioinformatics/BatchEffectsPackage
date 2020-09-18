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

package edu.mda.bioinfo.bevindex.algorithms;

import edu.mda.bioinfo.bevindex.gson.DropdownEntry;
import edu.mda.bioinfo.bevindex.gson.MBatchIndexData;
import java.io.File;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author Tod-Casasent
 */
public class BoxPlot
{
	public File mAlgDir = null;
	public MBatchIndexData mMID = null;
	
	public BoxPlot(File theAlgDir, MBatchIndexData theMID)
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
			de.entry_label = "Boxplot";
			de.dropdown_label = "Diagram Type";
			// dropdown_entries for diagram types like All-Sample and Group-Mean
			File [] diaDirs = mAlgDir.listFiles();
			for (File diaDir : diaDirs)
			{
				DropdownEntry bpType = new DropdownEntry();
				bpType.entry_label = diaDir.getName();
				bpType.dropdown_label = "Batch Type";
				de.addDropdownEntry(bpType);
				// dropdown_entries for batch types
				TreeSet<String> btSet = new TreeSet<>();
				String filePattern = "BoxPlot_" + diaDir.getName() + "_Diagram-" + "(.*?)" + "\\.png";
				Pattern pattern = Pattern.compile(filePattern);
				File [] btFiles = diaDir.listFiles();
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
					btDe.diagram_type = "boxplot";
					btDe.box_annotations = "/BoxPlot/" + diaDir.getName() + "/BoxPlot_" + diaDir.getName() + "_Annotations-" + bt + ".tsv";
					btDe.box_data = "/BoxPlot/" + diaDir.getName() + "/BoxPlot_" + diaDir.getName() + "_BoxData-" + bt + ".tsv";
					btDe.box_histogram = "/BoxPlot/" + diaDir.getName() + "/BoxPlot_" + diaDir.getName() + "_Histogram-" + bt + ".tsv";
					btDe.diagram_image = "/BoxPlot/" + diaDir.getName() + "/BoxPlot_" + diaDir.getName() + "_Diagram-" + bt + ".png";
					bpType.addDropdownEntry(btDe);
					added = true;
				}
			}
			if (added)
			{
				mMID.addDropdownEntry(de);
			}
		}
	}
}
