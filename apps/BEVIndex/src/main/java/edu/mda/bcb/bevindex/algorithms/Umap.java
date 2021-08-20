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

/**
 *
 * @author Tod-Casasent
 */
public class Umap
{
	public File mAlgDir = null;
	public MBatchIndexData mMID = null;
	
	public Umap(File theAlgDir, MBatchIndexData theMID)
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
			de.entry_label = "UMAP";
			de.dropdown_label = "Batch Type";
			// dropdown_entries for batch types (batch id, ship date, etc)
			File [] diaDirs = mAlgDir.listFiles();
			for (File diaDir : diaDirs)
			{
				DropdownEntry batchType = new DropdownEntry();
				batchType.entry_label = diaDir.getName();
				batchType.diagram_type = "umap";
				batchType.umap_batches = "/UMAP/" + diaDir.getName() + "/UMAP_Data-batc.tsv";
				batchType.umap_samples = "/UMAP/" + diaDir.getName() + "/UMAP_Data-samp.tsv";
				batchType.umap_umapdat = "/UMAP/" + diaDir.getName() + "/UMAP_Data-umap.tsv";
				File imageFile = new File(diaDir, "UMAP_Diagram.png");
				File legendFile = new File(diaDir, "UMAP_Legend.png");
				if (imageFile.exists())
				{
					batchType.diagram_image = "/UMAP/" + diaDir.getName() + "/UMAP_Diagram.png";
				}
				if (legendFile.exists())
				{
					batchType.legend_image = "/UMAP/" + diaDir.getName() + "/UMAP_Legend.png";
				}
				de.addDropdownEntry(batchType);
				added = true;
			}
			if (added)
			{
				mMID.addDropdownEntry(de);
			}
		}
	}
}
