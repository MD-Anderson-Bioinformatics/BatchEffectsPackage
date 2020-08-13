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

/**
 *
 * @author linux
 */
public class DSC
{
	public File mAlgDir = null;
	public MBatchIndexData mMID = null;
	
	public DSC(File theAlgDir, MBatchIndexData theMID)
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
			de.entry_label = "Dispersion Separability Criterion";
			de.diagram_type = "DSC";
			de.dsc_values = "DSCOverview.tsv";
			if (new File(mAlgDir, "DSCOverview.tsv").exists())
			{
				added = true;
			}
			if (added)
			{
				mMID.addDropdownEntry(de);
			}
		}
	}
}
