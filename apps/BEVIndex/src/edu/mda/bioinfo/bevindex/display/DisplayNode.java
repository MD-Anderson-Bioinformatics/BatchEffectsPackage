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

package edu.mda.bioinfo.bevindex.display;

import edu.mda.bioinfo.bevindex.utils.BEVtooltips;
import java.nio.file.Path;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 *
 * @author linux
 */
public class DisplayNode extends DisplayRun
{
	// all member variables are declared here, so JSON output is same for this and descendants
	public int mLevel;
	public boolean mDiagramFlag;
	public String mZipFile;
	public String mInternalLocation;
	public String mAlgorithm;
	public TreeSet<String> mOtherFiles;
	transient protected String [] mTransientLevelLabels;
	transient protected Path mTransientLocation;

	public DisplayNode(Path theLocation, int theLevel, String [] theLevelLabels, BEVtooltips theTooltips, 
			String theArchiveMarker, String theArchiveName)
	{
		super(getLevelLabel(theLevel, theLevelLabels), convertFilenameToDisplay(theLocation.getFileName().toString()), null, theTooltips,  theArchiveMarker, theArchiveName);
		//super(getLevelLabel(theLevel, theLevelLabels), convertFilenameToDisplay(theLocation.getFileName().toString()), null, theTooltips,  theArchiveMarker, theArchiveName);
		mLevel = theLevel;		// set here for real
		mDiagramFlag = false;	// set false here or true/false in descendants
		mZipFile = "";			// set in descendants
		mInternalLocation = "";	// set in descendants
		mAlgorithm = "";		// set in descendants
		mOtherFiles = null;		// set in descendants
		mTransientLocation = theLocation;	// set here
		mTransientLevelLabels = theLevelLabels;		// set here
		//mNodeType = "DisplayNode";
	}

	public void init(BEVtooltips theTooltips) throws Exception
	{
		if (mTransientLocation.toFile().isDirectory())
		{
			mChildren = getDisplayNodes(mTransientLocation, mLevel+1, mTransientLevelLabels, theTooltips);
			for (DisplayNode dn : mChildren)
			{
				dn.init(theTooltips);
			}
		}
	}
	
	static protected String getLevelLabel(int theLevel, String [] theLevelLabels)
	{
		if (theLevel>=theLevelLabels.length)
		{
			return "Diagram";
		}
		return theLevelLabels[theLevel];
	}
	
	static protected String convertFilenameToDisplay(String theString)
	{
		//return theString.replaceAll("([A-Z,-])", " $1").replaceAll("([_])", "");
		//return new StringBuilder(new StringBuilder(theString).reverse().toString().replaceFirst("([_])", " ")).reverse().toString().replaceAll("([_])", "-").replace(".tsv", "").replace(".TSV", "").replace(".png", "").replace(".PNG", "");
		return theString.replace(".tsv", "").replace(".TSV", "").replace(".png", "").replace(".PNG", "");
	}

	public void dataRelations(TreeMap<String, String> theDataRelations, String thePath)
	{
		if ((null!=this.mZipFile)&&(!"".equals(this.mZipFile)))
		{
			theDataRelations.put(thePath + " - " + this.mName, this.mZipFile);
		}
		else
		{
			for (DisplayNode dn : this.mChildren)
			{
				dn.dataRelations(theDataRelations, thePath + " - " + this.mName);
			}
		}
	}
}
