/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.mda.bioinfo.bevindex.display;

import edu.mda.bioinfo.bevindex.BEVtooltips;
import java.nio.file.Path;
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
		return new StringBuilder(new StringBuilder(theString).reverse().toString().replaceFirst("([_])", " ")).reverse().toString().replaceAll("([_])", "-").replace(".tsv", "").replace(".TSV", "").replace(".png", "").replace(".PNG", "");
	}

}
