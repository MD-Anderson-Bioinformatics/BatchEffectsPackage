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

package edu.mda.bcb.bevindex.display;

import edu.mda.bcb.bevindex.utils.BEVtooltips;
import java.io.File;
import java.io.IOException;
import java.nio.file.FileVisitOption;
import java.nio.file.FileVisitResult;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.SimpleFileVisitor;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.EnumSet;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 *
 * @author Tod-Casasent
 */
public class DisplayRun implements Comparable<DisplayRun>
{
	//public String mNodeType = "DisplayRun";
	public String mLabel;
	public String mName;
	public String [] mDefaultDiagram = null;
	public TreeSet<DisplayNode> mChildren;
	public String mTooltip;
	public String mNotice;
	transient public String mArchiveMarker;
	transient public String mArchiveName;
	
	public DisplayRun(String theArchiveMarker, String theArchiveName, String theNotice)
	{
		mLabel = null;
		mName = null;
		mDefaultDiagram = null;
		mChildren = null;
		mArchiveMarker = theArchiveMarker;
		mArchiveName = theArchiveName;
		mNotice = theNotice;
	}
	
	public DisplayRun(String theLabel, String theName, TreeSet<DisplayNode> theChildren, BEVtooltips theTooltips, 
			String theArchiveMarker, String theArchiveName)
	{
		mLabel = theLabel;
		mName = theName; //DisplayNode.convertFilenameToDisplay(theName);
		mTooltip = theTooltips.getTooltip(theLabel, theName);
		mDefaultDiagram = null;
		mChildren = theChildren;
		mArchiveMarker = theArchiveMarker;
		mArchiveName = theArchiveName;
		mNotice = null;
	}
	
	public void init(Path theBaseLocation, String theName, String [] theLevelLabels, String theLabel, String [] theDefaultDiagram, 
			Path theDatasetDir, BEVtooltips theTooltips) throws Exception
	{
		// TODO: write dataset init
		System.out.println("Display Run (dataset) " + theBaseLocation);
		System.out.println("dataset " + theDatasetDir);
		mLabel = theLabel;
		mName = theName;
		mTooltip = theTooltips.getTooltip(theLabel, theName);
		mDefaultDiagram = theDefaultDiagram;
		System.out.println("Find Children " + theBaseLocation);
		// theCurrentPath, theLevel, theGDCFlag, null
		mChildren = new TreeSet<>();
		getDisplayNodesWithInit(theBaseLocation, theLevelLabels, theDatasetDir, theTooltips);
	}
	
	public void initStart(Path theSearchStart, String theName, String [] theLevelLabels, String theLabel, 
			String [] theDefaultDiagram, BEVtooltips theTooltips) throws Exception
	{
		System.out.println("Display Run (entire) theSearchStart=" + theSearchStart);
		mLabel = theLabel;
		mName = theName; // DisplayNode.convertFilenameToDisplay(theName);
		mTooltip = theTooltips.getTooltip(theLabel, theName);
		mDefaultDiagram = theDefaultDiagram;
		System.out.println("Start from " + theSearchStart);
		mChildren = new TreeSet<>();
		DisplayNode newChild = new DisplayNode(theSearchStart, 0, theLevelLabels, theTooltips, mArchiveMarker, mArchiveName);
		mChildren.add(newChild);
		newChild.init(theTooltips);
	}
	
	public void init(Path theSearchStart, String theName, String [] theLevelLabels, String theLabel, 
			String [] theDefaultDiagram, BEVtooltips theTooltips) throws Exception
	{
		System.out.println("Display Run (entire) theSearchStart=" + theSearchStart);
		mLabel = theLabel;
		mName = theName;
		mTooltip = theTooltips.getTooltip(theLabel, theName);
		mDefaultDiagram = theDefaultDiagram;
		System.out.println("Find Children " + theSearchStart);
		// theCurrentPath, theLevel, theGDCFlag, null
		mChildren = getDisplayNodes(theSearchStart, 0, theLevelLabels, theTooltips);
		for (DisplayNode dn : mChildren)
		{
			dn.init(theTooltips);
		}
	}
	
	public void getDisplayNodesWithInit(Path theCurrentPath, String [] theLevelLabels, Path theDatasetDir, BEVtooltips theTooltips) throws IOException, Exception
	{
		int currentCount = theCurrentPath.getNameCount()+1;
		int datasetCount = theDatasetDir.getNameCount();
		//System.out.println("currentCount=" + currentCount);
		//System.out.println("datasetCount=" + datasetCount);
		TreeSet<DisplayNode> currentChildren = mChildren;
		for(int index=currentCount;index<=datasetCount; index++)
		{
			Path levelPath = Paths.get("/", theDatasetDir.subpath(0, index).toString());
			int levelNumber = index - currentCount;
			System.out.println("levelPath[" + levelNumber + "]=" + levelPath);
			if (levelPath.toFile().isDirectory())
			{
				if (new File(levelPath.toString(), mArchiveMarker).exists())
				{
					for(int defIndex=0;defIndex<mDefaultDiagram.length;defIndex++)
					{
						if ("*".equals(mDefaultDiagram[defIndex]))
						{
							mDefaultDiagram[defIndex] = theDatasetDir.getName(currentCount+defIndex).toString();
						}
					}
					Path zipFile = new File(levelPath.toString(), mArchiveName).toPath();
					ZipDisplayNode idc = new ZipDisplayNode(levelPath, levelNumber, theLevelLabels, zipFile, theTooltips, mArchiveMarker, mArchiveName);
					currentChildren.add(idc);
					idc.init(theTooltips);
				}
				else
				{
					DisplayNode newChild = new DisplayNode(levelPath, levelNumber, theLevelLabels, theTooltips, mArchiveMarker, mArchiveName);
					newChild.mChildren = new TreeSet<>();
					currentChildren.add(newChild);
					currentChildren = newChild.mChildren;
				}
			}
		}
	}
	
	public TreeSet<DisplayNode> getDisplayNodes(Path theCurrentPath, int theLevel, String [] theLevelLabels, BEVtooltips theTooltips) throws IOException, Exception
	{
		TreeSet<DisplayNode> subNodes = new TreeSet<>();
		System.out.println("DirDisplayName::getDisplayNodes theCurrentPath=" + theCurrentPath);
		Files.walkFileTree(theCurrentPath, EnumSet.noneOf(FileVisitOption.class), 1, new SimpleFileVisitor<Path>()
		{
			@Override
			public FileVisitResult visitFile(Path path, BasicFileAttributes attrs) throws IOException
			{
				if (attrs.isDirectory())
				{
					if (new File(path.toString(), mArchiveMarker).exists())
					{
						Path zipFile = new File(path.toString(), mArchiveName).toPath();
						ZipDisplayNode idc = new ZipDisplayNode(path, theLevel, theLevelLabels, zipFile, theTooltips, mArchiveMarker, mArchiveName);
						subNodes.add(idc);
					}
					else 
					{
						// Path theLocation, int theLevel, String theLabel, String theName, Path theZipFile
						subNodes.add(new DisplayNode(path, theLevel, theLevelLabels, theTooltips, mArchiveMarker, mArchiveName));
					}
				}
				return FileVisitResult.CONTINUE;
			}
		});
		return subNodes;
	}

	@Override
	public int compareTo(DisplayRun t)
	{
		return mName.compareTo(t.mName);
	}
	
	//
	
	public void removeZipPaths()
	{
		removeZipPaths(mChildren);
	}
	
	private void removeZipPaths(TreeSet<DisplayNode> theNodes)
	{
		for (DisplayNode dn : theNodes)
		{
			if (!"".equals(dn.mZipFile))
			{
				dn.mZipFile = new File(dn.mZipFile).getName();
			}
			if (null!=dn.mChildren)
			{
				removeZipPaths(dn.mChildren);
			}
		}
	}
	
	public void updateZipPaths(String theOrig, String theNew)
	{
		updateZipPaths(theOrig, theNew, mChildren);
	}
	
	private void updateZipPaths(String theOrig, String theNew, TreeSet<DisplayNode> theNodes)
	{
		for (DisplayNode dn : theNodes)
		{
			if (!"".equals(dn.mZipFile))
			{
				dn.mZipFile = dn.mZipFile.replace(theOrig, theNew);
			}
			if (null!=dn.mChildren)
			{
				updateZipPaths(theOrig, theNew, dn.mChildren);
			}
		}
	}
	
	public String getZipLocation(String thePath1, String thePath2, String thePath3, String thePath4)
	{
		String zipLocation = null;
		for (DisplayNode dn1 : this.mChildren)
		{
			if (null==zipLocation)
			{
				if (dn1.mName.equals(thePath1))
				{
					for (DisplayNode dn2 : dn1.mChildren)
					{
						if (dn2.mName.equals(thePath2))
						{
							for (DisplayNode dn3 : dn2.mChildren)
							{
								if (dn3.mName.equals(thePath3))
								{
									for (DisplayNode dn4 : dn3.mChildren)
									{
										if (dn4.mName.equals(thePath4))
										{
											zipLocation = dn4.mZipFile;
										}
									}
								}
							}
						}
					}
				}
			}
		}
		return zipLocation;
	}
	
	public TreeMap<String, String> dataRelations(String theIndexFile)
	{
		TreeMap<String, String> dataRelations = new TreeMap<>();
		for (DisplayNode dn : this.mChildren)
		{
			dn.dataRelations(dataRelations, theIndexFile);
		}
		return dataRelations;
	}
}
