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
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.TreeSet;

/**
 *
 * @author Tod-Casasent
 */
public class ZipDisplayNode extends DisplayNode
{

	public ZipDisplayNode(Path theLocation, int theLevel, String[] theLevelLabels, Path theZipFile, BEVtooltips theTooltips, String theArchiveMarker, String theArchiveName)
	{
		super(theLocation, theLevel, theLevelLabels, theTooltips, theArchiveMarker, theArchiveName);
		mDiagramFlag = false;	// set in descendants
		mZipFile = theZipFile.toFile().getAbsolutePath();
		mInternalLocation = "";	// set in descendants
		mAlgorithm = "";		// set in descendants
		mOtherFiles = null;		// set in descendants
		//mNodeType = "ZipDisplayNode";
	}

	protected String trimToSubDirs(Path theSubPath, boolean addTrailing)
	{
		// use get parent, since zips include base directory of parent
		String trimmed = theSubPath.toFile().getAbsolutePath().replaceFirst(new File(mZipFile).getParentFile().getParentFile().getAbsolutePath(), "");
		trimmed = trimmed.replaceFirst("/", "");
		if (addTrailing)
		{
			trimmed = trimmed + "/";
		}
		return trimmed;
	}

	@Override
	public void init(BEVtooltips theTooltips) throws Exception
	{
		System.out.println("Switch to internal at:" + mTransientLocation);
		System.out.println("mZipFile=" + mZipFile);
		// use DirectoryStream, since it is faster when we have hundreds of files
		if (null==mChildren)
		{
			mChildren = new TreeSet<>();
		}
		try (DirectoryStream<Path> stream = Files.newDirectoryStream(mTransientLocation))
		{
			stream.forEach((Path subPath) ->
			{
				if (subPath.toFile().isDirectory())
				{
					// String theDir, int theLevel, String [] theLevelLabels, Path theZipFile, BEVtooltips theTooltips,
					// String theArchiveMarker, String theArchiveName
					InternalDisplayNode idn = new InternalDisplayNode(subPath.toFile().getAbsolutePath(), mLevel+1, mTransientLevelLabels, Paths.get(mZipFile), theTooltips, mArchiveMarker, mArchiveName);
					mChildren.add(idn);
				}
				else
				{
					if (subPath.toFile().isFile())
					{
						if (subPath.getFileName().toString().equalsIgnoreCase("RunInfo.PNG"))
						{
							// For normal data, look for RunInfo.PNG
							InternalDisplayNode idn = new InternalDisplayNode(subPath.toFile().getAbsolutePath(), mLevel+1, mTransientLevelLabels, Paths.get(mZipFile), theTooltips,
									trimToSubDirs(subPath, false), "Standardized Data", new TreeSet<String>(), mArchiveMarker, mArchiveName);
							mChildren.add(idn);
						}
						else if (subPath.getFileName().toString().equalsIgnoreCase("DSCOverview.tsv"))
						{
							// handle DSC-Overview, which does not have RunInfo.PNG
							InternalDisplayNode idn = new InternalDisplayNode(subPath.toFile().getAbsolutePath(), mLevel+1, mTransientLevelLabels, Paths.get(mZipFile), theTooltips,
									trimToSubDirs(subPath, false), "DSC", new TreeSet<String>(), mArchiveMarker, mArchiveName);
							mChildren.add(idn);
						}
					}					
				}
			});
		}
		for (DisplayNode dn : mChildren)
		{
			dn.init(theTooltips);
		}
	}

}
