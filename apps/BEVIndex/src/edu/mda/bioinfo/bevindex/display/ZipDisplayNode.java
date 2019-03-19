/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.mda.bioinfo.bevindex.display;

import edu.mda.bioinfo.bevindex.BEVtooltips;
import java.io.File;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.TreeSet;

/**
 *
 * @author linux
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
							InternalDisplayNode idn = new InternalDisplayNode(subPath.toFile().getAbsolutePath(), mLevel+1, mTransientLevelLabels, Paths.get(mZipFile), theTooltips,
									trimToSubDirs(subPath, false), "Standardized Data", new TreeSet<String>(), mArchiveMarker, mArchiveName);
							mChildren.add(idn);
						}
					}					
				}
			});
		};
		for (DisplayNode dn : mChildren)
		{
			dn.init(theTooltips);
		}
	}

}
