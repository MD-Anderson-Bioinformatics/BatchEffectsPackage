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
public class InternalDisplayNode extends ZipDisplayNode
{
	transient private String mCurrentDir;

	// internal node only, not a diagram
	public InternalDisplayNode(String theDir, int theLevel, String [] theLevelLabels, Path theZipFile, BEVtooltips theTooltips,
			String theArchiveMarker, String theArchiveName)
	{
		super(Paths.get(theDir), theLevel, theLevelLabels, theZipFile, theTooltips, theArchiveMarker, theArchiveName);
		mCurrentDir = theDir;
		mChildren = new TreeSet<>();
	}

	// diagram location
		//, TreeSet<String> theDownloadableFiles
	public InternalDisplayNode(String theFilePath, int theLevel, String [] theLevelLabels, Path theZipFile, BEVtooltips theTooltips, 
			String theInternalLocation, String theAlgorithm, TreeSet<String> theOtherFiles, String theArchiveMarker, String theArchiveName)
	{
		super(Paths.get(theFilePath), theLevel, theLevelLabels, theZipFile, theTooltips, theArchiveMarker, theArchiveName);
		mLabel = "Diagram";
		mDiagramFlag = true;
		mInternalLocation = theInternalLocation;
		mAlgorithm = theAlgorithm;
		mOtherFiles = theOtherFiles;
		mChildren = new TreeSet<>();
	}
	
	@Override
	public void init(BEVtooltips theTooltips) throws Exception
	{
		if (mTransientLocation.toFile().isDirectory())
		{
			////////////////////////////////////////////////////////////////////////
			// check for diagrams in current dir (stop for all but DSC)
			Boolean stopSearch = null;
			mChildren = new TreeSet<>();
			// check for diagrams, if so process diagrams
			// use DirectoryStream, since it is faster when we have hundreds of files
			try (DirectoryStream<Path> stream = Files.newDirectoryStream(mTransientLocation))
			{
				for(Path subPath : stream)
				{
					if (subPath.toFile().isFile())
					{
						//System.out.println("IDN:init=" + subPath);
						String filename = subPath.getFileName().toString();
						File[] directories = subPath.toFile().getParentFile().listFiles(File::isDirectory);
						boolean hasSubDirs = directories.length>0;
						String trimmedFilename = trimToSubDirs(subPath, false);
						String algorithm = null;
						TreeSet<String> otherFiles = null;
						if (filename.matches("(BoxPlot_)(.*)(_BoxData-)(.*)(\\.tsv)"))
						{
							// add Boxplot DIAGRAM DisplayNode
							algorithm = "Boxplot";
							otherFiles = otherFilesBoxplot(trimmedFilename);
							// stop search here, "if" is so that DSCOverview.tsv "false" is not overwritten
							if (null==stopSearch)
							{
								stopSearch = true;
							}
						}
						else if (filename.matches("(FullMutCounts_)(.*)(\\.PNG)"))
						{
							// add MutBatch DisplayNode
							algorithm = "MutBatch";
							otherFiles = otherFilesMutBatch(trimmedFilename);
							// stop search here, "if" is so that DSCOverview.tsv "false" is not overwritten
							if (null==stopSearch)
							{
								stopSearch = true;
							}
						}
						else if (filename.equals("PCAValues.tsv"))
						{
							// add PCA DisplayNode
							algorithm = "PCA";
							otherFiles = otherFilesPCA(trimmedFilename);
							// stop search here, "if" is so that DSCOverview.tsv "false" is not overwritten
							if (null==stopSearch)
							{
								stopSearch = true;
							}
						}
						else if (filename.equals("HCData.tsv"))
						{
							// add HierarchicalClustering DisplayNode
							algorithm = "HierarchicalClustering";
							otherFiles = otherFilesHierarchicalClustering(trimmedFilename);
							// stop search here, "if" is so that DSCOverview.tsv "false" is not overwritten
							if (null==stopSearch)
							{
								stopSearch = true;
							}
						}
						else if (filename.matches("(SupervisedClust_Diagram-)(.*)(\\.png)"))
						{
							// add SupervisedClustering DisplayNode
							algorithm = "SupervisedClustering";
							otherFiles = otherFilesSupervisedClustering(trimmedFilename);
							// stop search here, "if" is so that DSCOverview.tsv "false" is not overwritten
							if (null==stopSearch)
							{
								stopSearch = true;
							}
						}
						else if ((filename.matches("(CDP_Plot_)(.*)(_Diagram\\.png)")) || (filename.matches("(CDP_Plot_)(.*)(_Diagram\\.PNG)")))
						{
							// add CDP DisplayNode
							algorithm = "CDP";
							otherFiles = otherFilesCDP(trimmedFilename);
							// stop search here, "if" is so that DSCOverview.tsv "false" is not overwritten
							if (null==stopSearch)
							{
								stopSearch = true;
							}
						}
						else if (filename.matches("(.*)(_ngchm\\.ngchm)"))
						{
							// add NGCHM DisplayNode
							algorithm = "NGCHM";
							otherFiles = otherFilesNGCHM(trimmedFilename);
							// stop search here, "if" is so that DSCOverview.tsv "false" is not overwritten
							if (null==stopSearch)
							{
								stopSearch = true;
							}
						}
						else if (filename.equalsIgnoreCase("RunInfo.png"))
						{
							// add CDP DisplayNode
							algorithm = "Standardized Data";
							otherFiles = otherFilesStandardizedData(trimmedFilename);
							// stop search here, "if" is so that DSCOverview.tsv "false" is not overwritten
							if (null==stopSearch)
							{
								stopSearch = true;
							}
						}
						else if (filename.equals("DSCOverview.tsv"))
						{
							// add DSC Overview Node
							algorithm = "DSC";
							otherFiles = otherFilesDSC(trimmedFilename);
							// DO NOT stop search here if there are subdirs
							if (hasSubDirs)
							{
								stopSearch = false;
							}
							else
							{
								stopSearch = true;
							}
						}
						if (null!=algorithm)
						{
							// 	public InternalDisplayNode(String theFilePath, int theLevel, String [] theLevelLabels, Path theZipFile, BEVtooltips theTooltips, 
							//  String theInternalLocation, String theAlgorithm, TreeSet<String> theOtherFiles, String theArchiveMarker, String theArchiveName)
							InternalDisplayNode idn = new InternalDisplayNode(trimmedFilename, mLevel+1, mTransientLevelLabels, 
									Paths.get(mZipFile), theTooltips, trimToSubDirs(subPath, false),
									algorithm, otherFiles, mArchiveMarker, mArchiveName);
							mChildren.add(idn);
						}
					}
				}
			}

			// if no diagrams found, and not the DSCOverview setup, then continue into subdirs
			if ((null==stopSearch)||(false==stopSearch))
			{
				try (DirectoryStream<Path> stream = Files.newDirectoryStream(mTransientLocation))
				{
					for(Path subPath : stream)
					{
						if (subPath.toFile().isDirectory())
						{
							InternalDisplayNode idn = new InternalDisplayNode(subPath.toFile().getAbsolutePath(), mLevel+1, mTransientLevelLabels, 
									Paths.get(mZipFile), theTooltips, mArchiveMarker, mArchiveName);
							mChildren.add(idn);
						}
					}
				}
				for (DisplayNode dn : mChildren)
				{
					dn.init(theTooltips);
				}
			}
		}
	}
	
	protected TreeSet<String> otherFilesStandardizedData(String theFile)
	{
		TreeSet<String> others = new TreeSet<>();
		others.add(new File(new File(theFile).getParentFile(), "RunInfo.tsv").getAbsolutePath());
		others.add(new File(new File(theFile).getParentFile(), "release.tsv").getAbsolutePath());
		return others;
	}

	protected TreeSet<String> otherFilesBoxplot(String theFile)
	{
		TreeSet<String> others = new TreeSet<>();
		// using TSV file, so replace "_BoxData-" with "_Annotations-"
		others.add(theFile.replace("_BoxData-", "_Annotations-"));
		// using TSV file, so replace "_BoxData-" with "_Histogram-"
		others.add(theFile.replace("_BoxData-", "_Histogram-"));
		return others;
	}

	protected TreeSet<String> otherFilesMutBatch(String theFile)
	{
		TreeSet<String> others = new TreeSet<>();
		others.add(theFile.replace("_Diagram.PNG", ".tsv"));
		others.add(Paths.get(Paths.get(theFile).getParent().toString(), "callReference.tsv").toString());
		return others;
	}

	protected TreeSet<String> otherFilesPCA(String theFile)
	{
		TreeSet<String> others = new TreeSet<>();
		// use PCAAnnotations.tsv as other file
		others.add(Paths.get(Paths.get(theFile).getParent().toString(), "PCAAnnotations.tsv").toString());
		return others;
	}

	protected TreeSet<String> otherFilesHierarchicalClustering(String theFile)
	{
		TreeSet<String> others = new TreeSet<>();
		// use HCOrder.tsv as other file
		others.add(Paths.get(Paths.get(theFile).getParent().toString(), "HCOrder.tsv").toString());
		return others;
	}

	protected TreeSet<String> otherFilesSupervisedClustering(String theFile)
	{
		TreeSet<String> others = new TreeSet<>();
		// using PNG file, so replace "_Diagram-" with "_Legend-"
		others.add(theFile.replace("_Diagram-", "_Legend-"));
		return others;
	}

	protected TreeSet<String> otherFilesCDP(String theFile)
	{
		TreeSet<String> others = new TreeSet<>();
		// no files to add
		return others;
	}

	protected TreeSet<String> otherFilesNGCHM(String theFile)
	{
		TreeSet<String> others = new TreeSet<>();
		// no files to add
		return others;
	}

	protected TreeSet<String> otherFilesDSC(String theFile)
	{
		TreeSet<String> others = new TreeSet<>();
		// no files to add
		return others;
	}
}
