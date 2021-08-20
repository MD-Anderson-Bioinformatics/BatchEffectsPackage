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

package edu.mda.bcb.bevindex;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import edu.mda.bcb.bevindex.algorithms.BoxPlot;
import edu.mda.bcb.bevindex.algorithms.CDP;
import edu.mda.bcb.bevindex.algorithms.DSC;
import edu.mda.bcb.bevindex.algorithms.HierarchicalClustering;
import edu.mda.bcb.bevindex.algorithms.Discrete;
import edu.mda.bcb.bevindex.algorithms.MutBatch;
import edu.mda.bcb.bevindex.algorithms.Ngchm;
import edu.mda.bcb.bevindex.algorithms.PCAmtm;
import edu.mda.bcb.bevindex.algorithms.PCAotm;
import edu.mda.bcb.bevindex.algorithms.SupervisedClustering;
import edu.mda.bcb.bevindex.algorithms.TriNova;
import edu.mda.bcb.bevindex.algorithms.Umap;
import edu.mda.bcb.bevindex.gson.MBatchIndex;
import edu.mda.bcb.bevindex.gson.MBatchIndexData;
import edu.mda.bcb.bevindex.gson.OriginalData;
import edu.mda.bcb.bevindex.utils.ZipList;
import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;
import java.util.TreeSet;
import org.apache.commons.io.filefilter.RegexFileFilter;

/**
 *
 * @author Tod-Casasent
 */
public class BEVIndex
{
	static public String M_VERSION = "BEVIndex BEA_VERSION_TIMESTAMP";
	static public String M_AUTOCORRECT_NOTICE = "This dataset has been corrected using an automated system without human input. The correction does not imply the presence or absence of batch effects in the original data. The user is solely responsible for assessing batch effects (e.g. by using our assessment tools) and deciding whether or not to use the corrected data, which may or may not have mitigated some useful biological information along with any technical artifacts.";
	
	/** 
	 * 
	 * @param args 
	 */
	public static void main(String[] args)
	{
		System.out.println(M_VERSION);
		try
		{
			String mbatchID = args[0];
			String resultDir = args[1];
			String dataDir = args[2];
			String zipDir = args[3];

			runBEVIndex(mbatchID, resultDir, dataDir, zipDir);
		}
		catch (Exception exp)
		{
			exp.printStackTrace(System.err);
		}
	}
	
	/** 
	 * 
	 * @param theMbatchID
	 * @param theResultDir
	 * @param theDataDir
	 * @param theZipDir
	 * @throws FileNotFoundException
	 * @throws IOException 
	 */
	static public void runBEVIndex(String theMbatchID, String theResultDir, String theDataDir, String theZipDir) throws FileNotFoundException, IOException
	{
		System.out.println("runBEVIndex starting");
		System.out.println(M_VERSION);
		// read original data
		GsonBuilder builder = new GsonBuilder();
		builder.setPrettyPrinting();
		Gson gson = builder.create();
		OriginalData od = null;
		File odJson = new File(theResultDir, "original_data.json");
		if (odJson.exists())
		{
			System.out.println("runBEVIndex odJson found:"+ odJson.getAbsolutePath());
			try(FileReader fr = new FileReader(odJson))
			{
				od = gson.fromJson(fr, OriginalData.class);
			}
		}
		// read source id
		File sourceFile = new File(theResultDir, "source_id.txt");
		String sourceID =  "none";
		if (sourceFile.exists())
		{
			System.out.println("runBEVIndex sourceFile found:"+ sourceFile.getAbsolutePath());
			sourceID = new String(Files.readAllBytes(sourceFile.toPath()));
		}
		// String versionStamp = "2020_03_13_1313";
		File versionStampFile = new File(theResultDir, "version_stamp.txt");
		String versionStamp = "none";
		if (versionStampFile.exists())
		{
			System.out.println("runBEVIndex versionStampFile found:"+ versionStampFile.getAbsolutePath());
			versionStamp = new String(Files.readAllBytes(versionStampFile.toPath()));
		}
		// String versionType = "All-original";
		File versionTypeFile = new File(theResultDir, "version_type.txt");
		String versionType = "none";
		if (versionTypeFile.exists())
		{
			System.out.println("runBEVIndex versionTypeFile found:"+ versionTypeFile.getAbsolutePath());
			versionType = new String(Files.readAllBytes(versionTypeFile.toPath()));
			System.out.println("runBEVIndex versionType:"+ versionType);
		}
		//
		String notice = null;
		if ((versionType.equals("EB_withPara"))||(versionType.equals("EB_withNonpara"))||
			(versionType.equals("MP_overall"))||(versionType.equals("MP_batch"))||
			(versionType.equals("ANOVA_adj"))||(versionType.equals("ANOVA_unadj"))||
			(versionType.equals("RBN_Replicates"))||(versionType.equals("RBN_Pseudoreps"))||
			(versionType.equals("EBN_Plus")) )
		{
			notice = M_AUTOCORRECT_NOTICE;
		}
		// process dir
		System.out.println("runBEVIndex call buildIndexObject");
		MBatchIndex mi = buildIndexObject(sourceID, theMbatchID, versionStamp, versionType, od, new File(theResultDir), notice);
		// have index object, write it
		System.out.println("runBEVIndex write index file");
		File indexFile = new File(theResultDir, "index.json");
		Files.write(indexFile.toPath(), gson.toJson(mi).getBytes(), StandardOpenOption.WRITE, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.CREATE);
		// ZIP data dir
		System.out.println("runBEVIndex theDataDir:"+ theDataDir);
		if (null!=theDataDir)
		{
			System.out.println("runBEVIndex not null theDataDir:"+ theDataDir);
			File dataDir = new File(theDataDir);
			if (dataDir.exists())
			{
				System.out.println("runBEVIndex dir exists theDataDir:"+ theDataDir);
				ZipList zl = new ZipList(dataDir);
				zl.lsInclude();
				zl.zipList(new File(theZipDir, theMbatchID+"-data.zip").toPath());
			}
		}
		// ZIP analysis dir
		System.out.println("runBEVIndex theResultDir:"+ theResultDir);
		File resultDir = new File(theResultDir);
		if (resultDir.exists())
		{
			System.out.println("runBEVIndex dir exists theResultDir:"+ theResultDir);
			ZipList zl = new ZipList(resultDir);
			zl.lsInclude();
			zl.zipList(new File(theZipDir, theMbatchID+"-results.zip").toPath());
		}
	}
	
	
	static public String checkForPattern(File theDir, String thePattern)
	{
		System.out.println("checkForPattern theDir=" + theDir.getAbsolutePath());
		System.out.println("checkForPattern thePattern=" + thePattern);
		String file = "";
		FileFilter filter = new RegexFileFilter(thePattern);
		File[] files = theDir.listFiles(filter);
		if ((null!=files)&&(files.length>0))
		{
			file = files[0].getName();
		}
		return file;
	}
	
	static public String checkForFile(File theFile, String thePath)
	{
		String file = "";
		if (theFile.exists())
		{
			file = thePath;
		}
		return file;
	}
	
	static public MBatchIndex buildIndexObject(String theSourceId, String theMbatchId,
			String theVersionStamp, String theVersionType, 
			OriginalData theOD, File theResultsDir, String theNotice)
	{
		MBatchIndex mi = new MBatchIndex();
		if (null==theNotice)
		{
			mi.notice = "";
		}
		else
		{
			mi.notice = theNotice;
		}
		if (null!=theOD)
		{
			mi.algorithm = theOD.algorithm;
			mi.category = theOD.category;
			mi.data = theOD.data;
			mi.details = theOD.details;
			mi.platform = theOD.platform;
			mi.project = theOD.project;
			mi.source = theOD.source;
			mi.subProject = theOD.subProject;
			mi.variant = theOD.variant;
			mi.version = theOD.version;
		}
		else
		{
			mi.algorithm = "unknown - user data";
			mi.category = "unknown - user data";
			mi.data = "unknown - user data";
			mi.details = "unknown - user data";
			mi.platform = "unknown - user data";
			mi.project = "unknown - user data";
			mi.source = "unknown - user data";
			mi.subProject = "unknown - user data";
			mi.variant = "unknown - user data";
			mi.version = "unknown - user data";
		}
		mi.SourceID = theSourceId;
		// mbatch info
		mi.mbatch = new MBatchIndexData();
		mi.mbatch.MBatchID = theMbatchId;
		mi.mbatch.mbatch_run = theVersionStamp;
		mi.mbatch.dataset_type = theVersionType;
		//"ANY_Corrections-ANAdjusted.tsv"
		mi.mbatch.matrix = checkForPattern(theResultsDir, "ANY_Corrections-.+\\.tsv");
		mi.mbatch.batchdata = checkForFile(new File(theResultsDir, "BatchData.tsv"), "BatchData.tsv");
		mi.mbatch.dropdown_label = "Algorithm";
		mi.mbatch.dropdown_entries = new TreeSet<>();
		// check for each subdir algorithm
		// BoxPlot Boxplot diagram_type=box
		BoxPlot bp = new BoxPlot(new File(theResultsDir, "BoxPlot"), mi.mbatch);
		bp.checkAndAddAlgorithm();
		// CDP diagram_type=cdp
		CDP cdp = new CDP(new File(theResultsDir, "CDP"), mi.mbatch);
		cdp.checkAndAddAlgorithm();
		// UMAP diagram_type=umap
		Umap umap = new Umap(new File(theResultsDir, "UMAP"), mi.mbatch);
		umap.checkAndAddAlgorithm();
		// TriNova diagram_type=tri-man
		TriNova tm = new TriNova(new File(theResultsDir, "TRINOVA"), mi.mbatch);
		tm.checkAndAddAlgorithm();
		// HierarchicalClustering diagram_type=hc
		HierarchicalClustering hc = new HierarchicalClustering(new File(theResultsDir, "HierarchicalClustering"), mi.mbatch);
		hc.checkAndAddAlgorithm();
		// NGCHM diagram_type=ngchm
		Ngchm chm = new Ngchm(new File(theResultsDir, "NGCHM"), mi.mbatch);
		chm.checkAndAddAlgorithm();
		// DSC diagram_type=dsc
		DSC dsc = new DSC(theResultsDir, mi.mbatch);
		dsc.checkAndAddAlgorithm();
		// PCA mtm diagram_type=pca
		PCAmtm pcamtm = new PCAmtm(new File(theResultsDir, "PCA"), mi.mbatch);
		pcamtm.checkAndAddAlgorithm();
		// PCA otm diagram_type=pca
		PCAotm pcaotm = new PCAotm(new File(theResultsDir, "PCA"), mi.mbatch);
		pcaotm.checkAndAddAlgorithm();
		// SupervisedClustering diagram_type=sc
		SupervisedClustering sc = new SupervisedClustering(new File(theResultsDir, "SupervisedClustering"), mi.mbatch);
		sc.checkAndAddAlgorithm();
		// Discrete diagram_type=discrete
		Discrete discrete = new Discrete(new File(theResultsDir, "Discrete"), mi.mbatch);
		discrete.checkAndAddAlgorithm();
		// MutBatch diagram_type=mutbatch
		MutBatch mutbat = new MutBatch(new File(theResultsDir, "MutBatch"), mi.mbatch);
		mutbat.checkAndAddAlgorithm();
		//
		return mi;
	}
}
