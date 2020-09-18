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

package edu.mda.bioinfo.bevindex;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import edu.mda.bioinfo.bevindex.algorithms.BoxPlot;
import edu.mda.bioinfo.bevindex.algorithms.CDP;
import edu.mda.bioinfo.bevindex.algorithms.DSC;
import edu.mda.bioinfo.bevindex.algorithms.HierarchicalClustering;
import edu.mda.bioinfo.bevindex.algorithms.Discrete;
import edu.mda.bioinfo.bevindex.algorithms.MutBatch;
import edu.mda.bioinfo.bevindex.algorithms.Ngchm;
import edu.mda.bioinfo.bevindex.algorithms.PCAmtm;
import edu.mda.bioinfo.bevindex.algorithms.PCAotm;
import edu.mda.bioinfo.bevindex.algorithms.SupervisedClustering;
import edu.mda.bioinfo.bevindex.gson.MBatchIndex;
import edu.mda.bioinfo.bevindex.gson.MBatchIndexData;
import edu.mda.bioinfo.bevindex.gson.OriginalData;
import edu.mda.bioinfo.bevindex.utils.ZipChildren;
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
	static public String M_VERSION = "BEVIndex 2020-09-11-1000";
	static public String M_AUTOCORRECT_NOTICE = "This dataset has been corrected using an automated system without human input. The correction does not imply the presence or absence of batch effects in the original data. The user is solely responsible for assessing batch effects (e.g. by using our assessment tools) and deciding whether or not to use the corrected data, which may or may not have mitigated some useful biological information along with any technical artifacts.";
	
	////////////////////////////////////////////////////////////////////////////
	//// main
	////////////////////////////////////////////////////////////////////////////
	public static void main(String[] args)
	{
		System.out.println(M_VERSION);
		try
		{
			// "MBATCH_ID_0000";
			//String mbatchID = "1585774011149";
			String mbatchID = args[0];
			// ".../2020_03_12_1022/original_data.json";
			//String originalDataJsonFile = "//code/dvlp-data/1585774011149/original_data.json";
			String originalDataJsonFile = args[1];
			// ".../2020_03_12_1022/MBatch/";
			//String mbatchResultsDir = "/code/dvlp-data/1585774011149/MBatch/1585774011149/";
			String mbatchResultsDir = args[2];
			// ".../2020_03_12_1022/";
			//String zipDir = "/code/dvlp-data/1585774011149/";
			String zipDir = args[3];
			runBEVIndex(mbatchID, originalDataJsonFile, mbatchResultsDir, zipDir);
		}
		catch (Exception exp)
		{
			exp.printStackTrace(System.err);
		}
	}
	
	static public void runBEVIndex(String mbatchID, String originalDataJsonFile, String mbatchResultsDir, String zipDir) throws FileNotFoundException, IOException
	{
		System.out.println(M_VERSION);
		// read original data
		GsonBuilder builder = new GsonBuilder();
		builder.setPrettyPrinting();
		Gson gson = builder.create();
		OriginalData od = null;
		if (null!=originalDataJsonFile)
		{
			if (new File(originalDataJsonFile).exists())
			{
				try(FileReader fr = new FileReader(originalDataJsonFile))
				{
					od = gson.fromJson(fr, OriginalData.class);
				}
			}
		}
		// read source id
		File sourceFile = new File(zipDir, "source_id.txt");
		String sourceID =  "none";
		if (sourceFile.exists())
		{
			sourceID = new String(Files.readAllBytes(sourceFile.toPath()));
		}
		// String versionStamp = "2020_03_13_1313";
		File versionStampFile = new File(zipDir, "version_stamp.txt");
		String versionStamp = "none";
		if (versionStampFile.exists())
		{
			versionStamp = new String(Files.readAllBytes(versionStampFile.toPath()));
		}
		// String versionType = "All-original";
		File versionTypeFile = new File(zipDir, "version_type.txt");
		String versionType = "none";
		if (versionTypeFile.exists())
		{
			versionType = new String(Files.readAllBytes(versionTypeFile.toPath()));
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
		MBatchIndex mi = buildIndexObject(sourceID, mbatchID, versionStamp, versionType, od, new File(mbatchResultsDir), notice);
		// have index object, write it
		File indexFile = new File(mbatchResultsDir, "index.json");
		Files.write(indexFile.toPath(), gson.toJson(mi).getBytes(), StandardOpenOption.WRITE, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.CREATE);
		// zip directory contents
		ZipChildren zd = new ZipChildren();
		zd.zipChildren(new File(mbatchResultsDir).toPath(), new File(zipDir, mbatchID+".zip").toPath());
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
