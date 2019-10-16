/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.mda.bioinfo.bevmanifest;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeSet;

/**
 *
 * @author linux
 */
public class BEVManifest
{

	static protected String mVersion = "BEVManifest 2019-10-08-1313";

	// MBATCH_SUCCESS.txt
	// Boxplot/ AllSample-Data  AllSample-RLE  Group-Mean -> *BoxData*.tsv
	// CDP CDP_*_Diagram.PNG
	// HierarchicalClustering HCData.tsv
	// NGCHM *.ngchm >=1
	// PCA/ BatchId  PlateId  ShipDate  TSS -> PCAValues.tsv
	// SupervisedClustering/ Batches -> SupervisedClust_Diagram-*.png
	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args)
	{
		System.out.println("BEVManifest");
		System.out.println(mVersion);
		// ARGUMENTS
		String searchDir = args[0];
		String fileName = args[1];
		String manifestRun = args[2];
		String manifestName = args[3];
		try
		{
			System.out.println("search:" + searchDir);
			System.out.println("for:" + fileName);
			FileFind ff = new FileFind();
			ff.find(Paths.get(searchDir), fileName);
			TreeSet<Path> toZipSet = ff.mMatches;
			List<Exception> errors = Collections.synchronizedList(new ArrayList<>());
			try (BufferedWriter runWriter = java.nio.file.Files.newBufferedWriter(Paths.get(manifestRun), Charset.availableCharsets().get("UTF-8")))
			{
				runWriter.write("Dataset\tBoxplot\tCDP\tHierClust\tNGCHM\tPCA\tSuperClust");
				runWriter.newLine();
				//System.out.println("PARALLEL-THREADS = 5");
				//System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", "5");
				//toZipSet.parallelStream()
				toZipSet.stream()
						.forEach((Path marker) ->
						{
							String dataDir = marker.getParent().toString();
							System.out.println("dataDir:" + dataDir);
							String dataset = dataDir.replace(searchDir, "");
							if (dataset.startsWith(File.separator))
							{
								dataset = dataset.substring(1);
							}
							String cntBoxplot = "NA";
							String cntCDP = "NA";
							String cntHierClust = "NA";
							String cntNGCHM = "NA";
							String cntPCA = "NA";
							String cntSuperClust = "NA";
							try
							{
								cntBoxplot = countFiles(dataDir, "Boxplot", ".*BoxData.*.tsv");
							}
							catch (Exception exp)
							{
								errors.add(new Exception("Error processing " + "Boxplot", exp));
							}
							try
							{
								cntCDP = countFiles(dataDir, "CDP", "CDP_.*_Diagram.PNG");
							}
							catch (Exception exp)
							{
								errors.add(new Exception("Error processing " + "CDP", exp));
							}
							try
							{
								cntHierClust = countFiles(dataDir, "HierarchicalClustering", "HCData.tsv");
							}
							catch (Exception exp)
							{
								errors.add(new Exception("Error processing " + "HierarchicalClustering", exp));
							}
							try
							{
								cntNGCHM = countFiles(dataDir, "NGCHM", ".*.ngchm");
							}
							catch (Exception exp)
							{
								errors.add(new Exception("Error processing " + "NGCHM", exp));
							}
							try
							{
								cntPCA = countFiles(dataDir, "PCA", "PCAValues.tsv");
							}
							catch (Exception exp)
							{
								errors.add(new Exception("Error processing " + "PCA", exp));
							}
							try
							{
								cntSuperClust = countFiles(dataDir, "SupervisedClustering", "SupervisedClust_Diagram-.*.png");
							}
							catch (Exception exp)
							{
								errors.add(new Exception("Error processing " + "SupervisedClustering", exp));
							}
							try
							{
								runWriter.write(dataset + "\t" + cntBoxplot + "\t" + cntCDP + "\t" + cntHierClust + "\t" + cntNGCHM + "\t" + cntPCA + "\t" + cntSuperClust);
								runWriter.newLine();
							}
							catch (Exception exp)
							{
								errors.add(new Exception("Error writting to " + manifestRun, exp));
							}
							File datasetManifest = new File(dataDir, manifestName);
							try (BufferedWriter datasetWriter = java.nio.file.Files.newBufferedWriter(datasetManifest.toPath(), Charset.availableCharsets().get("UTF-8")))
							{
								datasetWriter.write("Dataset\tBoxplot\tCDP\tHierClust\tNGCHM\tPCA\tSuperClust");
								datasetWriter.newLine();
								datasetWriter.write(dataset + "\t" + cntBoxplot + "\t" + cntCDP + "\t" + cntHierClust + "\t" + cntNGCHM + "\t" + cntPCA + "\t" + cntSuperClust);
								datasetWriter.newLine();
							}
							catch (Exception exp)
							{
								errors.add(new Exception("Error writting to " + manifestRun, exp));
							}
						});
			}
			for (Exception exp : errors)
			{
				exp.printStackTrace(System.err);
				exp.getCause().printStackTrace(System.err);
			}
		}
		catch (Exception exp)
		{
			exp.printStackTrace(System.err);
		}
	}

	static public String countFiles(String theDataDir, String theSubDir, String thePattern) throws IOException
	{
		String result = "NA";
		File subDir = new File(theDataDir, theSubDir);
		if (subDir.exists())
		{
			long count = Files.find(subDir.toPath(), Integer.MAX_VALUE,
					(thePath, theBasicFileAttributes) -> thePath.toFile().getName().matches(thePattern)).count();
			result = Long.toString(count);
		}
		return result;
	}
}
