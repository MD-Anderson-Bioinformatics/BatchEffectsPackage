/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.mda.bioinfo.bevmanifest;

import static edu.mda.bioinfo.bevmanifest.BEVManifest.mVersion;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.nio.charset.Charset;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Scanner;
import java.util.TreeSet;

/**
 *
 * @author linux
 */
public class SDBManifest
{

	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args)
	{
		System.out.println("SDBManifest");
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
				runWriter.write("Dataset\tMatrix-Samples\tMatrix-Features\tBatches-Samples\tBatches-Types\tAnnotations-Samples\tAnnotations-Values\tClinical-Samples\tClinical-Values");
				runWriter.newLine();
				//System.out.println("PARALLEL-THREADS = 5");
				//System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", "5");
				//toZipSet.parallelStream()
				toZipSet.stream()
						.forEach((Path marker) ->
						{
							String dataDir = marker.getParent().toString();
							System.out.println("dataDir:" + dataDir);
							File matrixFile = new File(dataDir, "matrix_data.tsv");
							File batchesFile = new File(dataDir, "batches.tsv");
							File clinicalFile = new File(dataDir, "clinical.tsv");
							File annotationsFile = new File(dataDir, "annotations.tsv");
							String dataset = dataDir.replace(searchDir, "");
							if (dataset.startsWith(File.separator))
							{
								dataset = dataset.substring(1);
							}
							String matrixSamples = "NA";
							String matrixFeatures = "NA";
							String batchesSamples = "NA";
							String batchesBatchTypes = "NA";
							String clinicalSamples = "NA";
							String clinicalValues = "NA";
							String annotationSamples = "NA";
							String annotationValues = "NA";
							try
							{
								String colRows = getMatrixString(matrixFile);
								String [] splitted = colRows.split(",");
								matrixSamples = splitted[0];
								matrixFeatures = splitted[1];
							}
							catch (Exception exp)
							{
								errors.add(new Exception("Error processing " + matrixFile, exp));
							}
							try
							{
								String colRows = getDFString(batchesFile);
								String [] splitted = colRows.split(",");
								batchesBatchTypes = splitted[0];
								batchesSamples = splitted[1];
								
							}
							catch (Exception exp)
							{
								errors.add(new Exception("Error processing " + batchesFile, exp));
							}
							try
							{
								String colRows = getDFString(clinicalFile);
								String [] splitted = colRows.split(",");
								clinicalValues = splitted[0];
								clinicalSamples = splitted[1];
							}
							catch (Exception exp)
							{
								errors.add(new Exception("Error processing " + clinicalFile, exp));
							}
							try
							{
								String colRows = getDFString(annotationsFile);
								String [] splitted = colRows.split(",");
								annotationValues = splitted[0];
								annotationSamples = splitted[1];
							}
							catch (Exception exp)
							{
								errors.add(new Exception("Error processing " + annotationsFile, exp));
							}
							try
							{
								runWriter.write(dataset + "\t" + matrixSamples + "\t" + matrixFeatures + "\t" + batchesSamples + "\t" + batchesBatchTypes + "\t" + 
										clinicalSamples + "\t" + clinicalValues + "\t" + annotationSamples + "\t" + annotationValues);
								runWriter.newLine();
								runWriter.flush();
							}
							catch (Exception exp)
							{
								errors.add(new Exception("Error writting to " + manifestRun, exp));
							}
							File datasetManifest = new File(dataDir, manifestName);
							try(BufferedWriter datasetWriter = java.nio.file.Files.newBufferedWriter(datasetManifest.toPath(), Charset.availableCharsets().get("UTF-8")))
							{
								datasetWriter.write("Dataset\tMatrix-Samples\tMatrix-Features\tBatches-Samples\tBatches-Types\tAnnotations-Samples\tAnnotations-Values\tClinical-Samples\tClinical-Values");
								datasetWriter.newLine();
								datasetWriter.write(dataset + "\t" + matrixSamples + "\t" + matrixFeatures + "\t" + batchesSamples + "\t" + batchesBatchTypes + "\t" + 
										clinicalSamples + "\t" + clinicalValues + "\t" + annotationSamples + "\t" + annotationValues);
								datasetWriter.newLine();
							}
							catch (Exception exp)
							{
								errors.add(new Exception("Error writing to " + datasetManifest, exp));
							}
						});
			}
			for (Exception exp : errors)
			{
				System.err.println(exp.getMessage());
				exp.printStackTrace(System.err);
				exp.getCause().printStackTrace(System.err);
			}
		}
		catch (Exception exp)
		{
			exp.printStackTrace(System.err);
		}
	}

	static public String getMatrixString(File theMatrixFile) throws FileNotFoundException
	{
		String result = "NA,NA";
		if (theMatrixFile.exists())
		{
			// count columns from first line minus one
			int columns = -1;
			// then count rows after first
			int rows = 0;
			try (Scanner scanner = new Scanner(theMatrixFile))
			{
				while (scanner.hasNext())
				{
					String line = scanner.nextLine();
					if (columns<0)
					{
						columns = line.split("\t").length - 1;
					}
					else
					{
						rows = rows + 1;
					}
				}
			}
			result = columns + "," + rows;
		}
		return result;
	}
	
	static public String getDFString(File theDFFile) throws FileNotFoundException
	{
		String result = "NA,NA";
		if (theDFFile.exists())
		{
			// count columns from first line minus one
			int columns = -1;
			// then count rows after first, that contain more than one column without "Unknown"
			int rows = 0;
			try (Scanner scanner = new Scanner(theDFFile))
			{
				while (scanner.hasNext())
				{
					String line = scanner.nextLine();
					if (columns<0)
					{
						columns = line.split("\t", -1).length - 1;
					}
					else
					{
						// check if more than 1 column has non "Unknown" entries
						int cnt = 0;
						for (String val : line.split("\t", -1))
						{
							if (!val.equalsIgnoreCase("Unknown"))
							{
								cnt = cnt + 1;
							}
						}
						if (cnt>1)
						{
							rows = rows + 1;
						}
					}
				}
			}
			result = columns + "," + rows;
			
		}
		return result;
	}

}
