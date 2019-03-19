package edu.mda.bioinfo.boxplotjava;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.GregorianCalendar;

/**
 *
 * @author tdcasasent
 */
public class BoxplotJava
{
	public static String getVersion()
	{
		return "BoxplotJava 2018-07-31-1330";
	}

	public static SimpleDateFormat mLogDate = new SimpleDateFormat("yyyy_MM_dd_HHmm.ss.SSS");

	public static void log(String theString)
	{
		System.out.println(" " + mLogDate.format(GregorianCalendar.getInstance().getTime()) + "\t" + "DEBUG" + "\t" + Thread.currentThread().getId() + "\t" + theString);
		System.out.flush();
	}

	public static void logData(double[] theMatrixData, String[] theBatchData,
			ReadMatrixFile rmf, ReadBatchFile rbf)
	{
		if (null!=theMatrixData)
		{
			BoxplotJava.log("theMatrixData.length=" + theMatrixData.length);
			System.out.print("theMatrixData=");
			for (int y = 0; y < 10; y++)
			{
				System.out.print(theMatrixData[y] + ", ");
			}
			System.out.println("-");
		}
		for (int x = 0; x < 2; x++)
		{
			System.out.print("x=" + x + " --> ");
			for (int y = 0; y < 10; y++)
			{
				System.out.print(rmf.mCombinedData[x][y] + ", ");
			}
			System.out.println("-");
		}
		if (null!=theBatchData)
		{
			System.out.println("theBatchData.length=" + theBatchData.length);
			System.out.print("theBatchData=");
			for (int y = 0; y < 10; y++)
			{
				System.out.print(theBatchData[y] + ", ");
			}
			System.out.println("-");
		}
		for (int x = 0; x < rbf.mBatches.length; x++)
		{
			System.out.print("x=" + x + " --> ");
			for (int y = 0; y < 10; y++)
			{
				System.out.print(rbf.mBatches[x][y] + ", ");
			}
			System.out.println("-");
		}
	}

	public static boolean allSampleData(double[] theMatrixData, String[] theBarcodes, String[] theGenes,
			String[] theBatchData, String[] theBatchTypes, String theOutputDir, String theTitle)
	{
		boolean success = false;
		try
		{
			BoxplotJava.log(getVersion());
			BoxplotJava.log("Prepare data for Boxplot_AllSampleData");
			BoxplotJava.log("theOutputDir="+theOutputDir);
			BoxplotJava.log("theTitle="+theTitle);
			ReadMatrixFile rmf = new ReadMatrixFile(theMatrixData, theBarcodes, theGenes);
			ReadBatchFile rbf = new ReadBatchFile(theBatchData, theBatchTypes, theBarcodes);
			//logData(theMatrixData, theBatchData, rmf, rbf);
			BoxplotJava.log("Before Boxplot_AllSampleData");
			try
			{
				Boxplot_AllSampleData.run(rmf, rbf, new File(theOutputDir), theTitle);
			}
			catch(Exception ex)
			{
				if (!ex.getMessage().startsWith("Excessive outliers."))
				{
					throw ex;
				}
			}

			BoxplotJava.log("After Boxplot_AllSampleData");
			success = true;
		}
		catch (Exception exp)
		{
			exp.printStackTrace(System.err);
		}
		return success;
	}

	public static boolean allSampleRLE(double[] theMatrixData, String[] theBarcodes, String[] theGenes,
			String[] theBatchData, String[] theBatchTypes, String theOutputDir, String theTitle)
	{
		boolean success = false;
		try
		{
			BoxplotJava.log(getVersion());
			BoxplotJava.log("Prepare data for Boxplot_AllSampleRLE");
			BoxplotJava.log("theOutputDir="+theOutputDir);
			BoxplotJava.log("theTitle="+theTitle);
			BoxplotJava.log("theBarcodes.length="+theBarcodes.length);
			BoxplotJava.log("theGenes.length="+theGenes.length);
			BoxplotJava.log("theMatrixData.length="+theMatrixData.length);
			BoxplotJava.log("theBatchTypes.length="+theBatchTypes.length);
			BoxplotJava.log("theBatchData.length="+theBatchData.length);
			ReadMatrixFile rmf = new ReadMatrixFile(theMatrixData, theBarcodes, theGenes);
			BoxplotJava.log("rmf.mGeneEqMap.size()="+rmf.mGeneEqMap.size());
			ReadBatchFile rbf = new ReadBatchFile(theBatchData, theBatchTypes, theBarcodes);
			BoxplotJava.log("rbf.mBatchTypes.length="+rbf.mBatchTypes.length);
			
			BoxplotJava.log("Before Boxplot_AllSampleRLE");
			try
			{
				Boxplot_AllSampleRLE.run(rmf, rbf, new File(theOutputDir), theTitle);
			}
			catch(Exception ex)
			{
				if (!ex.getMessage().startsWith("Excessive outliers."))
				{
					throw ex;
				}
			}
			BoxplotJava.log("After Boxplot_AllSampleRLE");
			success = true;
		}
		catch (Exception exp)
		{
			exp.printStackTrace(System.err);
		}
		return success;
	}

	public static boolean groupFunction(double[] theMatrixData, String[] theBarcodes, String[] theGenes,
			String[] theBatchData, String[] theBatchTypes, String theOutputDir, int theId, String theTitle)
	{
		// 2 is mean
		boolean success = false;
		try
		{
			BoxplotJava.log(getVersion());
			BoxplotJava.log("Prepare data for Boxplot_GroupFunction");
			BoxplotJava.log("theOutputDir="+theOutputDir);
			BoxplotJava.log("theTitle="+theTitle);
			ReadMatrixFile rmf = new ReadMatrixFile(theMatrixData, theBarcodes, theGenes);
			ReadBatchFile rbf = new ReadBatchFile(theBatchData, theBatchTypes, theBarcodes);
			GroupFunctions gf = GroupFunctions.getGroupFunctionFromId(theId);
			if (null == gf)
			{
				throw new Exception("Id '" + theId + "' not recognized for Group Function");
			}
			BoxplotJava.log("Before Boxplot_GroupFunction");
			try
			{
				Boxplot_GroupFunction.run(rmf, rbf, new File(theOutputDir), gf, theTitle);
			}
			catch(Exception ex)
			{
				if (!ex.getMessage().startsWith("Excessive outliers."))
				{
					throw ex;
				}
			}
			BoxplotJava.log("After Boxplot_GroupFunction");
			success = true;
		}
		catch (Exception exp)
		{
			exp.printStackTrace(System.err);
		}
		return success;
	}

	public static String[] readArrayString(String theFile) throws IOException
	{
		String line = java.nio.file.Files.readAllLines(Paths.get(theFile), Charset.availableCharsets().get("ISO-8859-1")).get(0);
		return line.split("\t", -1);
	}

	public static double[] readArrayDouble(String theFile) throws IOException
	{
		String[] line = readArrayString(theFile);
		double[] result = new double[line.length];
		for (int x = 0; x < line.length; x++)
		{
			if ("NA".equalsIgnoreCase(line[x]))
			{
				result[x] = Double.NaN;
			}
			else
			{
				result[x] = Double.parseDouble(line[x]);
			}
		}
		return result;
	}

	public static void mainSimulatedRData(String[] args)
	{
		//mainSimulatedRData
		try
		{
			System.out.println(getVersion());
			String baseDir = args[0];
			System.out.println("theBarcodes");
			String[] barcodes = readArrayString(new File(baseDir, "theBarcodes.txt").getAbsolutePath());
			System.out.println("theBatchData");
			String[] batchData = readArrayString(new File(baseDir, "theBatchData.txt").getAbsolutePath());
			System.out.println("theBatchTypes");
			String[] batchTypes = readArrayString(new File(baseDir, "theBatchTypes.txt").getAbsolutePath());
			System.out.println("theGenes");
			String[] genes = readArrayString(new File(baseDir, "theGenes.txt").getAbsolutePath());
			System.out.println("theMatrixData");
			double[] matrixData = readArrayDouble(new File(baseDir, "theMatrixData.txt").getAbsolutePath());
			System.out.println("output");
			File outputDir = new File(new File(baseDir, "output"), "BoxPlot");
			System.out.println("allSampleData");
			String title = "All Sample Data";
			allSampleData(matrixData, barcodes, genes, batchData, batchTypes, outputDir.getAbsolutePath(), title);
		}
		catch (IOException exp)
		{
			exp.printStackTrace(System.err);
		}
	}


	/**
	 *
	 * @param args arg[0] = base directory from which to read data and write output
	 */
//	public static void mainFile(String[] args)
	public static void mainFile(String[] args) throws IOException, Exception
	{
		System.out.println(getVersion());
		String dataDir = args[0];
		File outputDir = new File(args[1]);

		File matrixData = new File(dataDir, "matrix_data.tsv");
		File batchData = new File(dataDir, "batches.tsv");
		ReadMatrixFile rmf = new ReadMatrixFile(matrixData);
		rmf.readFile();
		ReadBatchFile rbf = new ReadBatchFile(batchData, rmf.mSampleList);
		rbf.readFile();
		System.out.println("-------------------------------------------");
		System.out.println("---- Boxplot_AllSampleData");
		System.out.println("-------------------------------------------");
		//logData(null, null, rmf, rbf);
		Boxplot_AllSampleData.run(rmf, rbf, outputDir, "Boxplot All Sample Data");
		System.out.println("-------------------------------------------");
		System.out.println("---- Boxplot_AllSampleRLE");
		System.out.println("-------------------------------------------");
		Boxplot_AllSampleRLE.run(rmf, rbf, outputDir, "Boxplot All Sample RLE");
		System.out.println("-------------------------------------------");
		System.out.println("---- Boxplot_GroupFunction");
		System.out.println("-------------------------------------------");
		Boxplot_GroupFunction.run(rmf, rbf, outputDir, GroupFunctions.MEAN, "Boxplot Group Function");
		System.out.println("-------------------------------------------");
		System.out.println("-------------------------------------------");
	}
}
