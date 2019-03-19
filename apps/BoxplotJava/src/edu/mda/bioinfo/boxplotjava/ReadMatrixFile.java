/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.mda.bioinfo.boxplotjava;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;

/**
 *
 * @author tdcasasent
 */
public class ReadMatrixFile
{

	protected int M_BARCODES;
	protected int M_GENES;
	public double[][] mCombinedData;
	public HashMap<String, Integer> mGeneEqMap = new HashMap<>();
	public ArrayList<String> mSampleList = new ArrayList<>();
	public int mCurrentColumn = 0;
	//
	public File mMatrixFile = null;

	public ReadMatrixFile(double[] theValues, String[] theBarcodes, String[] theGenes)
	{
		mMatrixFile = null;
		M_BARCODES = theBarcodes.length;
		M_GENES = theGenes.length;
		mCombinedData = convertToMatrix(theValues, M_BARCODES, M_GENES);
		mCurrentColumn = 0;
		mGeneEqMap = new HashMap<>();
		for (int x = 0; x < theGenes.length; x++)
		{
			mGeneEqMap.put(theGenes[x], x);
		}
		mSampleList = new ArrayList<>();
		mSampleList.addAll(Arrays.asList(theBarcodes));
	}

	public ReadMatrixFile(File theMatrixFile) throws FileNotFoundException, IOException
	{
		mMatrixFile = theMatrixFile;
		M_BARCODES = getColumnCount(theMatrixFile);
		M_GENES = getRowCount(theMatrixFile);
		mCombinedData = new double[M_BARCODES][M_GENES];
		mCurrentColumn = 0;
		for (int x = 0; x < M_BARCODES; x++)
		{
			for (int y = 0; y < M_GENES; y++)
			{
				mCombinedData[x][y] = Double.NaN;
			}
		}
		mGeneEqMap = new HashMap<>();
		mSampleList = new ArrayList<>();
	}

	public String[] geneEqInIndexOrder()
	{
		String[] geneList = new String[mGeneEqMap.size()];
		for (Entry<String, Integer> entry : mGeneEqMap.entrySet())
		{
			geneList[entry.getValue().intValue()] = entry.getKey();
		}
		return geneList;
	}

	static protected double[][] convertToMatrix(double[] theValues, int theBarcodeLength, int theGeneLength)
	{
		BoxplotJava.log("ReadMatrixFile::convertToMatrix theBarcodeLength=" + theBarcodeLength);
		BoxplotJava.log("ReadMatrixFile::convertToMatrix theGeneLength=" + theGeneLength);
		double[][] matrix = new double[theBarcodeLength][theGeneLength];
		int index = 0;
		for (int x = 0; x < theBarcodeLength; x++)
		{
			for (int y = 0; y < theGeneLength; y++)
			{
				matrix[x][y] = theValues[index];
				index = index + 1;
			}
		}
		return matrix;
	}

	static public String beforeTab(String theString)
	{
		int index = theString.indexOf("\t");
		String subString = theString.substring(0, index);
		return subString;
	}

	static public String afterTab(String theString)
	{
		int index = theString.indexOf("\t");
		String subString = theString.substring(index + 1);
		return subString;
	}

	public void readFile() throws IOException, Exception
	{
		BoxplotJava.log("ReadMatrixFile::process - Start");
		try (BufferedReader br = Files.newBufferedReader(
						Paths.get(mMatrixFile.getAbsolutePath()),
						Charset.availableCharsets().get("ISO-8859-1")))
		{
			// first line samples
			String line = br.readLine();
			BoxplotJava.log("ReadMatrixFile::process - populateSampleLists");
			populateSampleLists(afterTab(line).split("\t", -1));
			// do rest
			int nextIndex = 0;
			line = br.readLine();
			int lineCnt = 0;
			BoxplotJava.log("ReadMatrixFile::process - before lines");
			while (null != line)
			{
				String geneEq = beforeTab(line);
				String data = afterTab(line);
				nextIndex = populateGeneAndData(geneEq, data.split("\t", -1), mCurrentColumn);
				line = br.readLine();
				lineCnt = lineCnt + 1;
				if (0 == (lineCnt % 1000))
				{
					System.out.print(" " + lineCnt);
				}
				if (0 == (lineCnt % 10000))
				{
					BoxplotJava.log("");
				}
			}
			BoxplotJava.log(" -");
			BoxplotJava.log("ReadMatrixFile::process - after lines");
			mCurrentColumn = nextIndex;
		}
		BoxplotJava.log("ReadMatrixFile::process - Finish");
	}

	protected int populateGeneAndData(String theGeneEq, String[] theData, int theStart)
	{
		//BoxplotJava.log("ReadMatrixFile::populateGeneAndData - Start");
		Integer intGE = mGeneEqMap.get(theGeneEq);
		if (null == intGE)
		{
			int newIndex = mGeneEqMap.size();
			mGeneEqMap.put(theGeneEq, newIndex);
			intGE = mGeneEqMap.get(theGeneEq);
		}
		int indexGE = intGE.intValue();
		for (String value : theData)
		{
			double dVal = Double.NaN;
			if (!"NA".equalsIgnoreCase(value))
			{
				if (!"".equals(value))
				{
					dVal = Double.parseDouble(value);
				}
			}
			mCombinedData[theStart][indexGE] = dVal;
			theStart = theStart + 1;
		}
		//BoxplotJava.log("ReadMatrixFile::populateGeneAndData - Finish");
		return theStart;
	}

	protected void populateSampleLists(String[] theSamples)
	{
		BoxplotJava.log("ProcessFile::populateSampleLists - Start");
		mSampleList.addAll(Arrays.asList(theSamples));
		BoxplotJava.log("ProcessFile::populateSampleLists - Finish");
	}

	static protected int getRowCount(File theMatrixFile) throws FileNotFoundException, IOException
	{
		BoxplotJava.log("Matrix::getLines start");
		int lines = 0;
		ArrayList<String> rows = new ArrayList<>();
		try (BufferedReader br = new BufferedReader(new FileReader(theMatrixFile)))
		{
			// skip header
			br.readLine();
			for (String line = br.readLine(); line != null; line = br.readLine())
			{
				lines = lines + 1;
			}
		}
		BoxplotJava.log("Matrix::getLines finished");
		return lines;
	}

	static protected int getColumnCount(File theMatrixFile) throws FileNotFoundException, IOException
	{
		BoxplotJava.log("Matrix::getColumns start");
		int columns;
		try (BufferedReader br = new BufferedReader(new FileReader(theMatrixFile)))
		{
			// skip first blank column for row headers
			columns = br.readLine().split("\t", -1).length - 1;
		}
		BoxplotJava.log("Matrix::getColumns finished");
		return columns;
	}
}
