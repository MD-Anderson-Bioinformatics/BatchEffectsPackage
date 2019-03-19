/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.mda.bioinfo.boxplotjava;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author tdcasasent
 */
public class ReadBatchFile
{

	public ArrayList<String> mMatrixSampleList = null;
	public String[] mBatchTypes = null;
	// mBatches[batchTypeIndex][barcodeIndex]
	public String[][] mBatches = null;
	//
	public File mBatchFile = null;

	public ReadBatchFile(File theBatchFile, ArrayList<String> theMatrixSampleList) throws FileNotFoundException, IOException
	{
		mBatchFile = theBatchFile;
		// subtract 1, for Type (Sample column was already subtracted)
		int batchTypeCount = ReadMatrixFile.getColumnCount(mBatchFile);
		int sampleCount = theMatrixSampleList.size();
		mBatchTypes = new String[batchTypeCount];
		mBatches = new String[batchTypeCount][sampleCount];
		mMatrixSampleList = theMatrixSampleList;
	}

	public ReadBatchFile(String[] theBatchData, String[] theBatchTypes, String[] theBarcodes)
	{
		mBatchFile = null;
		int batchTypeCount = theBatchTypes.length;
		int sampleCount = theBarcodes.length;
		mBatchTypes = theBatchTypes;
		mBatches = convertToMatrix(theBatchData, batchTypeCount, sampleCount);
		mMatrixSampleList = new ArrayList<>();
		mMatrixSampleList.addAll(Arrays.asList(theBarcodes));
	}

	static protected String[][] convertToMatrix(String[] theValues, int theBatchTypeCount, int theBarcodeCount)
	{
		String[][] matrix = new String[theBatchTypeCount][theBarcodeCount];
		int index = 0;
		for (int x = 0; x < theBatchTypeCount; x++)
		{
			for (int y = 0; y < theBarcodeCount; y++)
			{
				matrix[x][y] = theValues[index];
				index = index + 1;
			}
		}
		return matrix;
	}

	public void readFile() throws IOException, Exception
	{
		BoxplotJava.log("ReadBatchFile::process - Start");
		try (BufferedReader br = Files.newBufferedReader(
						Paths.get(mBatchFile.getAbsolutePath()),
						Charset.availableCharsets().get("ISO-8859-1")))
		{
			// first line is batch types
			String line = br.readLine();
			{
				String[] splitted = line.split("\t");
				for (int x = 1; x < splitted.length; x++)
				{
					mBatchTypes[x - 1] = splitted[x];
				}
			}
			// do rest
			line = br.readLine();
			// mBatches[batchTypeIndex][barcodeIndex]
			BoxplotJava.log("ReadBatchFile::process - before lines");
			while (null != line)
			{
				String[] splitted = line.split("\t");
				int sampleIndex = mMatrixSampleList.indexOf(splitted[0]);
				// skip 0=barcode, 1=type
				for (int x = 1; x < splitted.length; x++)
				{
					// mBatches[batchTypeIndex][barcodeIndex]
					mBatches[x - 1][sampleIndex] = splitted[x];
				}
				//sampleIndex = sampleIndex + 1;
				line = br.readLine();
			}
		}
		BoxplotJava.log("ReadBatchFile::process - Finish");
	}
}
