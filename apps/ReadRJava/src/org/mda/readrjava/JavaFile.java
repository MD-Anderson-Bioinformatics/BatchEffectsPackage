/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.mda.readrjava;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.TreeSet;

/**
 *
 * @author tdcasasent
 */
public class JavaFile
{
	public static String M_DELIMITER = "\t";

	public String mFile = null;
	public ArrayList<String> mColumns = null;
	public ArrayList<String> mRows = null;
	public ArrayList<Double> mDoubleData = null;
	public ArrayList<String> mStringData = null;
	public int mHeaderCount = -1;
	public int mRowColumnCount = -1;
	public boolean mSkipFirstLineFlag = false;
	public boolean mPadDataFlag = false;

	public void addStringData(HashMap<String, String> theNewData) throws Exception
	{
		for(String col : mColumns)
		{
			String val = theNewData.get(col);
			if (null==val)
			{
				throw new Exception("New data missing value for column " + col);
			}
			else
			{
				mStringData.add(val);
			}
		}
	}

	public String getStringDataByColumn(int theRow, String theColumn) throws Exception
	{
		int index = mColumns.indexOf(theColumn);
		if(-1==index)
		{
			throw new Exception("Column '" + theColumn + "' not found");
		}
		//System.out.println("length=" + mStringData.size() + " row=" + theRow + " col=" + theColumn + " index=" + index + " math=" + (mColumns.size()*theRow) + " total=" + ((mColumns.size()*theRow)+index));
		//System.out.flush();
		return(mStringData.get((mHeaderCount*theRow)+index));
	}

	public JavaFile(String [] theColumns, boolean theCols, boolean theRows, boolean theData, boolean theSkipFirstLineFlag, boolean thePadDataFlag)
	{
		mPadDataFlag = thePadDataFlag;
		mSkipFirstLineFlag = theSkipFirstLineFlag;
		if (true==theCols)
		{
			mColumns = new ArrayList<>();
			mColumns.addAll(Arrays.asList(theColumns));
		}
		if (true==theRows)
		{
			mRows = new ArrayList<>();
		}
		if (true==theData)
		{
			mDoubleData = new ArrayList<>();
			mStringData = new ArrayList<>();
		}
		mHeaderCount = mColumns.size();
	}

	public JavaFile(String theFile, boolean theCols, boolean theRows, boolean theData, boolean theSkipFirstLineFlag, boolean thePadDataFlag)
	{
		mPadDataFlag = thePadDataFlag;
		mSkipFirstLineFlag = theSkipFirstLineFlag;
		mFile = theFile;
		if (true==theCols)
		{
			mColumns = new ArrayList<>();
		}
		if (true==theRows)
		{
			mRows = new ArrayList<>();
		}
		if (true==theData)
		{
			mDoubleData = new ArrayList<>();
			mStringData = new ArrayList<>();
		}
	};

	protected void checkAndIncrementRowCount(int theCount) throws Exception
	{
		if (-1==mRowColumnCount)
		{
			mRowColumnCount = theCount;
		}
		else
		{
			if (mRowColumnCount!=theCount)
			{
				throw new Exception("New row count (" + theCount + ") does not equal old row count (" + mRowColumnCount + ")");
			}
		}
	}

	public void loadDoubleData() throws IOException, Exception
	{
		boolean firstString = true;
		try(BufferedReader br = Files.newBufferedReader(Paths.get(mFile), ReadRJava.getCharset()))
		{
			//int cnt = 0;
			String line = br.readLine();
			if (true==mSkipFirstLineFlag)
			{
				line = br.readLine();
			}
			while(null!=line)
			{
				if (false==line.startsWith("#"))
				{
					ArrayList<String> strList = new ArrayList<>(Arrays.asList(line.split(M_DELIMITER, -1)));
					if ((true==firstString)&&(null!=mColumns))
					{
						// process first line
						if (null!=mRows)
						{
							strList.remove(0);
						}
						mColumns.addAll(strList);
						mHeaderCount = mColumns.size();
						firstString = false;
					}
					else
					{
						// process regular line
						if (null!=mRows)
						{
							mRows.add(strList.get(0));
							strList.remove(0);
							if (null!=mDoubleData)
							{
								checkAndIncrementRowCount(strList.size());
								for(String dat : strList)
								{
									if ("NA".equals(dat))
									{
										mDoubleData.add(Double.NaN);
									}
									else
									{
										mDoubleData.add(new Double(dat));
									}
								}
							}
						}
						else if (null!=mDoubleData)
						{
							checkAndIncrementRowCount(strList.size());
							for(String dat : strList)
							{
								if ("NA".equals(dat))
								{
									mDoubleData.add(Double.NaN);
								}
								else
								{
									mDoubleData.add(new Double(dat));
								}
							}
						}
					}
				}
				line = br.readLine();
				//cnt = cnt + 1;
				//if (0==(cnt % 100))
				//{
				//	System.out.println("lines processed " + cnt);
				//}
			}
		}
	}

	public void loadStringData() throws IOException, Exception
	{
		boolean firstString = true;
		try(BufferedReader br = Files.newBufferedReader(Paths.get(mFile), ReadRJava.getCharset()))
		{
			mRows = null;
			String line = br.readLine();
			if (true==mSkipFirstLineFlag)
			{
				//System.out.println("loadStringData - skip first line");
				line = br.readLine();
			}
			while(null!=line)
			{
				if (false==line.startsWith("#"))
				{
					ArrayList<String> strList = new ArrayList<>(Arrays.asList(line.split(M_DELIMITER, -1)));
					if ((true==firstString)&&(null!=mColumns))
					{
						//System.out.println("loadStringData - column headers");
						//System.out.println("loadStringData - process first line");
						//System.out.println("loadStringData - if first entry is empty, then skip it and set read rows flag");
						if ("".equals(strList.get(0)))
						{
							mRows = new ArrayList<>();
							mColumns.addAll(strList.subList(1, strList.size()));
						}
						else
						{
							mColumns.addAll(strList);
						}
						mHeaderCount = mColumns.size();
						firstString = false;
					}
					else
					{
						//System.out.println("loadStringData - data line");
						// process regular line
						if (null!=mRows)
						{
							//System.out.println("loadStringData - Row Header Data");
							mRows.add(strList.get(0));
							strList.remove(0);
							if (null!=mStringData)
							{
								checkAndIncrementRowCount(strList.size());
								for(String dat : strList)
								{
									mStringData.add(dat);
								}
							}
						}
						else if (null!=mStringData)
						{
							//System.out.println("loadStringData - String Data");
							if (strList.size()<mHeaderCount)
							{
								//System.out.println("WARN-READ_PAD-DATA - padding string data for " + mFile);
								while(strList.size()<mHeaderCount)
								{
									strList.add("");
								}
							}
							checkAndIncrementRowCount(strList.size());
							//System.out.println("loadStringData - no row line");
							for(String dat : strList)
							{
								mStringData.add(dat);
							}
						}
					}
				}
				else
				{
					//System.out.println("loadStringData - skip comment");
				}
				//System.out.println("loadStringData - next line");
				line = br.readLine();
				//cnt = cnt + 1;
				//if (0==(cnt % 100))
				//{
				//	System.out.println("lines processed " + cnt);
				//}
			}
		}
	}

	public int getmHeaderCount()
	{
		return mHeaderCount;
	}

	public int getmRowCount()
	{
		return mRowColumnCount;
	}

	public String [] getmColumns()
	{
		if (null==mColumns)
		{
			return null;
		}
		else
		{
			return mColumns.toArray(new String[0]);
		}
	}

	public String [] getmRows()
	{
		if (null==mRows)
		{
			return null;
		}
		else
		{
			return mRows.toArray(new String[0]);
		}
	}

	public double [] getmDoubleData()
	{
		double [] result = null;
		if (null==mRows)
		{
			result = null;
		}
		else
		{
			result = new double[mDoubleData.size()];
			for (int x=0; x<mDoubleData.size(); x++)
			{
				result[x] = mDoubleData.get(x).doubleValue();
			}
		}
		return result;
	}

	public String [] getmStringData()
	{
		if (null==mStringData)
		{
			return null;
		}
		else
		{
			return mStringData.toArray(new String[0]);
		}
	}

	public static boolean writeFile(String theFile, String [] theCols, String [] theRows,
			String [] theStringData, double [] theDoubleData) throws IOException
	{
		System.out.println("writeFile - start");
		// theCols is required, to determine size for arraylist
		boolean result = false;
		try(BufferedWriter bw = Files.newBufferedWriter(Paths.get(theFile), ReadRJava.getCharset()))
		{
			//////////////////////////////////////////
			//System.out.println("writeFile - write header line");
			//////////////////////////////////////////
			if (null!=theRows)
			{
				//System.out.println("writeFile - write header line with blank for rows column");
				bw.write("\t");
			}
			//System.out.println("writeFile - write header line without blank for rows column");
			boolean first = true;
			for(String col : theCols)
			{
				if (false==first)
				{
					bw.write("\t");
				}
				else
				{
					first = false;
				}
				bw.write(col);
			}
			bw.newLine();
			//////////////////////////////////////////
			//System.out.println("writeFile - write data lines");
			//////////////////////////////////////////
			int dataSize = 0;
			if (null!=theDoubleData)
			{
				dataSize = theDoubleData.length;
			}
			else
			{
				dataSize = theStringData.length;
			}
			int indexData = 0;
			int indexRow = 0;
			while(indexData<dataSize)
			{
				// write a data line - row name if provided
				if (null!=theRows)
				{
					//System.out.println("writeFile - write a data line - row name if provided");
					bw.write(theRows[indexRow]);
					bw.write("\t");
				}
				// write a data line - data
				first = true;
				for(int x=0;x<theCols.length;x++)
				{
					if (false==first)
					{
						bw.write("\t");
					}
					else
					{
						first = false;
					}
					if (null!=theDoubleData)
					{
						bw.write(Double.toString(theDoubleData[indexData]));
					}
					else
					{
						bw.write(theStringData[indexData]);
					}
					indexData = indexData + 1;
				}
				bw.newLine();
				indexRow = indexRow + 1;
			}
			result = true;
		}
		System.out.println("writeFile - done");
		return result;
	}

	public void checkNames()
	{
		mColumns = checkNamesList(mColumns);
	}

	static public ArrayList<String> checkNamesList(ArrayList<String> theList)
	{
		ArrayList<String> newColumns = new ArrayList<>();
		for(String name: theList)
		{
			newColumns.add(checkname_contents(name));
		}
		// for each duplicate value, append (to second and subsequent occurances) dot N where N starts at 1
		TreeSet<String> tree = new TreeSet<>(newColumns);
		theList = newColumns;
		for(String name: tree)
		{
			if (Collections.frequency(newColumns, name)>1)
			{
				newColumns = new ArrayList<>();
				int count = 0;
				for(String testName : theList)
				{
					if (name.equals(testName))
					{
						if (count>0)
						{
							testName = testName + "." + count;
						}
						count = count + 1;
					}
					newColumns.add(testName);
				}
				theList = newColumns;
			}
		}
		return(theList);
	}

	static public String checkname_contents(String theName)
	{
		// duplicates R make.names functionality
		// only contains letters, numbers, dot, underline -- replace everything else with dot
		theName = theName.replaceAll("[^a-zA-Z0-9_]", ".");
		// if name starts with dot number, prepend with "X"
		if (theName.matches("^\\.\\d"))
		{
			theName = "X" + theName;
		}
		return(theName);
	}

	public String getFirstSortedValue(String theColumn) throws Exception
	{
		TreeSet<String> uniqueOnlyList = new TreeSet<>();
		for(int index=0;index<(this.mStringData.size()/this.mHeaderCount);index++)
		{
			String col = this.getStringDataByColumn(index, theColumn);
			uniqueOnlyList.add(col);
		}
		return(uniqueOnlyList.first());
	}

	public String getLastSortedValue(String theColumn) throws Exception
	{
		TreeSet<String> uniqueOnlyList = new TreeSet<>();
		for(int index=0;index<(this.mStringData.size()/this.mHeaderCount);index++)
		{
			String col = this.getStringDataByColumn(index, theColumn);
			uniqueOnlyList.add(col);
		}
		return(uniqueOnlyList.last());
	}
}
