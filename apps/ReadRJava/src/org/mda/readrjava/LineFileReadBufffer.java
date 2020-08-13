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

package org.mda.readrjava;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 *
 * @author tdcasasent
 */
public class LineFileReadBufffer implements AutoCloseable
{
	protected ArrayList<String> mHeaders = null;
	protected String mFile = null;
	protected BufferedReader mReader = null;
	protected boolean mSkipFirstLineFlag = false;
	protected boolean mSkipSecondLineFlag = false;
	protected boolean mChecknamesFlag = false;
	protected boolean mIsWGBSFlag = false;

	public LineFileReadBufffer(String theFile, boolean theSkipFirstLineFlag, boolean theSkipSecondLineFlag,
			boolean theChecknamesFlag, boolean theIsWGBSFlag)
	{
		mHeaders = null;
		mFile = theFile;
		mSkipFirstLineFlag = theSkipFirstLineFlag;
		mSkipSecondLineFlag = theSkipSecondLineFlag;
		mChecknamesFlag = theChecknamesFlag;
		mIsWGBSFlag = theIsWGBSFlag;
	}

	public void init() throws IOException
	{
		mHeaders = new ArrayList<>();
		long rowValueCount = 0;
		boolean firstString = true;
		try(BufferedReader br = Files.newBufferedReader(Paths.get(mFile), ReadRJava.getCharset()))
		{
			String line = br.readLine();
			if (true==mSkipFirstLineFlag)
			{
				line = br.readLine();
			}
			while(null!=line)
			{
				if (false==line.startsWith("#"))
				{
					ArrayList<String> strList = null;
					if ((true==firstString)&&(true==mIsWGBSFlag))
					{
						// this is a WGBS with fields="chrom,chromStart,chromEnd,name,score,strand,percentMeth,numCTreads" at the end
						// extract the "fields" value and part that for column headers
						// find "fields" index
						int fieldIndex = line.indexOf("fields");
						// find first " after fields
						int firstQuote = line.indexOf("\"", fieldIndex);
						// find next " after taht quote
						int nextQuote = line.indexOf("\"", firstQuote+1);
						// get string in beween
						String substring = line.substring(firstQuote+1, nextQuote);
						// split on commas
						String [] columnHeaders = substring.split(",", -1);
						// copy
						strList = new ArrayList<>(Arrays.asList(columnHeaders));
					}
					else
					{
						strList = new ArrayList<>(Arrays.asList(line.split("\t", -1)));
					}
					if (true==firstString)
					{
						if (mChecknamesFlag)
						{
							strList = JavaFile.checkNamesList(strList);
						}
						mHeaders.addAll(strList);
						rowValueCount = mHeaders.size();
						firstString = false;
						if (true==mSkipSecondLineFlag)
						{
							line = br.readLine();
						}
					}
					else
					{
						if (strList.size()>rowValueCount)
						{
							rowValueCount = strList.size();
						}
					}
				}
				line = br.readLine();
			}
		}
		if (rowValueCount>mHeaders.size())
		{
			for(int x=0;x<(rowValueCount-mHeaders.size());x++)
			{
				mHeaders.add("unnamed_" +  String.format("%07d", x));
			}
		}
	}

	public HashMap<String, String> getNextRow() throws IOException
	{
		HashMap<String, String> nextRow = null;
		String line = null;
		if (null==mHeaders)
		{
			init();
			mReader = Files.newBufferedReader(Paths.get(mFile), ReadRJava.getCharset());
			// skip first line if requested
			if (true==mSkipFirstLineFlag)
			{
				mReader.readLine();
			}
			if (true==mSkipSecondLineFlag)
			{
				mReader.readLine();
			}
			// skip comment lines if present
			line = mReader.readLine();
			while (line.startsWith("#"))
			{
				line = mReader.readLine();
			}
			// skip headers
			//mReader.readLine();
		}
		line = mReader.readLine();
		// skip comment lines if present
		while ((null!=line)&&(line.startsWith("#")))
		{
			line = mReader.readLine();
		}
		if (null!=line)
		{
			nextRow = new HashMap<>();
			String [] strList = line.split("\t", -1);
			for(int x=0;x<mHeaders.size();x++)
			{
				String header = mHeaders.get(x);
				if (x<strList.length)
				{
					nextRow.put(header, strList[x]);
				}
				else
				{
					nextRow.put(header, "");
				}
			}
		}
		return(nextRow);
	}

	@Override
	public void close() throws Exception
	{
		mReader.close();
	}


}
