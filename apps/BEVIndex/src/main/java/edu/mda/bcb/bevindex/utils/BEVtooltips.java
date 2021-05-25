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

package edu.mda.bcb.bevindex.utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

/**
 *
 * @author Tod-Casasent
 */
public class BEVtooltips
{

	protected HashMap<String, HashMap<String, String>> mTooltipMap = new HashMap<>();

	public BEVtooltips()
	{
		mTooltipMap = new HashMap<>();
	}

	public void readTsv(String thePath) throws FileNotFoundException, IOException, Exception
	{
		mTooltipMap = new HashMap<>();
		if (null!=thePath)
		{
			ArrayList<String> headers = new ArrayList<>();
			try (Scanner scanner = new Scanner(new File(thePath)))
			{
				while (scanner.hasNext())
				{
					String line = scanner.nextLine();
					String [] splitted = line.split("\t", -1);
					if (0==headers.size())
					{
						for (int x=0;x<splitted.length;x++)
						{
							headers.add(splitted[x]);
						}
					}
					else
					{
						String header_mLabel = splitted[headers.indexOf("mLabel")];
						String header_mName = splitted[headers.indexOf("mName")];
						String header_mTooltip = splitted[headers.indexOf("mTooltip")];
						addTooltip(header_mLabel, header_mName, header_mTooltip);
					}
				}
			}
		}
	}
	
	protected void addTooltip(String theLabel, String theName, String theTooltip) throws Exception
	{
		HashMap<String, String> nameToTooltip = mTooltipMap.get(theLabel);
		if (null==nameToTooltip)
		{
			nameToTooltip = new HashMap<>();
		}
		String tooltip = nameToTooltip.get(theName);
		if (null!=tooltip)
		{
			throw new Exception("Found duplicate tooltip entries for same label/name pair "+ theLabel + " and " + theName );
		}
		else
		{
			nameToTooltip.put(theName, theTooltip);
			//nameToTooltip.put(theName, StringEscapeUtils.escapeJson(theTooltip));
			mTooltipMap.put(theLabel, nameToTooltip);
		}
	}
	
	public String getTooltip(String theLabel, String theName)
	{
		String result = "";
		HashMap<String, String> nameToTooltip = mTooltipMap.get(theLabel);
		if (null!=nameToTooltip)
		{
			result = nameToTooltip.get(theName);
			if (null==result)
			{
				result = nameToTooltip.get("*");
			}
			if (null==result)
			{
				result = "";
			}
		}
		return result;
	}
}
