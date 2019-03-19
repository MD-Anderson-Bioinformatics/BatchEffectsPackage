/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.mda.bioinfo.bevindex;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

/**
 *
 * @author linux
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
