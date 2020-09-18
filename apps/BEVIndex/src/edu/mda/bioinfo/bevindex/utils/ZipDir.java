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

package edu.mda.bioinfo.bevindex.utils;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Path;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

/**
 *
 * @author Tod-Casasent
 */
public class ZipDir
{
	public ZipDir()
	{
		
	}
	
	public void zipContents(Path theDir, Path theZip) throws FileNotFoundException, IOException
	{
		// remove ZIP if it exists
		if (theZip.toFile().exists())
		{
			System.out.println("Remove old ZIP");
			theZip.toFile().delete();
		}
		// first collect contents to add (to prevent ZIP being self-referential
		File [] dirList = theDir.toFile().listFiles();
		// then create ZIP
		ZipOutputStream zip = null;
		try
		{
			zip = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(theZip.toFile())));
			// top level entry
			zip.putNextEntry(new ZipEntry(theDir.toFile().getName() + File.separator));
			// and then contents to top dirs
			for (File entry : dirList)
			{
				addEntry(zip, theDir.toFile().getName() + File.separator, entry);
			}
		}
		finally
		{
			if (null != zip)
			{
				zip.flush();
				zip.close();
			}
		}
	}
	
	protected void addEntry(ZipOutputStream theZip, String thePath, File theEntry) throws IOException
	{
		if (theEntry.isDirectory())
		{
			System.out.println("Add directory " + theEntry);
			String newPath = thePath + theEntry.getName() + File.separator;
			// add directory
			theZip.putNextEntry(new ZipEntry(newPath));
			// add contents
			File [] dirList = theEntry.listFiles();
			for(File entry : dirList)
			{
				addEntry(theZip, newPath, entry);
			}
		}
		else
		{
			//System.out.println("Add file " + theEntry);
			// add file
			// write the file to the output
			byte[] buf = new byte[1024];
			int len;
			if (theEntry.getName().equals("index.json"))
			{
				theZip.putNextEntry(new ZipEntry(theEntry.getName()));
			}
			else
			{
				theZip.putNextEntry(new ZipEntry(thePath + theEntry.getName()));
			}
			try(FileInputStream in = new FileInputStream(theEntry))
			{
				while ((len = in.read(buf)) > 0)
				{
					// Write the Result
					theZip.write(buf, 0, len);
				}
			}
		}
	}
	
}
