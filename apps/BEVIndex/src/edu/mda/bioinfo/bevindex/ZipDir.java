/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.mda.bioinfo.bevindex;

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
 * @author linux
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
