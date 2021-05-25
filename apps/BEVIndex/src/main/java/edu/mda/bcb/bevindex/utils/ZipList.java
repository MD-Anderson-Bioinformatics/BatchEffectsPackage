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

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.TreeSet;
import static java.util.stream.Collectors.toList;
import java.util.stream.Stream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

/**
 *
 * @author Tod-Casasent
 */
public class ZipList
{
	private TreeSet<File> mIncludeFiles = null;
	private File mBaseDir = null;
	public ZipList(File theBaseDir)
	{
		mBaseDir = theBaseDir;
		mIncludeFiles = new TreeSet<>();
	}
	
	public void lsInclude() throws IOException
	{
		lsInclude(mBaseDir.toPath());
	}
	
	private void lsInclude(Path thePath) throws IOException
	{
		System.out.println("ZipList::lsInclude check path: " + thePath.toString());
		try (Stream<Path> myFiles = Files.list(thePath))
		{
			for (Path file : myFiles.filter(Files::isRegularFile).collect(toList()))
			{
				System.out.println("ZipList::lsInclude add file: " + file.toString());
				mIncludeFiles.add(file.toFile());
			}
		}
		try (Stream<Path> myFiles = Files.list(thePath))
		{
			for (Path dir : myFiles.filter(Files::isDirectory).collect(toList()))
			{
				lsInclude(dir);
			}
		}
	}
	
	public void zipList(Path theZip) throws FileNotFoundException, IOException
	{
		// remove ZIP if it exists
		if (theZip.toFile().exists())
		{
			System.out.println("Remove old ZIP");
			theZip.toFile().delete();
		}
		// then create ZIP
		ZipOutputStream zip = null;
		try
		{
			System.out.println("ZipList::zipList zip: " + theZip.toString());
			zip = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(theZip.toFile())));
			for (File file : mIncludeFiles)
			{
				addEntry(zip, file, file.getAbsolutePath().replaceFirst(mBaseDir.getAbsolutePath(), ""));
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
	
	protected void addEntry(ZipOutputStream theZip, File theFile, String theEntryName) throws IOException
	{
		if (theEntryName.startsWith("/"))
		{
			theEntryName = theEntryName.substring(1);
		}
		if (theEntryName.startsWith("\\"))
		{
			theEntryName = theEntryName.substring(1);
		}
		System.out.println("ZipList::addEntry theEntryName: " + theEntryName + " for " + theFile.getAbsolutePath());
		// add file
		// write the file to the output
		byte[] buf = new byte[1024];
		int len;
		System.out.println("Add entry " + theFile + " as " + theEntryName);
		theZip.putNextEntry(new ZipEntry(theEntryName));
		try(FileInputStream in = new FileInputStream(theFile))
		{
			while ((len = in.read(buf)) > 0)
			{
				// Write the Result
				theZip.write(buf, 0, len);
			}
		}
	}
}
