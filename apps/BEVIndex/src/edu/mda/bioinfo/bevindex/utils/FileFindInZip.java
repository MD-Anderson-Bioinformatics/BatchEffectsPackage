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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Enumeration;
import java.util.TreeSet;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

/**
 *
 * @author linux
 */
public class FileFindInZip
{
	public TreeSet<Path> mMatches = null;
	public int mColumns = -1;
	public int mRows = -1;

	public FileFindInZip()
	{
		mMatches = new TreeSet<>();
	}

	public void find(Path theZipPath, String theFilename) throws IOException
	{
		if (theZipPath.toFile().exists())
		{
			try(ZipFile zf = new ZipFile(theZipPath.toFile()))
			{
				Enumeration<? extends ZipEntry> entries = zf.entries();
				while(entries.hasMoreElements())
				{
					ZipEntry entry = entries.nextElement();
					if (new File(entry.getName()).getName().equals(theFilename))
					{
						try(BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(entry))))
						{
							mColumns = br.readLine().split("\t", -1).length;
							mRows = 0;
							while ((br.readLine()) != null)
							{
								mRows += 1;
							}
						}
						return;
					}
				}
			}
		}
	}

	public void findAndCopy(Path theZipPath, String theFilename, File theCopyTo) throws IOException
	{
		if (theZipPath.toFile().exists())
		{
			try(ZipFile zf = new ZipFile(theZipPath.toFile()))
			{
				Enumeration<? extends ZipEntry> entries = zf.entries();
				while(entries.hasMoreElements())
				{
					ZipEntry entry = entries.nextElement();
					if (new File(entry.getName()).getName().equals(theFilename))
					{
						try(BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(entry))))
						{
							try (BufferedWriter bw = java.nio.file.Files.newBufferedWriter(Paths.get(theCopyTo.getAbsolutePath()), Charset.availableCharsets().get("UTF-8")))
							{
								String line = br.readLine();
								while (line != null)
								{
									bw.write(line);
									bw.newLine();
									line = br.readLine();
								}
							}
						}
						return;
					}
				}
			}
		}
	}
}
