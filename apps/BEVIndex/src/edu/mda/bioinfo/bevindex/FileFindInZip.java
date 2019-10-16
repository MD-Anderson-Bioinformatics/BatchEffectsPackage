/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.mda.bioinfo.bevindex;

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
