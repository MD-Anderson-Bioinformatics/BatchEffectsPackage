/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.mda.bioinfo.bevmanifest;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.TreeSet;
import static java.util.stream.Collectors.toList;
import java.util.stream.Stream;

/**
 *
 * @author linux
 */
public class FileFind
{
	public TreeSet<Path> mMatches = null;

	FileFind()
	{
		mMatches = new TreeSet<>();
	}

	void find(Path thePath, String theFilename) throws IOException
	{
		if (thePath.toFile().isDirectory())
		{
			//System.out.println("dir: " + thePath);
			// check files in this directory
			boolean match = false;
			try (Stream<Path> myFiles = Files.list(thePath))
			{
				for (Path file : myFiles.filter(Files::isRegularFile).collect(toList()))
				{
					//System.out.println("file: " + file);
					if (theFilename.equals(file.getFileName().toString()))
					{
						//System.out.println("matched file: " + file);
						mMatches.add(file);
						match = true;
						// return exist for loop, not function
						return;
					}
				}
			}
			if (false==match)
			{
				// if we reach this spot, no file was found, so continue down
				// check direcotories in this directory
				//System.out.println("no match found. Check sub-dirs of: " + thePath);
				try (Stream<Path> myFiles = Files.list(thePath))
				{
					for (Path dir : myFiles.filter(Files::isDirectory).collect(toList()))
					{
						find(dir, theFilename);
					}
				}
			}
		}
		// should only reach here if nothing found above
	}
}
