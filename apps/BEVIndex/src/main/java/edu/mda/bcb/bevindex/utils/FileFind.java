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

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.TreeSet;
import static java.util.stream.Collectors.toList;
import java.util.stream.Stream;

/**
 *
 * @author Tod-Casasent
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
