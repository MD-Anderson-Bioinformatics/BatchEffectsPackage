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

import edu.mda.bcb.bevindex.BEVIndex;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 *
 * @author Tod-Casasent
 */
public class CheckNumbers
{

	////////////////////////////////////////////////////////////////////////////
	//// main
	////////////////////////////////////////////////////////////////////////////
	public static void main(String[] args)
	{
		System.out.println(BEVIndex.M_VERSION);
		System.out.println("Start CheckNumbers");
		try
		{
			////////////////////////////////////////////////////////////////////
			String baseDir = args[0];
			String zip = args[1];
			String file1 = args[2];
			String file2 = args[3];
			////////////////////////////////////////////////////////////////////
			FileFind ff = new FileFind();
			ff.find(Paths.get(baseDir), zip);
			System.out.println("found " + ff.mMatches.size() + " zip files");
			for (Path path : ff.mMatches)
			{
				FileFindInZip ffiz = new FileFindInZip();
				ffiz.find(path, file1);
				int matrixSamples = ffiz.mColumns-1;
				ffiz = new FileFindInZip();
				ffiz.find(path, file2);
				int batchSamples = ffiz.mRows;
				if (Math.abs(matrixSamples-batchSamples)>5)
				{
					System.out.println(path + " samples were " + matrixSamples + " for matrix and " + batchSamples + " for batch");
				}
			}
		}
		catch (Exception exp)
		{
			exp.printStackTrace(System.err);
		}
		System.out.println("Done CheckNumbers");
	}
}
