/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.mda.bioinfo.bevindex;

import static edu.mda.bioinfo.bevindex.BEVIndex.mVersion;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 *
 * @author linux
 */
public class CheckNumbers
{

	////////////////////////////////////////////////////////////////////////////
	//// main
	////////////////////////////////////////////////////////////////////////////
	public static void main(String[] args)
	{
		System.out.println(mVersion);
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
