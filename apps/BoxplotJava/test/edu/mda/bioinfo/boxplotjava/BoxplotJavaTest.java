/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.mda.bioinfo.boxplotjava;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Comparator;
import java.util.TreeSet;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author linux
 */
public class BoxplotJavaTest
{
	
	public BoxplotJavaTest()
	{
	}
	
	@BeforeClass
	public static void setUpClass()
	{
	}
	
	@AfterClass
	public static void tearDownClass()
	{
	}
	
	@Before
	public void setUp()
	{
	}
	
	@After
	public void tearDown()
	{
	}

	protected void addToList(TreeSet<String> theList, Path thePath, Path theParent)
	{
		String a = thePath.toAbsolutePath().toString();
		String b = theParent.toAbsolutePath().toString();
		theList.add(a.replace(b, ""));
	}
	
	protected void compareFiles(String theFile1, String theFile2) throws IOException, Exception
	{
		System.out.println("Checking file1: " + theFile1);
		System.out.println("Checking file2: " + theFile2);
		File file1 = new File(theFile1);
		File file2 = new File(theFile2);
		try(BufferedReader br1 = Files.newBufferedReader(Paths.get(file1.getAbsolutePath()), StandardCharsets.UTF_8))
		{
			try(BufferedReader br2 = Files.newBufferedReader(Paths.get(file2.getAbsolutePath()), StandardCharsets.UTF_8))
			{
				String line1 = br1.readLine();
				String line2 = br2.readLine();
				while(null!=line2)
				{
					if (!line1.equals(line2))
					{
						throw new Exception("Files do not match: " + theFile1);
					}
					line1 = br1.readLine();
					line2 = br2.readLine();
				}
			}
		}
	}

	/**
	 * Test of mainFile method, of class BoxplotJava.
	 */
	@Test
	public void testMainFile()
	{
		System.out.println("mainFile");
		String[] args = { "../../data/testing_static/Boxplot/data", "../../data/testing_dynamic/Boxplot"};
		String testOut = "../../data/testing_static/Boxplot/output";
		try
		{
			if (new File(args[1]).exists())
			{
				Files.walk(Paths.get(args[1]))
					.sorted(Comparator.reverseOrder())
					.map(Path::toFile)
					.forEach(File::delete);
			}
			BoxplotJava.mainFile(args);
			////////////////////////////////////////////////////////
			TreeSet<String> staticOutput = new TreeSet<>();
			Path staticOutPath = Paths.get(testOut);
			Files.walk(staticOutPath)
				.filter(Files::isRegularFile)
				.forEach(p -> addToList(staticOutput, p, staticOutPath));
			////////////////////////////////////////////////////////
			TreeSet<String> testOutput = new TreeSet<>();
			Path testOutPath = Paths.get(args[1]);
			Files.walk(testOutPath)
				.filter(Files::isRegularFile)
				.forEach(p -> addToList(testOutput, p, testOutPath));
			////////////////////////////////////////////////////////
			// 3224 - updated since it was missing histogram files.
			if (testOutput.size()!=3224)
			{
				fail("Expected 3224 files but found " + testOutput.size());
			}
			for(String compareFile : staticOutput)
			{
				compareFiles(new File(args[1], compareFile).getAbsolutePath(), 
						new File(testOut, compareFile).getAbsolutePath());
			}
		}
		catch(Exception exp)
		{
			exp.printStackTrace(System.err);
			fail(exp.getMessage());
		}
	}
	
}
