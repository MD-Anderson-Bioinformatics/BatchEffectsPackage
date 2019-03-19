/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.mda.readrjava;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import org.junit.After;
import org.junit.AfterClass;
import static org.junit.Assert.fail;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 *
 * @author linux
 */
public class ReadRJavaTest
{
	
	public ReadRJavaTest()
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

	/**
	 * Test of getCharset method, of class ReadRJava.
	 */
//	@Test
//	public void testGetCharset()
//	{
//		System.out.println("getCharset");
//		Charset expResult = null;
//		Charset result = ReadRJava.getCharset();
//		assertEquals(expResult, result);
//		// TODO review the generated test code and remove the default call to fail.
//		fail("The test case is a prototype.");
//	}

	/**
	 * Test of loadDoubleData method, of class ReadRJava.
	 */
	@Test
	public void testDoubleData()
	{
		try
		{
			System.out.println("loadDoubleData");
			String inputFile = "../../data/testing_static/MATRIX_DATA/brca_agi4502_matrix_data.tsv";
			String outputFile = "../../data/testing_dynamic/ReadRJava/brca_agi4502_matrix_data.tsv";
			new File(outputFile).delete();
			new File(outputFile).getParentFile().mkdirs();
			boolean theCols = true;
			boolean theRows = true;
			boolean theData = true;
			JavaFile result = ReadRJava.loadDoubleData(inputFile, theCols, theRows, theData);
			ReadRJava.writeDoubleData_All(outputFile, result.getmColumns(), result.getmRows(), result.getmDoubleData());
			compareFiles(outputFile, inputFile);
		}
		catch(Exception exp)
		{
			exp.printStackTrace(System.err);
			System.err.flush();
			fail(exp.getMessage());
		}
	}

	@Test
	public void testDataframeData()
	{
		try
		{
			System.out.println("loadDoubleData");
			String inputFile = "../../data/testing_static/MATRIX_DATA/brca_agi4502_batches.tsv";
			String outputFile = "../../data/testing_dynamic/ReadRJava/brca_agi4502_batches.tsv";
			new File(outputFile).delete();
			new File(outputFile).getParentFile().mkdirs();
			JavaFile result = ReadRJava.loadStringData(inputFile);
			ReadRJava.writeStringData_All(outputFile, result.getmColumns(), result.getmRows(), result.getmStringData());
			compareFiles(outputFile, inputFile);
		}
		catch(Exception exp)
		{
			exp.printStackTrace(System.err);
			System.err.flush();
			fail(exp.getMessage());
		}
	}


	public void compareFiles(String theNewFile, String theCompareFile) throws IOException, Exception
	{
		System.out.println("Checking file: " + theNewFile);
		File newFile = new File(theNewFile);
		String compareFile = newFile.getName();
		try(BufferedReader br1 = Files.newBufferedReader(Paths.get(newFile.getAbsolutePath()), StandardCharsets.UTF_8))
		{
			try(BufferedReader br2 = Files.newBufferedReader(Paths.get(theCompareFile), ReadRJava.getCharset()))
			{
				String line1 = br1.readLine();
				String line2 = br2.readLine();
				while(null!=line2)
				{
					if (!line1.equals(line2))
					{
						throw new Exception("Files do not match: " + theNewFile);
					}
					line1 = br1.readLine();
					line2 = br2.readLine();
				}
			}
		}
	}
//	/**
//	 * Test of loadStringDataSkipFirstLine method, of class ReadRJava.
//	 */
//	@Test
//	public void testLoadStringDataSkipFirstLine()
//	{
//		System.out.println("loadStringDataSkipFirstLine");
//		String theFile = "";
//		JavaFile expResult = null;
//		JavaFile result = ReadRJava.loadStringDataSkipFirstLine(theFile);
//		assertEquals(expResult, result);
//		// TODO review the generated test code and remove the default call to fail.
//		fail("The test case is a prototype.");
//	}
//
//	/**
//	 * Test of loadStringData method, of class ReadRJava.
//	 */
//	@Test
//	public void testLoadStringData()
//	{
//		System.out.println("loadStringData");
//		String theFile = "";
//		JavaFile expResult = null;
//		JavaFile result = ReadRJava.loadStringData(theFile);
//		assertEquals(expResult, result);
//		// TODO review the generated test code and remove the default call to fail.
//		fail("The test case is a prototype.");
//	}
//
//	/**
//	 * Test of loadStringData_PadHeaders method, of class ReadRJava.
//	 */
//	@Test
//	public void testLoadStringData_PadHeaders()
//	{
//		System.out.println("loadStringData_PadHeaders");
//		String theFile = "";
//		JavaFile expResult = null;
//		JavaFile result = ReadRJava.loadStringData_PadHeaders(theFile);
//		assertEquals(expResult, result);
//		// TODO review the generated test code and remove the default call to fail.
//		fail("The test case is a prototype.");
//	}
//
//	/**
//	 * Test of writeDoubleData_All method, of class ReadRJava.
//	 */
//	@Test
//	public void testWriteDoubleData_All()
//	{
//		System.out.println("writeDoubleData_All");
//		String theFile = "";
//		String[] theCols = null;
//		String[] theRows = null;
//		double[] theData = null;
//		boolean expResult = false;
//		boolean result = ReadRJava.writeDoubleData_All(theFile, theCols, theRows, theData);
//		assertEquals(expResult, result);
//		// TODO review the generated test code and remove the default call to fail.
//		fail("The test case is a prototype.");
//	}
//
//	/**
//	 * Test of writeDoubleData_Column method, of class ReadRJava.
//	 */
//	@Test
//	public void testWriteDoubleData_Column()
//	{
//		System.out.println("writeDoubleData_Column");
//		String theFile = "";
//		String[] theCols = null;
//		double[] theData = null;
//		boolean expResult = false;
//		boolean result = ReadRJava.writeDoubleData_Column(theFile, theCols, theData);
//		assertEquals(expResult, result);
//		// TODO review the generated test code and remove the default call to fail.
//		fail("The test case is a prototype.");
//	}
//
//	/**
//	 * Test of writeStringData_All method, of class ReadRJava.
//	 */
//	@Test
//	public void testWriteStringData_All()
//	{
//		System.out.println("writeStringData_All");
//		String theFile = "";
//		String[] theCols = null;
//		String[] theRows = null;
//		String[] theData = null;
//		boolean expResult = false;
//		boolean result = ReadRJava.writeStringData_All(theFile, theCols, theRows, theData);
//		assertEquals(expResult, result);
//		// TODO review the generated test code and remove the default call to fail.
//		fail("The test case is a prototype.");
//	}
//
//	/**
//	 * Test of writeStringData_Column method, of class ReadRJava.
//	 */
//	@Test
//	public void testWriteStringData_Column()
//	{
//		System.out.println("writeStringData_Column");
//		String theFile = "";
//		String[] theCols = null;
//		String[] theData = null;
//		boolean expResult = false;
//		boolean result = ReadRJava.writeStringData_Column(theFile, theCols, theData);
//		assertEquals(expResult, result);
//		// TODO review the generated test code and remove the default call to fail.
//		fail("The test case is a prototype.");
//	}
//	
}
