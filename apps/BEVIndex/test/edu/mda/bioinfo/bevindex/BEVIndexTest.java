/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.mda.bioinfo.bevindex;

import java.io.File;
import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.AfterClass;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 *
 * @author linux
 */
public class BEVIndexTest
{
	
	public BEVIndexTest()
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
		File staticDir = new File("/BatchEffectsPackage_data/testing_static/BEVIndex/archive");
		File dynamicDir = new File("/BatchEffectsPackage_data/testing_dynamic/BEVIndex");
		FileUtils.deleteQuietly(dynamicDir);
		try
		{
			FileUtils.copyDirectoryToDirectory(staticDir, dynamicDir);
		}
		catch(Exception exp)
		{
			exp.printStackTrace(System.err);
			fail(exp.getMessage());
		}
	}
	
	@After
	public void tearDown()
	{
	}

	@Test
	public void testAddExternalIndexToDataRun() throws Exception
	{
		try
		{
//			System.out.println("addExternalIndexToDataRun");
//			String theDataRunDir = "/BatchEffectsPackage_data/testing_dynamic/BEVIndex/archive";
//			String[] theDataRuns = { "2018-07-11-1200"};
//			String[] theDefault = { "PCA" };
//			String[] theNameRuns = { "Test Display Run"};
//			String[] theLabelRuns = { "Test Display Run Label"};
//			String[] theFilesRuns = { "/BatchEffectsPackage_data/testing_dynamic/BEVIndex/archive/GDC_2018-07-11-1200.json" };
//			String[] theLevelLabels = { "Test Run Options", "Algorithm", "Diagram Type", "Sub-Type" };
//			String theTooltipTsv = "/BatchEffectsPackage_data/testing_static/BEVIndex/data/tooltip_gdc.tsv";
//			BEVIndex.addExternalIndexToDataRun(theDataRunDir, theDataRuns, theDefault, theNameRuns, theLabelRuns, theFilesRuns, theLevelLabels, theTooltipTsv);
//			System.out.println(new File(theDataRunDir, "GDC_2018-07-11-1200.json").getAbsolutePath());
//			System.out.println(new File(theDataRunDir, "compare.json").getAbsolutePath());
//			assertTrue("The files differ!", FileUtils.contentEquals(new File(theDataRunDir, "GDC_2018-07-11-1200.json"), new File(theDataRunDir, "compare.json")));
		}
		catch(Exception exp)
		{
			exp.printStackTrace(System.err);
			fail(exp.getMessage());
		}
	}
}
