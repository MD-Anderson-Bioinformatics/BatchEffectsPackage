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

package edu.mda.bioinfo.bevindex;

import java.io.File;
import org.apache.commons.io.FileUtils;
import static org.junit.Assert.fail;
import org.junit.Before;
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
