/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.mda.dscjava;

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
public class DscJavaTest
{
	
	public DscJavaTest()
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
	 * Test of main method, of class DscJava.
	 */
	@Test
	public void test()
	{
		double dsc = 0;
		try
		{
			System.out.println("generic test");
			dsc = DscJava.test();
			System.out.println(dsc);
		}
		catch(Exception exp)
		{
			exp.printStackTrace(System.err);
			System.err.flush();
			fail(exp.getMessage());
		}
		assertEquals(0.45348045392686015, dsc, 0.000001);
	}
	
}
