/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.mda.readrjava;

import java.nio.charset.Charset;

/**
 *
 * @author tdcasasent
 */
public class ReadRJava
{
	public static String M_VERSION = "2018-07-30-1430";

	public static Charset getCharset()
	{
		return Charset.availableCharsets().get("ISO-8859-1");
	}

	/**
	 * @param args the command line arguments
	 */
	/*
	public static void main(String[] args)
	{
		try
		{
			// TODO code application logic here
			long start = System.currentTimeMillis();
			JavaFile jf = loadStringData_PadHeaders("/Users/tdcasasent/Desktop/sort/test.tsv");
			for(String column : jf.mColumns)
			{
				System.out.println("column " + column);
			}
			long stop = System.currentTimeMillis();
			System.out.println("Complete " + (stop-start));
		}
		catch(Exception exp)
		{
			exp.printStackTrace(System.out);
		}
	}
*/
	public static JavaFile loadDoubleData(String theFile, boolean theCols, boolean theRows, boolean theData)
	{
		System.out.println("ReadRJavaL::loadDoubleData " + M_VERSION);
		JavaFile result = null;
		try
		{
			result = new JavaFile(theFile, theCols, theRows, theData, false, false);
			result.loadDoubleData();
			if (result.mHeaderCount!=result.mRowColumnCount)
			{
				throw new Exception("result.mHeaderCount (" + result.mHeaderCount + ") does not equal result.mRowCount (" + result.mRowColumnCount + ")");
			}
		}
		catch(Exception exp)
		{
			result = null;
			System.err.println("ERROR loadDoubleData loading file " + theFile + " theCols=" + theCols + " theRows=" + theRows + " theData=" + theData);
			exp.printStackTrace(System.err);
		}
		System.out.println("ReadRJavaL::loadDoubleData done");
		return result;
	}

	public static JavaFile loadStringDataSkipFirstLine(String theFile)
	{
		System.out.println("ReadRJavaL::loadStringDataSkipFirstLine " + M_VERSION);
		JavaFile result = null;
		try
		{
			result = new JavaFile(theFile, true, false, true, true, false);
			result.loadStringData();
			if (result.mHeaderCount!=result.mRowColumnCount)
			{
				throw new Exception("result.mHeaderCount (" + result.mHeaderCount + ") does not equal result.mRowCount (" + result.mRowColumnCount + ")");
			}
		}
		catch(Exception exp)
		{
			result = null;
			System.err.println("ERROR loadStringDataSkipFirstLine loading file " + theFile);
			exp.printStackTrace(System.err);
		}
		System.out.println("ReadRJavaL::loadStringDataSkipFirstLine done");
		return result;
	}

	public static JavaFile loadStringData(String theFile)
	{
		System.out.println("ReadRJavaL::loadStringData " + M_VERSION);
		JavaFile result = null;
		try
		{
			// row reading is determined by file format, so pass false
			result = new JavaFile(theFile, true, false, true, false, false);
			result.loadStringData();
			if (result.mHeaderCount!=result.mRowColumnCount)
			{
				throw new Exception("result.mHeaderCount (" + result.mHeaderCount + ") does not equal result.mRowCount (" + result.mRowColumnCount + ")");
			}
		}
		catch(Exception exp)
		{
			result = null;
			System.err.println("ERROR loadStringData loading file " + theFile );
			exp.printStackTrace(System.err);
		}
		System.out.println("ReadRJavaL::loadStringData done");
		return result;
	}

	public static JavaFile loadStringData_PadHeaders(String theFile)
	{
		System.out.println("ReadRJavaL::loadStringData_PadHeaders " + M_VERSION);
		JavaFile result = null;
		try
		{
			// row reading is determined by file format, so pass false
			result = new JavaFile(theFile, true, false, true, false, true);
			result.loadStringData();
			if (result.mHeaderCount!=result.mRowColumnCount)
			{
				int counter = 1;
				for (int x=result.mHeaderCount; x<result.mRowColumnCount; x++)
				{
					System.out.println("WARN-READ_PAD-HEADER - padding headers for " + theFile);
					result.mColumns.add("unnamed_" +  String.format("%07d", counter));
					counter = counter + 1;
				}
			}
		}
		catch(Exception exp)
		{
			result = null;
			System.err.println("ERROR loadStringData_PadHeaders loading file " + theFile );
			exp.printStackTrace(System.err);
		}
		System.out.println("ReadRJavaL::loadStringData_PadHeaders done");
		return result;
	}

	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////

	public static boolean writeDoubleData_All(String theFile, String [] theCols, String [] theRows, double [] theData)
	{
		System.out.println("ReadRJava::writeDoubleData_All " + M_VERSION);
		boolean result = false;
		try
		{
			result = JavaFile.writeFile(theFile, theCols, theRows, null, theData);
		}
		catch(Exception exp)
		{
			result = false;
			System.err.println("ERROR writeDoubleData_All writing file " + theFile);
			exp.printStackTrace(System.err);
		}
		System.out.println("ReadRJava::writeDoubleData_All done");
		return result;
	}

	public static boolean writeDoubleData_Column(String theFile, String [] theCols, double [] theData)
	{
		System.out.println("ReadRJava::writeDoubleData_Column" + M_VERSION);
		boolean result = false;
		try
		{
			result = JavaFile.writeFile(theFile, theCols, null, null, theData);
		}
		catch(Exception exp)
		{
			result = false;
			System.err.println("ERROR writeDoubleData_Columnwriting file " + theFile);
			exp.printStackTrace(System.err);
		}
		System.out.println("ReadRJava::writeDoubleData_Columndone");
		return result;
	}

	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////

	public static boolean writeStringData_All(String theFile, String [] theCols, String [] theRows, String [] theData)
	{
		System.out.println("ReadRJava::writeStringData_All " + M_VERSION);
		boolean result = false;
		try
		{
			result = JavaFile.writeFile(theFile, theCols, theRows, theData, null);
		}
		catch(Exception exp)
		{
			result = false;
			System.err.println("ERROR writeStringData_All writing file " + theFile);
			exp.printStackTrace(System.err);
		}
		System.out.println("ReadRJava::writeStringData_All done");
		return result;
	}

	public static boolean writeStringData_Column(String theFile, String [] theCols, String [] theData)
	{
		System.out.println("ReadRJava::writeStringData_Column " + M_VERSION);
		boolean result = false;
		try
		{
			result = JavaFile.writeFile(theFile, theCols, null, theData, null);
		}
		catch(Exception exp)
		{
			result = false;
			System.err.println("ERROR writeStringData_Column writing file " + theFile);
			exp.printStackTrace(System.err);
		}
		System.out.println("ReadRJava::writeStringData_Column done");
		return result;
	}

}
