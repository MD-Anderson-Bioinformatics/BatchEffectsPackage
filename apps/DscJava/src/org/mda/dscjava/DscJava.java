/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.mda.dscjava;

import java.util.ArrayList;
import java.util.concurrent.ForkJoinPool;

/**
 *
 * @author tdcasasent
 */
public class DscJava
{

	public DscJava()
	{
		//
	}

	static public double test()
	{
		double[] data =
		{
			1.41642285, 0.10694894, -1.23408485, -0.66788459, -0.72350214, -0.39967497,
			-1.20617325, 2.02534118, 0.04824682, -0.23245758, 0.58790321, -1.62881050,
			-1.04373111, -0.80310432, 1.20987098, -1.56504974, -0.28991085, -0.18063447,
			0.50197599, 1.38438714, -0.28225272, -0.52935371, -0.79628805, -1.99388328,
			0.61984947, -0.30598282, 0.46818999, -1.07130663, -1.92119735, 0.09369464,
			-1.37164781, -1.54513036, 0.66259059, -0.09702015, -0.81402846, -0.19330046,
			-1.46182162, 1.29318353, 1.18095836, 1.77052505, 1.25097229, 1.51530726,
			0.44529141, -1.11088399, -2.28993315, -1.17506637, -0.25687152, 0.53933316,
			1.67527333, 1.74301675, -0.94686291, -1.46451769, 0.44926929, -0.73454580,
			-0.72579708, -1.33489991, -0.46994309, 0.51241297, 1.12108049, -1.36887916,
			1.66774028, 0.49402844, 0.36732663, 0.10993039, -0.07856253, -0.20510748,
			0.65800235, -0.03943472, -0.57510140, -0.68024807, 2.12295046, -1.03431131,
			0.24757499, 0.09078522, -0.47703802, 0.43342363, 0.29852667, -0.72552589,
			0.67071872, 0.68971619, -1.73853724, 0.96191800, 0.80947295, -1.06464023,
			-0.87065928, -0.54730635, 1.32059294, 0.36257119, 2.14670453, -1.05164755,
			-0.60255080, -0.86457967, 0.86925408, -0.27599802, 1.68193731, 0.20491488,
			-0.02002944, 0.64227494, -0.47946586, -0.15497779
		};
		int[] dim =
		{
			4, 25
		};
		String[] batchIds =
		{
			"a", "b", "c", "a", "b",
			"a", "b", "c", "a", "b",
			"a", "b", "c", "a", "b",
			"a", "b", "c", "a", "b",
			"a", "b", "c", "a", "b"
		};
		PcaDsc pcaDsc = new DscJava().getDSCwithExcerpt(data, dim, batchIds);
		System.out.println("getDSCwithExcerpt");
		System.out.println(pcaDsc.toString());
		PcaDsc [] pcaDscList = new DscJava().doDscPerms(data, dim, batchIds, 200, 5);
		ArrayList<Double> dscVals = new ArrayList<>();
		for(PcaDsc pca : pcaDscList)
		{
			System.out.println("doDscPerms");
			System.out.println(pca.toString());
			dscVals.add(pca.mDSC);
		}
		return pcaDsc.mDSC;
	}

	public void getDscJavaVersion()
	{
		System.out.println("DscJava 2018-07-31-1230");
	}

	protected double[][] convertToMatrix(double[] theValues, int[] theDim)
	{
		double[][] matrix = new double[theDim[0]][theDim[1]];
		int index = 0;
		for (int y = 0; y < theDim[1]; y++)
		{
			for (int x = 0; x < theDim[0]; x++)
			{
				matrix[x][y] = theValues[index];
				index = index + 1;
			}
		}
		return matrix;
	}

	public void matrixTest(double[] theValues, int[] theDim, String[] theBatchIds)
	{
		try
		{
			double[][] matrix = convertToMatrix(theValues, theDim);
			System.out.println("matrix");
			for (double x : matrix[0])
			{
				System.out.print(x + "\t");
			}
			System.out.println();
			//////////////
			System.out.println("theBatchIds");
			for (String x : theBatchIds)
			{
				System.out.print(x + "\t");
			}
			System.out.println();
		}
		catch (Exception exp)
		{
			System.out.println("matrixTest " + exp.getMessage());
			System.err.println("matrixTest " + exp.getMessage());
			exp.printStackTrace(System.out);
			exp.printStackTrace(System.err);
		}
	}

	// called from R
	public PcaDsc getDSCwithExcerpt(double[] theValues, int[] theDim, String[] theBatchIds)
	{
		PcaDsc result = null;
		try
		{
			getDscJavaVersion();
			System.out.println("Number of Components/Genes theDim[0]=" + theDim[0]);
			System.out.println("Number of Samples theDim[1]=" + theDim[1]);
			System.out.println("theBatchIds.length=" + theBatchIds.length);
			double[][] matrix = convertToMatrix(theValues, theDim);
			System.out.println("matrix.length=" + matrix.length);
			System.out.println("matrix[0].length=" + matrix[0].length);
			result = DSC.calc(matrix, theBatchIds);
			//System.out.println(result.toString());
		}
		catch (Exception exp)
		{
			System.out.println("getDSCwithExcerpt " + exp.getMessage());
			System.err.println("getDSCwithExcerpt " + exp.getMessage());
			exp.printStackTrace(System.out);
			exp.printStackTrace(System.err);
		}
		catch (Throwable thr)
		{
			System.out.println("Err: getDSCwithExcerpt " + thr.getMessage());
			System.err.println("Err: getDSCwithExcerpt " + thr.getMessage());
			thr.printStackTrace(System.out);
			thr.printStackTrace(System.err);
		}
		return result;
	}

	// called from R
	public PcaDsc[] doDscPerms(double[] theValues, int[] theDim, String[] theBatchIds, int thePerms, int theThreads)
	{
		PcaDsc[] results = null;
		try
		{
			getDscJavaVersion();
			System.out.println("Number of Components/Genes theDim[0]=" + theDim[0]);
			System.out.println("Number of Samples theDim[1]=" + theDim[1]);
			System.out.println("theBatchIds.length=" + theBatchIds.length);
			double[][] matrix = convertToMatrix(theValues, theDim);
			System.out.println("matrix.length=" + matrix.length);
			System.out.println("matrix[0].length=" + matrix[0].length);
			System.out.println("perms=" + thePerms);
			System.out.println("threads=" + theThreads);
			System.out.flush();
			ForkJoinPool pool = new ForkJoinPool(theThreads);
			ArrayList<PcaDsc> myresults = pool.invoke(new Perm(matrix, theBatchIds, thePerms, theThreads));
			//System.out.println("results contains " + myresults.size());
			results = myresults.toArray(new PcaDsc[0]);
			//System.out.println(results[0].toString());
		}
		catch (Exception exp)
		{
			System.out.println("pvalueDSCwithExcerpt " + exp.getMessage());
			System.err.println("pvalueDSCwithExcerpt " + exp.getMessage());
			exp.printStackTrace(System.out);
			exp.printStackTrace(System.err);
		}
		return results;
	}
}
