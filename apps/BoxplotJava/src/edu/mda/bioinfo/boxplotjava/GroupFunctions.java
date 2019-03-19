package edu.mda.bioinfo.boxplotjava;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 *
 * @author tdcasasent
 */
public enum GroupFunctions
{
	KURTOSIS("Kurtosis", 0),
	MAX("Max", 1),
	MEAN("Mean", 2),
	MIN("Min", 3),
	POPULATION_VARIANCE("PopulationVariance", 4),
	SKEWNESS("Skewness", 5),
	SUM("Sum", 6),
	SUMSQ("SumOfSquares", 7),
	SAMPLE_VARIANCE("SampleVariance", 8);

	private final String mName;
	private final int mId;

	GroupFunctions(String theName, int theId)
	{
		mName = theName;
		mId = theId;
	}

	public static GroupFunctions getGroupFunctionFromId(int theId)
	{
		GroupFunctions result = null;
		for(GroupFunctions gf : GroupFunctions.values())
		{
			if (theId==gf.mId)
			{
				result = gf;
			}
		}
		return result;
	}

	public static double doFunction(double [] theValues, GroupFunctions theFunction)
	{
		DescriptiveStatistics ds = new DescriptiveStatistics(theValues);
		double value = 0;
		switch(theFunction)
		{
			case KURTOSIS:
				value = ds.getKurtosis();
				break;
			case MAX:
				value = ds.getMax();
				break;
			case MEAN:
				value = ds.getMean();
				break;
			case MIN:
				value = ds.getMin();
				break;
			case POPULATION_VARIANCE:
				value = ds.getPopulationVariance();
				break;
			case SKEWNESS:
				value = ds.getSkewness();
				break;
			case SUM:
				value = ds.getSum();
				break;
			case SUMSQ:
				value = ds.getSumsq();
				break;
			case SAMPLE_VARIANCE:
				value = ds.getVariance();
				break;
		}
		return value;
	}
}
