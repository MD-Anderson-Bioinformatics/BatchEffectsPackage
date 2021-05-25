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

package edu.mda.bcb.dscjava;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;
import java.util.concurrent.RecursiveTask;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.stat.StatUtils;

/**
 *
 * @author tdcasasent
 */
public class DSC extends RecursiveTask<ArrayList<PcaDsc>>
{
	public Perm mPerm = null;

	public DSC(Perm thePerm)
	{
		super();
		mPerm = thePerm;
	}

	@Override
	protected ArrayList<PcaDsc> compute()
	{
		ArrayList<PcaDsc> results = new ArrayList<PcaDsc>();
		int currentPerm = mPerm.getAndIncrementCurrentPerm();
		int maxPerm = mPerm.getPerms();
		final double [][] matrix = mPerm.getMatrix();
		final String [] batchIds = mPerm.getBatchIds();
		while(currentPerm<maxPerm)
		{
			if (0 == (currentPerm % 100))
			{
				System.out.println("Thread " + Thread.currentThread().getId() + " is starting perm # " + currentPerm);
				System.out.flush();
			}
			//System.out.println("do calc");
			// TODO needs calculations here
			// first permute the data by row
			RandomDataGenerator rand = new RandomDataGenerator();
			double [][] permData = matrix.clone();
			for(int x=0;x<permData.length;x++)
			{
				double [] rowArray = permData[x];
				ArrayList<Double> rowList = new ArrayList<>();
				for(int c=0;c<rowArray.length;c++)
				{
					rowList.add(rowArray[c]);
				}
				Double[] sampleClass = Arrays.asList(rand.nextSample(rowList, rowArray.length)).toArray(new Double[0]);
				permData[x] = ArrayUtils.toPrimitive(sampleClass);
			}
			PcaDsc myObj = calc(permData, batchIds);
			results.add(myObj);
			currentPerm = mPerm.getAndIncrementCurrentPerm();
		}
		return results;
	}

	public static PcaDsc calc(double [][] permData, String [] batchIds)
	{
		// permData now holds the permutated data
		PcaDsc myObj = new PcaDsc();
		// number of genes/components
		int numberOfGenes = permData.length;
		int sampleCount = batchIds.length;
		//message("mean of across all the samples in the dataset")
		double [] meanListAcrossSamples = new double[numberOfGenes];
		//System.out.println("Mean across samples");
		for(int x=0;x<numberOfGenes;x++)
		{
			meanListAcrossSamples[x] = StatUtils.mean(removeNaN(permData[x]));
			//System.out.print(meanListAcrossSamples[x] + ", ");
		}
		//
		double dwValue = 0;
		double dbValue = 0;
		// this removes duplicates in batch ids (it also sorts)
		TreeSet<String> batchIdList = new TreeSet<>();
		batchIdList.addAll(Arrays.asList(batchIds));
		//
		ArrayList<Double> gDSC = new ArrayList<>();
		ArrayList<Double> gDW = new ArrayList<>();
		ArrayList<Double> gDB = new ArrayList<>();
		// for each gene
		for(int x=0;x<numberOfGenes;x++)
		{
			//System.out.println("loop 2");
			double geneDw = 0;
			double geneDb = 0;
			for (String batchName : batchIdList)
			{
				//System.out.println("loop 3");
				double [] compRow = permData[x];
				int numSamplesInBatch = countOccurances(batchIds, batchName);
				double [] valuesForBatch = removeNaN(getBatchOccurances(compRow, batchIds, batchName));
				double batchMean = StatUtils.mean(valuesForBatch);
				double sampleVariance = StatUtils.variance(valuesForBatch);
				geneDw = geneDw + ( ((double)(numSamplesInBatch-1)) / ((double)(sampleCount)) * sampleVariance );
				geneDb = geneDb + ( (((double)(numSamplesInBatch)) / (double)(sampleCount)) *
									( batchMean - meanListAcrossSamples[x] ) *
									( batchMean - meanListAcrossSamples[x] ) );
				//geneDw = geneDw + ((numSamplesInBatch-1)/numberOfSamples*sampleVariance);
				//geneDb = geneDb + ((numSamplesInBatch / numberOfSamples) * ( batchMean - meanListAcrossSamples[x] ) * ( batchMean - meanListAcrossSamples[x] ) );
				//System.out.println("b=" + batchName + " nsb=" + numSamplesInBatch + " ns=" + numberOfSamples + " bm=" + batchMean + " sv=" + sampleVariance + " geneDw=" + geneDw + " geneDb=" + geneDb + " m=" + meanListAcrossSamples[x]);
			}
			dwValue = dwValue + geneDw;
			dbValue = dbValue + geneDb;
			// take square root
			geneDw = Math.sqrt(geneDw);
			geneDb = Math.sqrt(geneDb);
			//
			double geneDSC = Double.NaN;
			//if ((geneDw!=BigDecimal.NaN)&&(1==geneDw.compareTo(BigDecimal.ZERO)))
			geneDw = epsilonZeroCheck(geneDw);
			if (0==geneDw)
			{
				geneDSC = Double.POSITIVE_INFINITY;
			}
			else if (0==epsilonZeroCheck(geneDb))
			{
				geneDSC = 0;
			}
			else if (geneDw!=0)
			{
				geneDSC = geneDb/geneDw;
			}
			gDSC.add(geneDSC);
			gDW.add(geneDw);
			gDB.add(geneDb);
		}
		// take square root
		dwValue = Math.sqrt(dwValue);
		dbValue = Math.sqrt(dbValue);
		// avoid divide by zero error
		// message("avoid divide by zero error")
		// message("dwValue ", dwValue)
		double DSC = Double.NaN;
		//if ((dwValue!=Double.NaN)&&(dwValue > 0))
		if (dwValue > 0)
		{
			DSC = dbValue/dwValue;
		}
		if (batchIds.length<2)
		{
			gDSC = new ArrayList<>();
			gDB = new ArrayList<>();
			for(int x=0;x<batchIds.length;x++)
			{
				gDSC.add(0.0); //nrow(permData)
				gDB.add(0.0); //nrow(permData))
			}
			DSC = 0;
			dbValue = 0;
		}
		myObj.mListOfGeneDSC =  epsilonZeroCheck(gDSC);
		myObj.mListOfGeneDW =  gDW;
		myObj.mListOfGeneDB =  gDB;
		myObj.mDSC = epsilonZeroCheck(DSC);
		myObj.mDB = dbValue;
		myObj.mDW = dwValue;
		//System.out.println("RESULTS myObj.mDSC =" + myObj.mDSC);
		//System.out.println("RESULTS myObj.mDB =" + myObj.mDB);
		//System.out.println("RESULTS myObj.mDW =" + myObj.mDW);
		return myObj;
	}

/*
	private PcaDsc computeBD()
	{
		PcaDsc myObj = new PcaDsc();
		// number of samples
		int numberOfSamples = permData[0].length;
		//message("mean of across all the samples in the dataset")
		double [] meanListAcrossSamples = new double[numberOfSamples];
		//System.out.println("Mean across samples");
		for(int x=0;x<numberOfSamples;x++)
		{
			meanListAcrossSamples[x] = StatUtils.mean(removeNaN(permData[x]));
			//System.out.print(meanListAcrossSamples[x] + ", ");
		}
		//
		BigDecimal dwValue = BigDecimal.ZERO;
		BigDecimal dbValue = BigDecimal.ZERO;
		// this removes duplicates in batch ids (it also sorts)
		TreeSet<String> batchIdList = new TreeSet<>();
		batchIdList.addAll(Arrays.asList(batchIds));
		//
		ArrayList<Double> gDSC = new ArrayList<>();
		ArrayList<Double> gDW = new ArrayList<>();
		ArrayList<Double> gDB = new ArrayList<>();
		// for each gene
		for(int x=0;x<batchIds.length;x++)
		{
			//System.out.println("loop 2");
			BigDecimal geneDw = BigDecimal.ZERO;
			BigDecimal geneDb = BigDecimal.ZERO;
			for (String batchName : batchIdList)
			{
				//System.out.println("loop 3");
				double [] compRow = permData[x];
				int numSamplesInBatch = countOccurances(batchIds, batchName);
				double [] valuesForBatch = removeNaN(getBatchOccurances(compRow, batchIds, batchName));
				double batchMean = StatUtils.mean(valuesForBatch);
				double sampleVariance = StatUtils.variance(valuesForBatch);
				geneDw = geneDw.add(((new BigDecimal(numSamplesInBatch-1)).divide(
						new BigDecimal(numberOfSamples), 60, RoundingMode.HALF_UP).multiply(
						new BigDecimal(sampleVariance))));
				geneDb = geneDb.add( (new BigDecimal(numSamplesInBatch).divide(
						new BigDecimal(numberOfSamples), 60, RoundingMode.HALF_UP)).multiply(
						( new BigDecimal(batchMean).subtract(new BigDecimal(meanListAcrossSamples[x] )))).multiply(
						( new BigDecimal(batchMean).subtract(new BigDecimal(meanListAcrossSamples[x] )))));
				//geneDw = geneDw + ((numSamplesInBatch-1)/numberOfSamples*sampleVariance);
				//geneDb = geneDb + ((numSamplesInBatch / numberOfSamples) * ( batchMean - meanListAcrossSamples[x] ) * ( batchMean - meanListAcrossSamples[x] ) );
				//System.out.println("b=" + batchName + " nsb=" + numSamplesInBatch + " ns=" + numberOfSamples + " bm=" + batchMean + " sv=" + sampleVariance + " geneDw=" + geneDw + " geneDb=" + geneDb + " m=" + meanListAcrossSamples[x]);
			}
			dwValue = dwValue.add(geneDw);
			dbValue = dbValue.add(geneDb);
			// take square root
			geneDw = new BigDecimal(Math.sqrt(geneDw.doubleValue()));
			geneDb = new BigDecimal(Math.sqrt(geneDb.doubleValue()));
			//
			double geneDSC = Double.NaN;
			//if ((geneDw!=BigDecimal.NaN)&&(1==geneDw.compareTo(BigDecimal.ZERO)))
			if (1==geneDw.compareTo(BigDecimal.ZERO))
			{
				geneDSC = geneDb.divide(geneDw, 60, RoundingMode.HALF_UP).doubleValue();
			}
			gDSC.add(geneDSC);
			gDW.add(geneDw.doubleValue());
			gDB.add(geneDb.doubleValue());
		}
		// take square root
		dwValue = new BigDecimal(Math.sqrt(dwValue.doubleValue()));
		dbValue = new BigDecimal(Math.sqrt(dbValue.doubleValue()));
		// avoid divide by zero error
		// message("avoid divide by zero error")
		// message("dwValue ", dwValue)
		double DSC = Double.NaN;
		//if ((dwValue!=Double.NaN)&&(dwValue > 0))
		if (1==dwValue.compareTo(BigDecimal.ZERO))
		{
			DSC = dbValue.divide(dwValue, 60, RoundingMode.HALF_UP).doubleValue();
		}
		if (batchIds.length<2)
		{
			gDSC = new ArrayList<>();
			gDB = new ArrayList<>();
			for(int x=0;x<batchIds.length;x++)
			{
				gDSC.add(0.0); //nrow(permData)
				gDB.add(0.0); //nrow(permData))
			}
			DSC = 0;
			dbValue = BigDecimal.ZERO;
		}
		myObj.mListOfGeneDSC =  epsilonZeroCheck(gDSC);
		myObj.mListOfGeneDW =  gDW;
		myObj.mListOfGeneDB =  gDB;
		myObj.mDSC = epsilonZeroCheck(DSC);
		myObj.mDB = dbValue.doubleValue();
		myObj.mDW = dwValue.doubleValue();
		//System.out.println("RESULTS myObj.mDSC =" + myObj.mDSC);
		//System.out.println("RESULTS myObj.mDB =" + myObj.mDB);
		//System.out.println("RESULTS myObj.mDW =" + myObj.mDW);
		return myObj;
	}
*/
	static public double epsilonZeroCheck(double theValue)
	{
		// if theValue is less than or equal 1x10^-7, return 0, otherwise return theValue
		if (theValue==Double.NaN)
		{
			return(theValue);
		}
		else if ((Math.abs(theValue)<=0.0000001)&&(Math.abs(theValue)>0))
		{
			return(0);
		}
		else
		{
			return(theValue);
		}
	}

	static public ArrayList<Double> epsilonZeroCheck(ArrayList<Double> theValues)
	{
		for(int x=0; x<theValues.size();x++)
		{
			theValues.set(x, epsilonZeroCheck(theValues.get(x)));
		}
		return theValues;
	}

	static public double epsilonInfinityCheck(double theValue)
	{
		if (theValue==Double.NaN)
		{
			return(theValue);
		}
		else if (theValue==Double.NEGATIVE_INFINITY)
		{
			return(theValue);
		}
		else if (theValue==Double.POSITIVE_INFINITY)
		{
			return(theValue);
		}
		else if (theValue>=1000000)
		{
			return(Double.POSITIVE_INFINITY);
		}
		else if (theValue<=-1000000)
		{
			return(Double.NEGATIVE_INFINITY);
		}
		else if ((theValue<=0.0000001)&&(theValue>0))
		{
			return(0);
		}
		else
		{
			return(theValue);
		}
	}

	static public ArrayList<Double> epsilonInfinityCheck(ArrayList<Double> theValues)
	{
		for(int x=0; x<theValues.size();x++)
		{
			theValues.set(x, epsilonInfinityCheck(theValues.get(x)));
		}
		return theValues;
	}

	static public int countOccurances(String [] theArray, String theString)
	{
		int count = 0;
		for (String myElement : theArray)
		{
			if (theString.equals(myElement))
			{
				count = count + 1;
			}
		}
		return count;
	}

	static public double [] getBatchOccurances(double [] theCompRow, String [] theBatchIds, String theBatchid)
	{
		ArrayList<Double> results = new ArrayList<>();
		for (int x = 0; x<theCompRow.length; x++)
		{
			if (theBatchid.equals(theBatchIds[x]))
			{
				results.add(theCompRow[x]);
			}
		}
		return ArrayUtils.toPrimitive(results.toArray(new Double[0]));
	}

	static public double [] removeNaN(double [] theOriginal)
	{
		ArrayList<Double> results = new ArrayList<>();
		for (double val : theOriginal)
		{
			if (val!= Double.NaN)
			{
				results.add(val);
			}
		}
		return ArrayUtils.toPrimitive(results.toArray(new Double[0]));
	}
}
