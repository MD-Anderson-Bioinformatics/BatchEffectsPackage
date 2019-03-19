/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.mda.dscjava;

import java.util.ArrayList;
import org.apache.commons.lang3.ArrayUtils;

/**
 *
 * @author tdcasasent
 */
public class PcaDsc
{
	public ArrayList<Double> mListOfGeneDSC = new ArrayList<>();
	public ArrayList<Double> mListOfGeneDB = new ArrayList<>();
	public ArrayList<Double> mListOfGeneDW = new ArrayList<>();
	public double mDSC = new Double(0);
	public double mDB = new Double(0);
	public double mDW = new Double(0);

	@Override
	public String toString()
	{
		return "PcaDsc{" + "mListOfGeneDSC=" + mListOfGeneDSC + ", mListOfGeneDB=" + mListOfGeneDB + ", mListOfGeneDW=" + mListOfGeneDW + ", mDSC=" + mDSC + ", mDB=" + mDB + ", mDW=" + mDW + '}';
	}

	public double getmDB()
	{
		return DSC.epsilonInfinityCheck(mDB);
	}

	public double getmDSC()
	{
		return DSC.epsilonInfinityCheck(mDSC);
	}

	public double getmDW()
	{
		return DSC.epsilonInfinityCheck(mDW);
	}

	public double [] getmListOfGeneDB()
	{
		return ArrayUtils.toPrimitive(DSC.epsilonInfinityCheck(mListOfGeneDB).toArray(new Double[0]));
	}

	public double [] getmListOfGeneDSC()
	{
		return ArrayUtils.toPrimitive(DSC.epsilonInfinityCheck(mListOfGeneDSC).toArray(new Double[0]));
	}

	public double [] getmListOfGeneDW()
	{
		return ArrayUtils.toPrimitive(DSC.epsilonInfinityCheck(mListOfGeneDW).toArray(new Double[0]));
	}
}
