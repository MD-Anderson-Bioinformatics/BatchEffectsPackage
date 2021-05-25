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
	public double mDSC = 0;
	public double mDB = 0;
	public double mDW = 0;

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
