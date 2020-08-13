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

package org.mda.dscjava;

import java.util.ArrayList;
import java.util.concurrent.RecursiveTask;

/**
 *
 * @author tdcasasent
 */
public class Perm extends RecursiveTask<ArrayList<PcaDsc>>
{
	private double [][] mMatrix = null;
	private String [] mBatchIds = null;
	private int mThreads = 0;
	private int mPerms = 0;
	private int mCurrentPerm = 0;

	public synchronized int getAndIncrementCurrentPerm()
	{
		mCurrentPerm = mCurrentPerm + 1;
		return(mCurrentPerm - 1);
	}

	public final synchronized double [][] getMatrix()
	{
		return(mMatrix);
	}

	public final int getPerms()
	{
		return(mPerms);
	}

	public final String [] getBatchIds()
	{
		return(mBatchIds);
	}

	public Perm(double [][] theMatrix, String [] theBatchIds, int thePerms, int theThreads)
	{
		super();
		mMatrix = theMatrix;
		mBatchIds = theBatchIds;
		mPerms = thePerms;
		mCurrentPerm = 0;
		mThreads = theThreads;
	}

	@Override
	protected ArrayList<PcaDsc> compute()
	{
		ArrayList<PcaDsc> results = new ArrayList<>();
		//System.out.println("from perm " + mPerms);
		System.out.println("Perm::compute - mPerms=" + mPerms);
		System.out.println("Perm::compute - mThreads=" + mThreads);
		System.out.flush();
		try
		{
			ArrayList<DSC> permList = new ArrayList<>();
			for (int threads=0;threads<mThreads;threads++)
			{
				DSC doCalc = new DSC(this);
				permList.add(doCalc);
				doCalc.fork();
				System.out.println("Perm::compute - fork");
				System.out.flush();
			}
			System.out.println("Perm::compute - do joins");
			System.out.flush();
			for(DSC myDSC : permList)
			{
				ArrayList<PcaDsc> dscList = myDSC.join();
				results.addAll(dscList);
			}
			System.out.println("after joins");
			System.out.flush();
		}
		catch(Exception exp)
		{
			System.out.println("Error in compute " + exp.getMessage());
			System.err.println("Error in compute " + exp.getMessage());
			exp.printStackTrace(System.out);
			exp.printStackTrace(System.err);
		}
		return results;
	}

}
