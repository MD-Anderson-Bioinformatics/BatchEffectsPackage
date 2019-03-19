package edu.mda.bioinfo.boxplotjava;

import static edu.mda.bioinfo.boxplotjava.BoxplotJava.getVersion;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;

/**
 *
 * @author tdcasasent
 */
public class Boxplot_GroupFunction extends Boxplot_Mixin
{
	public GroupFunctions mGroupFunction = GroupFunctions.MEAN;

	public static void run(ReadMatrixFile theMatrix, ReadBatchFile theBatch,
		File theOutputDir, GroupFunctions theFunction, String theTitle) throws IOException, Exception
	{
		Boxplot_GroupFunction bp;
		{
			BoxplotJava.log(getVersion());
			double [][] matrixData = theMatrix.mCombinedData;
			String [] geneList = theMatrix.geneEqInIndexOrder();
			String [] sampleList = theMatrix.mSampleList.toArray(new String[0]);
			String [] batchTypes = theBatch.mBatchTypes;
			String [][] batches = theBatch.mBatches;
			bp = new Boxplot_GroupFunction(matrixData, geneList, sampleList, batchTypes, batches,
					new File(theOutputDir, "Group-" + theFunction.name()), theFunction, theTitle);
		}
		bp.process();
	}

	public Boxplot_GroupFunction(
			double[][] theMatrixData,
			String[] theGeneList,
			String[] theSampleList,
			String[] theBatchTypeList,
			String[][] theBatchValues,
			File theOutputDir,
			GroupFunctions theFunction,
			String theTitle)
	{
		super(theMatrixData, theGeneList, theSampleList, theBatchTypeList, theBatchValues, theOutputDir, theTitle);
		mGroupFunction = theFunction;
	}

	@Override
	protected void processInternal() throws IOException, Exception
	{
		BoxplotJava.log("Boxplot_GroupFunction::processInternal - start");
		for (int batchTypeIndex=0;batchTypeIndex<mBatchTypeList.length;batchTypeIndex++)
		{
			mHistogramData.clear();
			TreeSet<BoxplotElement> elementWrappers = new TreeSet<>();
			String batchTypeName = mBatchTypeList[batchTypeIndex];
			BoxplotJava.log("Boxplot_GroupFunction::processInternal - batchTypeName=" + batchTypeName);
			// list of batches
			BoxplotJava.log("Boxplot_GroupFunction::processInternal - list of batches");
			String [] batchList;
			{
				TreeSet<String> batchSet = new TreeSet<>();
				batchSet.addAll(Arrays.asList(mBatchValues[batchTypeIndex]));
				batchList = batchSet.toArray(new String[0]);
			}
			// group function values for samples
			BoxplotJava.log("Boxplot_GroupFunction::processInternal - group function values for samples");
			double [] sampleGroupValues = new double[mSampleList.length];
			int maxPossible = mSampleList.length;
			for (int x=0;x<mSampleList.length;x++)
			{
				sampleGroupValues[x] = GroupFunctions.doFunction(
						getNonNaDoubles(getNonNAindexes(mMatrixData[x]), mMatrixData[x]), mGroupFunction);
			}
			// iterate through batches
			BoxplotJava.log("Boxplot_GroupFunction::processInternal - iterate through batches");
			for(int batchIndex=0;batchIndex<batchList.length;batchIndex++)
			{
				String batchId = batchList[batchIndex];
				ArrayList<String> batchLabelList = new ArrayList<>();
				ArrayList<Double> batchValueList = new ArrayList<>();
				BoxplotJava.log("Boxplot_GroupFunction::processInternal - iterate through batches - before inside loop");
				for(int valIndex=0;valIndex<sampleGroupValues.length;valIndex++)
				{
					if (batchId.equals(mBatchValues[batchTypeIndex][valIndex]))
					{
						batchLabelList.add(mSampleList[valIndex]);
						batchValueList.add(sampleGroupValues[valIndex]);
					}
				}
				BoxplotJava.log("Boxplot_GroupFunction::processInternal - iterate through batches - after inside loop");
				//
				double [] values = getDoubleData(batchValueList);
				String [] valLabels = batchLabelList.toArray(new String[0]);
				String elementLabel = batchId;
				String groupLabel = batchId;
				int groupId = batchIndex;
				BoxplotJava.log("================================");
				BoxplotJava.log("Boxplot_GroupFunction");
				BoxplotJava.log("================================");
				BoxplotJava.log("batchTypeIndex=" + batchTypeIndex);
				BoxplotJava.log("values.length=" + values.length);
				BoxplotJava.log("valLabels.length=" + valLabels.length);
				BoxplotJava.log("elementLabel=" + elementLabel);
				BoxplotJava.log("groupLabel=" + groupLabel);
				BoxplotJava.log("groupId=" + groupId);
				BoxplotJava.log("================================");
				elementWrappers.add(BoxplotImpl.getBoxplotData(values, valLabels, elementLabel,
						groupLabel, groupId, maxPossible));
				// Add building histogram data here
				// boxplot id: elementLabel
				// values: value
				buildHistogramData(elementLabel, values);
			}
			BoxplotJava.log("Boxplot_GroupFunction::processInternal - after iterate through batches");
			// processInternal elementWrappers, write file and write image
			BoxplotJava.log("Boxplot_GroupFunction::processInternal - call writeBoxplotFiles");
			BoxplotImpl.writeBoxplotFiles(elementWrappers,
					new File(mOutputDir,
					"BoxPlot_Group-" + mGroupFunction.name() + "_Diagram-" + batchTypeName + ".png"),
					new File(mOutputDir, "BoxPlot_Group-" + mGroupFunction.name() + "_Annotations-" + batchTypeName + ".tsv"),
					new File(mOutputDir, "BoxPlot_Group-" + mGroupFunction.name() + "_BoxData-" + batchTypeName + ".tsv"),
					new File(mOutputDir, "BoxPlot_Group-" + mGroupFunction.name() + "_CatData-" + batchTypeName + "-").getAbsolutePath(),
					"feature values", batchTypeName, mTitle, batchTypeName);
			BoxplotJava.log("Boxplot_GroupFunction::processInternal - call processHistogram");
			processHistogram(new File(mOutputDir, "BoxPlot_Group-" + mGroupFunction.name() + "_Histogram-" + batchTypeName).getAbsolutePath());
			BoxplotJava.log("Boxplot_GroupFunction::processInternal - after files");
		}
		BoxplotJava.log("Boxplot_GroupFunction::processInternal - finished");
	}
}
