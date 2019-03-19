package edu.mda.bioinfo.boxplotjava;

import static edu.mda.bioinfo.boxplotjava.BoxplotJava.getVersion;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeSet;

/**
 *
 * @author tdcasasent
 */
public class Boxplot_AllSampleData extends Boxplot_Mixin
{
	public static void run(ReadMatrixFile theMatrix, ReadBatchFile theBatch, File theOutputDir, String theTitle)
		throws IOException, Exception
	{
		Boxplot_AllSampleData bp;
		{
			BoxplotJava.log(getVersion());
			double [][] matrixData = theMatrix.mCombinedData;
			BoxplotJava.log("Boxplot_AllSampleData::run - matrixData.length=" + matrixData.length);
			BoxplotJava.log("Boxplot_AllSampleData::run - matrixData[0].length=" + matrixData[0].length);
			String [] geneList = theMatrix.geneEqInIndexOrder();
			BoxplotJava.log("Boxplot_AllSampleData::run - geneList.length=" + geneList.length);
			String [] sampleList = theMatrix.mSampleList.toArray(new String[0]);
			BoxplotJava.log("Boxplot_AllSampleData::run - sampleList.length=" + sampleList.length);
			String [] batchTypes = theBatch.mBatchTypes;
			BoxplotJava.log("Boxplot_AllSampleData::run - batchTypes.length=" + batchTypes.length);
			String [][] batches = theBatch.mBatches;
			BoxplotJava.log("Boxplot_AllSampleData::run - batches.length=" + batches.length);
			BoxplotJava.log("Boxplot_AllSampleData::run - batches[0].length=" + batches[0].length);
			bp = new Boxplot_AllSampleData(matrixData, geneList, sampleList, batchTypes, batches,
					new File(theOutputDir, "AllSample-Data"), theTitle);
		}
		bp.process();
	}

	public Boxplot_AllSampleData(
			double[][] theMatrixData,
			String[] theGeneList,
			String[] theSampleList,
			String[] theBatchTypeList,
			String[][] theBatchValues,
			File theOutputDir,
			String theTitle)
	{
		super(theMatrixData, theGeneList, theSampleList, theBatchTypeList, theBatchValues, theOutputDir, theTitle);
	}

	@Override
	protected void processInternal() throws IOException, Exception
	{
		BoxplotJava.log("Boxplot_AllSampleData::processInternal - start");
		for (int batchTypeIndex=0;batchTypeIndex<mBatchTypeList.length;batchTypeIndex++)
		{
			mHistogramData.clear();
			TreeSet<BoxplotElement> elementWrappers = new TreeSet<>();
			String batchTypeName = mBatchTypeList[batchTypeIndex];
			BoxplotJava.log("Boxplot_AllSampleData::processInternal - batchTypeName=" + batchTypeName);
			ArrayList<String> groups = new ArrayList<>();
			for(int sampleIndex=0;sampleIndex<mSampleList.length;sampleIndex++)
			{
				int [] nonNaIndexes = getNonNAindexes(mMatrixData[sampleIndex]);
				//
				double [] values = getNonNaDoubles(nonNaIndexes, mMatrixData[sampleIndex]);
				String [] valLabels = getNonNaStrings(nonNaIndexes, mGeneList);
				String elementLabel = mSampleList[sampleIndex];
				String groupLabel = mBatchValues[batchTypeIndex][sampleIndex];
				int groupId = groups.indexOf(groupLabel);
				if (groupId < 0)
				{
					groups.add(groupLabel);
					groupId = groups.indexOf(groupLabel);
				}
				BoxplotJava.log("================================");
				BoxplotJava.log("Boxplot_AllSampleData");
				BoxplotJava.log("================================");
				BoxplotJava.log("batchTypeIndex=" + batchTypeIndex);
				BoxplotJava.log("values.length=" + values.length);
				BoxplotJava.log("valLabels.length=" + valLabels.length);
				BoxplotJava.log("elementLabel=" + elementLabel);
				BoxplotJava.log("groupLabel=" + groupLabel);
				BoxplotJava.log("groupId=" + groupId);
				BoxplotJava.log("================================");
				elementWrappers.add(BoxplotImpl.getBoxplotData(values, valLabels, elementLabel, groupLabel, groupId,
						mMatrixData[sampleIndex].length));
				// Add building histogram data here
				// boxplot id: elementLabel
				// values: value
				buildHistogramData(elementLabel, values);
			}
			BoxplotJava.log("Boxplot_AllSampleData::processInternal - after inside loop");
			// processInternal elementWrappers, write file and write image
			BoxplotImpl.writeBoxplotFiles(elementWrappers,
					new File(mOutputDir, "BoxPlot_AllSample-Data_Diagram-" + batchTypeName + ".png"),
					new File(mOutputDir, "BoxPlot_AllSample-Data_Annotations-" + batchTypeName + ".tsv"),
					new File(mOutputDir, "BoxPlot_AllSample-Data_BoxData-" + batchTypeName + ".tsv"),
					new File(mOutputDir, "BoxPlot_AllSample-Data_CatData-" + batchTypeName + "-").getAbsolutePath(),
					"feature values", batchTypeName, mTitle, batchTypeName);
			processHistogram(new File(mOutputDir, "BoxPlot_AllSample-Data_Histogram-" + batchTypeName).getAbsolutePath());
			BoxplotJava.log("Boxplot_AllSampleData::processInternal - after files");
		}
		BoxplotJava.log("Boxplot_AllSampleData::processInternal - finished");
	}
}
