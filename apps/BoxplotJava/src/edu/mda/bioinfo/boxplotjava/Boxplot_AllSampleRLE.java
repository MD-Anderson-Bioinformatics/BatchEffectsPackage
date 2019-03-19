package edu.mda.bioinfo.boxplotjava;

import static edu.mda.bioinfo.boxplotjava.BoxplotJava.getVersion;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeSet;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 *
 * @author tdcasasent
 */
public class Boxplot_AllSampleRLE extends Boxplot_Mixin
{
	public static void run(ReadMatrixFile theMatrix, ReadBatchFile theBatch, File theOutputDir, String theTitle) 
		throws IOException, Exception
	{
		Boxplot_AllSampleRLE bp = null;
		{
			BoxplotJava.log(getVersion());
			double [][] matrixData = theMatrix.mCombinedData;
			String [] geneList = theMatrix.geneEqInIndexOrder();
			String [] sampleList = theMatrix.mSampleList.toArray(new String[0]);
			String [] batchTypes = theBatch.mBatchTypes;
			String [][] batches = theBatch.mBatches;
			bp = new Boxplot_AllSampleRLE(matrixData, geneList, sampleList, batchTypes, batches,
					new File(theOutputDir, "AllSample-RLE"), theTitle);
		}
		bp.process();
	}

	public Boxplot_AllSampleRLE(
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

	public double [] doColumnMedianSubtraction(double [] theColumnValues)
	{
		DescriptiveStatistics ds = new DescriptiveStatistics(theColumnValues);
		double median = ds.getPercentile(50);
		double [] values = new double[theColumnValues.length];
		for(int x=0;x<theColumnValues.length;x++)
		{
			values[x] = theColumnValues[x] - median;
		}
		return values;
	}

	@Override
	protected void processInternal() throws IOException, Exception
	{
		BoxplotJava.log("Boxplot_AllSampleRLE::processInternal - start");
		for (int batchTypeIndex=0;batchTypeIndex<mBatchTypeList.length;batchTypeIndex++)
		{
			mHistogramData.clear();
			TreeSet<BoxplotElement> elementWrappers = new TreeSet<>();
			String batchTypeName = mBatchTypeList[batchTypeIndex];
			BoxplotJava.log("Boxplot_AllSampleRLE::processInternal - batchTypeName=" + batchTypeName);
			ArrayList<String> groups = new ArrayList<>();
			for(int sampleIndex=0;sampleIndex<mSampleList.length;sampleIndex++)
			{
				int [] nonNaIndexes = getNonNAindexes(mMatrixData[sampleIndex]);
				//
				double [] values = doColumnMedianSubtraction(getNonNaDoubles(nonNaIndexes, mMatrixData[sampleIndex]));
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
				BoxplotJava.log("Boxplot_AllSampleRLE");
				BoxplotJava.log("================================");
				BoxplotJava.log("batchTypeIndex=" + batchTypeIndex);
				BoxplotJava.log("values.length=" + values.length);
				BoxplotJava.log("valLabels.length=" + valLabels.length);
				BoxplotJava.log("elementLabel=" + elementLabel);
				BoxplotJava.log("groupLabel=" + groupLabel);
				BoxplotJava.log("groupId=" + groupId);
				BoxplotJava.log("================================");
				elementWrappers.add(BoxplotImpl.getBoxplotData(values, valLabels, elementLabel,
						groupLabel, groupId, mMatrixData[sampleIndex].length));
				// Add building histogram data here
				// boxplot id: elementLabel
				// values: value
				buildHistogramData(elementLabel, values);
			}
			BoxplotJava.log("Boxplot_AllSampleRLE::processInternal - after inside loop");
			// processInternal elementWrappers, write file and write image
			BoxplotImpl.writeBoxplotFiles(elementWrappers,
					new File(mOutputDir, "BoxPlot_AllSample-RLE_Diagram-" + batchTypeName + ".png"),
					new File(mOutputDir, "BoxPlot_AllSample-RLE_Annotations-" + batchTypeName + ".tsv"),
					new File(mOutputDir, "BoxPlot_AllSample-RLE_BoxData-" + batchTypeName + ".tsv"),
					new File(mOutputDir, "BoxPlot_AllSample-RLE_CatData-" + batchTypeName + "-").getAbsolutePath(),
					"feature values", batchTypeName, mTitle, batchTypeName);
			processHistogram(new File(mOutputDir, "BoxPlot_AllSample-RLE_Histogram-" + batchTypeName).getAbsolutePath());
			BoxplotJava.log("Boxplot_AllSampleRLE::processInternal - after files");
		}
		BoxplotJava.log("Boxplot_AllSampleRLE::processInternal - finished");
	}
}
