package edu.mda.bioinfo.boxplotjava;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.HistogramType;

/**
 *
 * @author tdcasasent
 */
public abstract class Boxplot_Mixin
{
	/// MATRIX
	// mMatrixData[M_BARCODES][M_GENES]

	public double[][] mMatrixData = null;
	public String[] mGeneList = null;
	public String[] mSampleList = null;
	// BATCH TYPES
	public String[] mBatchTypeList = null;
	// mBatchValues[batchTypeIndex][barcodeIndex]
	public String[][] mBatchValues = null;
	public File mOutputDir = null;
	public String mTitle = "";
	// building histograms
	//public TreeMap<Integer, TreeMap<String, ArrayList<Double>>> mHistogramDataByBins = null;
	//static protected int [] mBinCounts = { 10, 20, 30, 40, 50 };
	public TreeMap<String, ArrayList<Double>> mHistogramData = null;
	
	public Boxplot_Mixin(
			double[][] theMatrixData,
			String[] theGeneList,
			String[] theSampleList,
			String[] theBatchTypeList,
			String[][] theBatchValues,
			File theOutputDir,
			String theTitle)
	{
		mOutputDir = theOutputDir;
		mMatrixData = theMatrixData;
		mGeneList = theGeneList;
		mSampleList = theSampleList;
		mBatchTypeList = theBatchTypeList;
		mBatchValues = theBatchValues;
		mTitle = theTitle;
		mHistogramData = new TreeMap<>();
	}

	public double [] getDoubleData(ArrayList<Double> theList)
	{
		double [] result = new double[theList.size()];
		for (int x = 0; x < theList.size(); x++)
		{
			result[x] = theList.get(x);
		}
		return result;
	}

	public int [] getIntData(ArrayList<Integer> theList)
	{
		int [] result = new int[theList.size()];
		for (int x = 0; x < theList.size(); x++)
		{
			result[x] = theList.get(x);
		}
		return result;
	}

	public int [] getNonNAindexes(double [] theValues)
	{
		ArrayList<Integer> nonNA = new ArrayList<>();
		for(int x=0;x<theValues.length;x++)
		{
			double value = theValues[x];
			if ( (false==Double.isInfinite(value)) &&
				 (false==Double.isNaN(value)) )
			{
				nonNA.add(x);
			}
		}
		return getIntData(nonNA);
	}

	public String [] getNonNaStrings(int [] theIndexes, String [] theValues)
	{
		ArrayList<String> values = new ArrayList<>();
		for(int index : theIndexes)
		{
			values.add(theValues[index]);
		}
		return values.toArray(new String[0]);
	}

	public double [] getNonNaDoubles(int [] theIndexes, double [] theValues)
	{
		ArrayList<Double> values = new ArrayList<>();
		for(int index : theIndexes)
		{
			values.add(theValues[index]);
		}
		return getDoubleData(values);
	}

	public void process() throws IOException, Exception
	{
		processInternal();
	}

	protected abstract void processInternal() throws IOException, Exception;
	
	protected void processHistogram(String theFileStub) throws IOException, Exception
	{
		BoxplotJava.log("Boxplot_Mixin::processHistogram - started");
		BoxplotJava.log("Boxplot_Mixin::processHistogram - theFileStub=" + theFileStub);
		String pngFile = theFileStub + ".png";
		String tsvFile = theFileStub + ".tsv";
		TreeMap<String, ArrayList<Double>> histMap = mHistogramData;
		BoxplotJava.log("Boxplot_Mixin::processHistogram - mHistogramData.size()=" + mHistogramData.size());
		// build list of histogram data
		ArrayList<JFreeChart> histList = new ArrayList<>();
		TreeMap<String, ArrayList<Double>> histFileOutput = new TreeMap<>();
		int maxLength = 0;
		for (Entry<String, ArrayList<Double>> entry : histMap.entrySet())
		{
			HistogramDataset hd = new HistogramDataset();
			hd.setType(HistogramType.FREQUENCY);
			double [] vals = entry.getValue().stream().mapToDouble(Double::doubleValue).toArray();
			int binCount = computeBinCount(vals);
			hd.addSeries("freq", vals, binCount);
			histList.add(ChartFactory.createHistogram(entry.getKey(), null, null, hd, PlotOrientation.VERTICAL, false, false, false));
			ArrayList<Double> values = new ArrayList<>();
			for(int x=0;x<binCount;x++)
			{
				values.add(hd.getX(0, x).doubleValue());
				values.add(hd.getY(0, x).doubleValue());
			}
			if (binCount>maxLength)
			{
				maxLength = binCount;
			}
			histFileOutput.put(entry.getKey(), values);
		}
		BoxplotJava.log("Boxplot_Mixin::processHistogram - histFileOutput.size()=" + histFileOutput.size());
		try(BufferedWriter bw = Files.newBufferedWriter(Paths.get(tsvFile), Charset.availableCharsets().get("ISO-8859-1")))
		{
			bw.write("entry");
			bw.write("\t");
			bw.write("size");
			for(int x=0;x<maxLength;x++)
			{
				bw.write("\tx" + x);
				bw.write("\ty" + x);
			}
			bw.newLine();
			//
			for (Entry<String, ArrayList<Double>> entry : histFileOutput.entrySet())
			{
				bw.write(entry.getKey());
				bw.write("\t");
				ArrayList<Double> vals = entry.getValue();
				bw.write(Integer.toString(vals.size()/2));
				for(Double val : vals)
				{
					bw.write("\t");
					bw.write(val.toString());
				}
				bw.newLine();
			}
		}
		// write png of histogram data
		final CombinedDomainXYPlot plot = new CombinedDomainXYPlot();
		for(int x=0;x<histList.size();x++)
		{
			plot.add(histList.get(x).getXYPlot());
		}
		ChartUtilities.saveChartAsPNG(new File(pngFile), new JFreeChart(plot), 1600, (200*histList.size()));
	}
	
	protected void buildHistogramData(String theId, double [] theValues) throws Exception
	{
		ArrayList<Double> myIdValues = mHistogramData.get(theId);
		if (null==myIdValues)
		{
			myIdValues = new ArrayList<>();
		}
		else
		{
			throw new Exception("Duplicate id for histogram data: " + theId);
		}
		for (double val : theValues)
		{
			myIdValues.add(val);
		}
		mHistogramData.put(theId, myIdValues);
	}
	

	protected int computeBinCount(double [] dValues)
	{
		// https://github.com/meyerjp3/psychometrics/blob/master/src/main/java/com/itemanalysis/psychometrics/histogram/FreedmanDiaconisBinCalculation.java
		// http://data.broadinstitute.org/igv/projects/downloads/IGVDistribution_2.1.11/src/org/broad/igv/tools/converters/ExpressionFormatter.java
		int binCount = 0;
		double iqr = StatUtils.percentile(dValues, 75) - StatUtils.percentile(dValues, 25);
		if (iqr == 0) 
		{
			double [] deviations = new double[dValues.length];
			for(int i=0;i<dValues.length;i++)
			{
				deviations[i] = Math.abs(dValues[i]);
			}
			// MAD, as defined at http://stat.ethz.ch/R-manual/R-devel/library/stats/html/mad.html
            iqr = 1.4826 * StatUtils.percentile(deviations, 50);
		}
		if (iqr > 0) 
		{
			// R code: ceiling(diff(range(x))/(2 * h * length(x)^(-1/3)))
			List<Double> dlist = Arrays.asList(ArrayUtils.toObject(dValues));
			double diffrange = Collections.max(dlist) - Collections.min(dlist);
			double pow = Math.pow(dValues.length,(1.0/3.0));
			double h = 2.0*(iqr/pow);
			double bins = diffrange/h;
			binCount = new Double(Math.ceil(bins)).intValue();
		}
		if (binCount<1)
		{
			binCount = 1;
		}
		return binCount;
	}

}
