package edu.mda.bioinfo.boxplotjava;

import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.TreeSet;
import org.apache.commons.lang3.ArrayUtils;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.LegendItemCollection;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.CategoryLabelPositions;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.renderer.category.BoxAndWhiskerRenderer;
import org.jfree.data.statistics.BoxAndWhiskerCalculator;
import org.jfree.data.statistics.BoxAndWhiskerItem;
import org.jfree.data.statistics.DefaultBoxAndWhiskerCategoryDataset;

import org.mda.legendjava.LegendJava;
/**
 *
 * @author tdcasasent
 * @author jmmelott
 */
public class BoxplotImpl
{
	static private double MAX_NUMBER_OF_OUTLIERS = 20000;

	BoxplotImpl()
	{
	}

	static public BoxplotElement getBoxplotData(
			double[] theValues,
			String[] theValLabels,
			String theElementLabel,
			String theGroupLabel,
			int theGroupId,
			int theMaxPossibleValues)
			throws Exception
	{
		//BoxplotJava.log("getBoxplotData start");
		// Statistics to compute
		BoxplotElement element;
		BoxAndWhiskerItem bawItem;

		//BoxplotJava.log("getBoxplotData pre stats");
		bawItem = BoxAndWhiskerCalculator.calculateBoxAndWhiskerStatistics(Arrays.asList(ArrayUtils.toObject(theValues)), true);
		//BoxplotJava.log("getBoxplotData post stats");
		List<Double> outliers = bawItem.getOutliers();
		//BoxplotJava.log("getBoxplotData after getOutliers");

		/*
		if (outliers.size() > MAX_NUMBER_OF_OUTLIERS)
		{
			String outliersErrorMessage = "Excessive outliers. Boxplot processing for batch stopped. "
				+ "Data contains more than " + MAX_NUMBER_OF_OUTLIERS + " outliers. "
				+ outliers.size() + " in " + theValues.length + " values.";
			BoxplotJava.log(outliersErrorMessage);
			throw new Exception(outliersErrorMessage);
		}
		*/

		//BoxplotJava.log("getBoxplotData descriptive");
		//DescriptiveStatistics ds = new DescriptiveStatistics(theValues);
		//double iqr = ds.getPercentile(75) - ds.getPercentile(25);

		double minRegularValue = bawItem.getMinRegularValue().doubleValue();
		double maxRegularValue = bawItem.getMaxRegularValue().doubleValue();

		ArrayList<String> lowerOutLabels = new ArrayList<>();
		ArrayList<Double> lowerOutValues = new ArrayList<>();
		ArrayList<String> upperOutLabels = new ArrayList<>();
		ArrayList<Double> upperOutValues = new ArrayList<>();

		//BoxplotJava.log("getBoxplotData values loop");
		//BoxplotJava.log("getBoxplotData theValues.length=" + theValues.length);
		//BoxplotJava.log("getBoxplotData outliers.size()=" + outliers.size());
		TreeSet<Double> setOutliers = new TreeSet<>();
		setOutliers.addAll(outliers);
		//BoxplotJava.log("getBoxplotData setOutliers.size()=" + setOutliers.size());
		for (int i = 0; i < theValues.length; i++)
		{
			double value = theValues[i];
			//BoxplotJava.log("getBoxplotData contains");
			if (setOutliers.contains(value))
			{
				//BoxplotJava.log("getBoxplotData found");
				if (value < minRegularValue)
				{
					lowerOutValues.add(theValues[i]);
					lowerOutLabels.add(theValLabels[i]);
				}
				if (value > maxRegularValue)
				{
					upperOutValues.add(theValues[i]);
					upperOutLabels.add(theValLabels[i]);
				}
			}
			//BoxplotJava.log("getBoxplotData bottom of for");
		}

		//BoxplotJava.log("getBoxplotData upper lower calc");
		int origLowerSize = lowerOutValues.size();
		int origUpperSize = upperOutValues.size();
		int lowerSize = lowerOutValues.size();
		int upperSize = upperOutValues.size();
		double increment = 1;
		String message = "";

		if (outliers.size() > MAX_NUMBER_OF_OUTLIERS)
		{
			float fractionDisplayed = (float) (MAX_NUMBER_OF_OUTLIERS / outliers.size());
			increment = outliers.size() / MAX_NUMBER_OF_OUTLIERS;
			lowerSize = Math.round(lowerSize * fractionDisplayed);
			upperSize = Math.round(upperSize * fractionDisplayed);
			if (lowerSize > 0)
			{
				message += lowerSize + " of " + origLowerSize + " lower outliers shown.";
			}
			if (upperSize > 0)
			{
				message += upperSize + " of " + origUpperSize + " upper outliers shown.";
			}
		}
		//BoxplotJava.log("getBoxplotData lowerSize=" + lowerSize);
		//BoxplotJava.log("getBoxplotData upperSize=" + upperSize);

		String lowerOutLabelsString[] = new String[lowerSize];
		double lowerOutValuesDouble[] = new double[lowerSize];
		String upperOutLabelsString[] = new String[upperSize];
		double upperOutValuesDouble[] = new double[upperSize];

		//BoxplotJava.log("getBoxplotData lowerOutLabels");
		String calcItem;
		int itemToGet;
		if (lowerOutLabels.size() > 0)
		{
			for (int i = 0; i < lowerSize; i++)
			{
				calcItem = (new Double(increment * i)).toString();
				itemToGet = Math.min(origLowerSize, Integer.parseInt(calcItem.substring(0, calcItem.indexOf("."))));
				lowerOutLabelsString[i] = lowerOutLabels.get(itemToGet);
				lowerOutValuesDouble[i] = lowerOutValues.get(itemToGet);
			}
		}

		//BoxplotJava.log("getBoxplotData upperOutLabels");
		if (upperOutLabels.size() > 0)
		{
			for (int i = 0; i < upperSize; i++)
			{
				calcItem = (new Double(increment * i)).toString();
				itemToGet = Math.min(origUpperSize, Integer.parseInt(calcItem.substring(0, calcItem.indexOf("."))));
				upperOutLabelsString[i] = upperOutLabels.get(itemToGet);
				upperOutValuesDouble[i] = upperOutValues.get(itemToGet);
			}
		}

		//BoxplotJava.log("getBoxplotData BoxplotElement");
		element = new BoxplotElement(bawItem,
				lowerOutLabelsString, lowerOutValuesDouble,
				upperOutLabelsString, upperOutValuesDouble,
				theElementLabel, theGroupLabel,
				theMaxPossibleValues, theValues.length,
				theValues, theValLabels);
		element.setMessage(message);
		//BoxplotJava.log("getBoxplotData finished");

		return element;
	}

	static public void writeBoxplotFiles(TreeSet<BoxplotElement> theElements, File thePNGFile,
			File theAnnotationsFile, File theBoxFile, String theCatPathPrefix,
			String theDataLabel, String theSampleLabel, String theTitle, String theGroup) throws IOException
	{
		BoxplotJava.log("writeBoxplotFiles - boxplotPng");
		boxplotPng(theElements, thePNGFile, theDataLabel, theSampleLabel, theTitle, theGroup);
		BoxplotJava.log("writeBoxplotFiles - boxplotAnnotations");
		boxplotAnnotations(theElements, theAnnotationsFile);
		BoxplotJava.log("writeBoxplotFiles - boxplotBox");
		boxplotBox(theElements, theBoxFile);
		BoxplotJava.log("writeBoxplotFiles - boxplotCat");
		boxplotCat(theElements, theCatPathPrefix);
		BoxplotJava.log("writeBoxplotFiles - done");
	}

	static public void boxplotCat(TreeSet<BoxplotElement> theElements, String theCatPathPrefix) throws IOException
	{
		for (BoxplotElement be : theElements)
		{
			try (BufferedWriter bw = Files.newBufferedWriter(
							Paths.get(theCatPathPrefix + be.getElementLabel() + ".tsv"),
							Charset.availableCharsets().get("ISO-8859-1")))
			{
				bw.write("id\tvalue");
				bw.newLine();
				// all values
				//for(int x=0;x<be.mValues.length;x++)
				//{
				//	String label = be.mValueLabels[x];
				//	double value = be.mValues[x];
				//	bw.write(label + "\t" + Double.toString(value));
				//	bw.newLine();
				//}
				// This does just outliers
				for(int x=0;x<be.getUpperOutLabel().length;x++)
				{
					String label = be.getUpperOutLabel()[x];
					double value = be.getUpperOutValues()[x];
					bw.write(label + "\t" + Double.toString(value));
					bw.newLine();
				}
				for(int x=0;x<be.getLowerOutLabels().length;x++)
				{
					String label = be.getLowerOutLabels()[x];
					double value = be.getLowerOutValues()[x];
					bw.write(label + "\t" + Double.toString(value));
					bw.newLine();
				}
			}
		}
	}

	static public void boxplotBox(TreeSet<BoxplotElement> theElements, File theBoxFile) throws IOException
	{
		try (BufferedWriter bw = Files.newBufferedWriter(
						Paths.get(theBoxFile.getAbsolutePath()),
						Charset.availableCharsets().get("ISO-8859-1")))
		{
			// headers
			bw.write("Id\tLowerOutMax\tLowerOutMin\tLowerNotch\tLowerWhisker\tLowerHinge\tMedian\tUpperHinge\tUpperWhisker\tUpperNotch\tUpperOutMin\tUpperOutMax");
			bw.newLine();
			// non-na data points for each element
			for (BoxplotElement be : theElements)
			{
				//double boxLength = Math.abs(be.getElement().getQ3().doubleValue() - be.getElement().getQ1().doubleValue());
				double upperWisker = be.getElement().getMaxRegularValue().doubleValue();
				double lowerWisker = be.getElement().getMinRegularValue().doubleValue();
				String lowerMinOut = "NA";
				String lowerMaxOut = "NA";
				if (be.getLowerOutValues().length>0)
				{
                                    List<Double> dblList = Arrays.asList(ArrayUtils.toObject(be.getLowerOutValues()));
                                    lowerMinOut =  Double.toString(Collections.min(dblList));
                                    dblList = Arrays.asList(ArrayUtils.toObject(be.getLowerOutValues()));
                                    lowerMaxOut =  Double.toString(Collections.max(dblList));
				}
				String upperMaxOut = "NA";
				String upperMinOut = "NA";
				if (be.getUpperOutValues().length>0)
				{
                                    List<Double> dblList = Arrays.asList(ArrayUtils.toObject(be.getUpperOutValues()));
                                    upperMinOut =  Double.toString(Collections.min(dblList));
                                    dblList = Arrays.asList(ArrayUtils.toObject(be.getUpperOutValues()));
                                    upperMaxOut =  Double.toString(Collections.max(dblList));
				}
				bw.write(be.getElementLabel()
						+ "\t" + lowerMinOut // lower min out
						+ "\t" + lowerMaxOut // lower max out
						+ "\t" + (be.getElement().getMedian().doubleValue() + (-1.58 * (be.getIQR() / Math.sqrt(be.getNonNADataPoints())))) // need to calculate lower notch
						+ "\t" + lowerWisker // lower whisker
						+ "\t" + be.getElement().getQ1().doubleValue() // lower hinge
						+ "\t" + be.getElement().getMedian().doubleValue() // median
						+ "\t" + be.getElement().getQ3().doubleValue() // upper hinge
						+ "\t" + upperWisker // upperer whisker
						+ "\t" + (be.getElement().getMedian().doubleValue() + (1.58 * (be.getIQR() / Math.sqrt(be.getNonNADataPoints())))) // need to calculate upper notch
						+ "\t" + upperMinOut // upper min out
						+ "\t" + upperMaxOut); // upper max out
				bw.newLine();
			}
		}
	}

	static public void boxplotAnnotations(TreeSet<BoxplotElement> theElements, File theAnnotationsFile) throws IOException
	{
		try (BufferedWriter bw = Files.newBufferedWriter(
						Paths.get(theAnnotationsFile.getAbsolutePath()),
						Charset.availableCharsets().get("ISO-8859-1")))
		{
			// headers
			bw.write("key\tvalue");
			bw.newLine();
			// total data points
			bw.write("Total-Data-Points\t" + theElements.first().getTotalPossibleDataPoints());
			bw.newLine();
			// non-na data points for each element
			for (BoxplotElement be : theElements)
			{
				bw.write("Non-NA-Points-" + be.getElementLabel() + "\t" + be.getNonNADataPoints());
				bw.newLine();
			}
		}
	}

	static public void boxplotPng(TreeSet<BoxplotElement> theElements, File thePNGFile,
			String theDataLabel, String theSampleLabel, String theTitle, String theGroup) throws IOException
	{
		BoxplotJava.log("BoxplotImpl::boxplotPng theElements.size()=" + theElements.size());
		BoxplotJava.log("BoxplotImpl::boxplotPng thePNGFile=" + thePNGFile);
		BoxplotJava.log("BoxplotImpl::boxplotPng theDataLabel=" + theDataLabel);
		BoxplotJava.log("BoxplotImpl::boxplotPng theSampleLabel=" + theSampleLabel);
		BoxplotJava.log("BoxplotImpl::boxplotPng theTitle=" + theTitle);
		BoxplotJava.log("BoxplotImpl::boxplotPng theGroup=" + theGroup);
		DefaultBoxAndWhiskerCategoryDataset dataset = new DefaultBoxAndWhiskerCategoryDataset();

		String categoryLabelText;
		for (BoxplotElement be : theElements)
		{
			categoryLabelText = be.getElementLabel() + "   " + be.getGroupLabel();
			dataset.add(be.getElement(), be.getGroupLabel(), categoryLabelText);
		}

		CategoryAxis xAxis = new CategoryAxis();
		xAxis.setCategoryLabelPositions(CategoryLabelPositions.DOWN_90);
		xAxis.setTickLabelFont(new Font("Monospace", Font.PLAIN, 9));
		xAxis.setCategoryMargin(0.0);
		xAxis.setUpperMargin(0.0);
		xAxis.setLowerMargin(0.0);
		NumberAxis yAxis = new NumberAxis(theDataLabel);
		yAxis.setAutoRange(true);

		BoxAndWhiskerRenderer renderer = new BoxAndWhiskerRenderer();
		renderer.setMeanVisible(false);
		renderer.setBaseOutlinePaint(Color.BLACK);

		/*
		JMM
		This is a completely idiotic hack that experimentation came up with that actually worked that probably shouldn't have.
		We are assigning a category to each boxandwhisker item to get to the labels to show up with the sample id.
		Using a positive number or 0 for the margin gave us useless lines instead of boxes.
		For some reason, using a negative number made it work better.
		Since the size of the margins are based on the number of samples I originally used the negative of the sample size
		Then we determined that we wanted the coloring of the boxplots based on the groups rather than the samples,
		so I changed the margin item to the negative of the number of groups rather than samples.
		At this point, the boxes were overlapping each other. Additional experimentation led to a margin values of
		.73 x the group size rather than 1 x the group size.
		*/
		ArrayList<String> groups = new ArrayList<>();
		HashMap<String, ArrayList<String>> groupMessages = new HashMap<>();
		ArrayList<String> tempMessageList;

		for (BoxplotElement be : theElements)
		{
			if (!groups.contains(be.getGroupLabel()))
			{
				groups.add(be.getGroupLabel());
				//System.out.println("Adding group " + be.getGroupLabel() + " and message:" + be.getMessage());
			}
			if (be.getMessage()!=null && be.getMessage().length() > 0)
			{
				tempMessageList = groupMessages.get(be.getGroupLabel());
				if (tempMessageList == null)
				{
					tempMessageList = new ArrayList<>();
					groupMessages.put(be.getGroupLabel(), tempMessageList);
				}
				tempMessageList.add("   " + be.getElementLabel() + " : " + be.getMessage());
			}
		}
		renderer.setItemMargin(-0.73 * groups.size());
		//End completely idiotic hack.

		CategoryPlot plot = new CategoryPlot(dataset, xAxis, yAxis, renderer);

		//Create legend
		LegendItemCollection legendItems = plot.getLegendItems();
		int legendItemCount = legendItems.getItemCount();

		ArrayList<String> legendNameList = new ArrayList<>();
		ArrayList<String> legendColorList = new ArrayList<>();
		String[] legendNames;
		String[] legendColors;

		Color color;
		String legendColor;
		String legendName;
		for (int i=0; i<legendItemCount; i++)
		{
			color = (Color) legendItems.get(i).getFillPaint();
			legendColor = String.format("#%02x%02x%02x", color.getRed(), color.getGreen(), color.getBlue());
			legendColorList.add(legendColor);
			//System.out.println( "color:" + legendColors[i] );

			legendName = legendItems.get(i).getLabel();
			legendNameList.add(legendName);
			tempMessageList = groupMessages.get(legendName);
			if (tempMessageList != null)
			{
				for (String msg : tempMessageList )
				{
					legendColorList.add(legendColor);
					legendNameList.add(msg);
					//System.out.println("Adding message " + msg + " to group " + legendName );
				}
			}
		}
		legendNames = legendNameList.toArray(new String[legendNameList.size()]);
		legendColors = legendColorList.toArray(new String[legendColorList.size()]);

		JFreeChart chart = new JFreeChart(theTitle + " - " + theSampleLabel, new Font("SansSerif", Font.BOLD, 14), plot, false);

		thePNGFile.getParentFile().mkdirs();
		ChartUtilities.saveChartAsPNG(thePNGFile, chart, 1600, 800);

		String legendPNGFileName = thePNGFile.getCanonicalPath().replace("_Diagram-", "_Legend-");
		//System.out.println("legendPNGFileName:" + legendPNGFileName);

		// Write the legend for the chart. Must be done here as the colors for the chart are dynamically created by JFreeChart
		LegendJava.writeLegend(theTitle + " - " + theSampleLabel, "", legendNames, legendColors, null, legendPNGFileName);
	}
}
