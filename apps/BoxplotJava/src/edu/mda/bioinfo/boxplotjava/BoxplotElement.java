package edu.mda.bioinfo.boxplotjava;

import org.jfree.data.statistics.BoxAndWhiskerItem;

/**
 * @author jmmelott
 */
public class BoxplotElement implements Comparable<BoxplotElement>
{
	private BoxAndWhiskerItem element;
	private String[] lowerOutLabels;
	private double[] lowerOutValues;
	private String[] upperOutLabels;
	private double[] upperOutValues;
	//public String[] mValueLabels;
	//public double[] mValues;
	private String elementLabel;
	private String groupLabel;
	private int mTotalPossibleDataPoints = -1;
	private int mNonNADataPoints = -1;
	private double mIQR = 0;
	private String message = "";

	public BoxplotElement(BoxAndWhiskerItem theElement, String[] theLowerOutLabels, double[] theLowerOutValues,
		String[] theUpperOutLabels, double[] theUpperOutValues, String theElementLabel, String theGroupLabel,
		int theTotalPossibleDataPoints, int theNonNADataPoints, double [] theValues, String [] theValueLabels)
	{
		this.element = theElement;
		this.lowerOutLabels = theLowerOutLabels;
		this.lowerOutValues = theLowerOutValues;
		this.upperOutLabels = theUpperOutLabels;
		this.upperOutValues = theUpperOutValues;
		this.elementLabel = theElementLabel;
		this.groupLabel = theGroupLabel;
		this.mTotalPossibleDataPoints = theTotalPossibleDataPoints;
		this.mNonNADataPoints = theNonNADataPoints;
		this.mIQR = theElement.getQ3().doubleValue() - theElement.getQ1().doubleValue();
		//this.mValueLabels = theValueLabels;
		//this.mValues = theValues;
	}

	public double getIQR()
	{
		return mIQR;
	}

	public int getTotalPossibleDataPoints()
	{
		return mTotalPossibleDataPoints;
	}

	public int getNonNADataPoints()
	{
		return mNonNADataPoints;
	}

	public BoxAndWhiskerItem getElement()
	{
		return element;
	}

	public void setElement(BoxAndWhiskerItem element)
	{
		this.element = element;
	}

	public String[] getLowerOutLabels()
	{
		return lowerOutLabels;
	}

	public void setLowerOutLabel(String[] theLowerOutLabels)
	{
		this.lowerOutLabels = theLowerOutLabels;
	}

	public double[] getLowerOutValues()
	{
		return lowerOutValues;
	}

	public void setLowerOutValue(double[] theLowerOutValue)
	{
		this.lowerOutValues = theLowerOutValue;
	}

	public String[] getUpperOutLabel()
	{
		return upperOutLabels;
	}

	public void setUpperOutLabel(String[] theUpperOutLabel)
	{
		this.upperOutLabels = theUpperOutLabel;
	}

	public double[] getUpperOutValues()
	{
		return upperOutValues;
	}

	public void setUpperOutValue(double[] theUpperOutValues)
	{
		this.upperOutValues = theUpperOutValues;
	}

	public String getElementLabel()
	{
		return elementLabel;
	}

	public void setElementLabel(String elementLabel)
	{
		this.elementLabel = elementLabel;
	}

	public String getGroupLabel()
	{
		return groupLabel;
	}

	public void setGroupLabel(String groupLabel)
	{
		this.groupLabel = groupLabel;
	}

	public String getMessage()
	{
		return message;
	}

	public void setMessage(String message)
	{
		this.message = message;
	}

	@Override
	public int compareTo(BoxplotElement o)
	{
		int compared = 0;
		if ( (this.groupLabel!=null) &&
			 (o.groupLabel!=null))
		{
			compared = this.groupLabel.compareTo(o.groupLabel);
		}
		if (0==compared)
		{
			compared = this.elementLabel.compareTo(o.elementLabel);
		}
		return compared;
	}
}
