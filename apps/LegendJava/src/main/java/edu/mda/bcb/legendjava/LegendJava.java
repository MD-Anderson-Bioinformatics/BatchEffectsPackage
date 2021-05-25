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

package edu.mda.bcb.legendjava;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import javax.imageio.ImageIO;

/**
 *
 * @author tdcasasent
 */
public class LegendJava
{

	public static void main(String[] args)
	{
		String path = args[0];
		try
		{
			String title = "my batch type id or other very long title my batch type id or other very long title my batch type id or other very long title";
			String version = "Version 00.11.22";
			String[] legendNames =
			{
				"00 A1 - UCSF (13)",
				"01 GM - MD Anderson (15)",
				"02 D8 - Greater Poland Cancer Center (63)",
				"03 HN - Ontario Institute for Cancer Research (OICR) (1)",
				"04 Lorem ipsum dolor sit amet, consectetur adipiscing elit",
				"05 Curabitur non eros eget nibh egestas mollis",
				"06 Ut ipsum sem, rutrum eget aliquam luctus, posuere at metus",
				"07 Nam nisi dui, pulvinar vitae dignissim ullamcorper, dignissim id metus",
				"08 Lorem ipsum dolor sit amet, consectetur adipiscing elit",
				"09 Mauris posuere arcu et leo fringilla rhoncus",
				"10 Donec nec neque tellus, sed accumsan magna",
				"11 In massa sem, pretium vitae aliquet vehicula, varius id risus",
				"12 Aenean nec cursus enim",
				"13 Suspendisse congue elementum facilisis",
				"14 Nulla non faucibus libero"
			};
			String[] legendColors =
			{
				"#FF0000", //0
				"#00FF00", //1
				"#0000FF", //2
				"#FA0000", //3
				"#AD0000", //4
				"#00FA00", //5
				"#00AF00", //6
				"#0000FA", //7
				"#0000AF", //8
				"#FAFAFA", //9
				"#AFAFAF", //10
				"#FF00FF", //11
				"#00FFFF", //12
				"#FFFF00", //13
				"#FAFFAF"  //14
			};

			int[] legendSymbols =
			{
				0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14
			};
			writeLegend(title, version, legendNames, legendColors, legendSymbols, new File(path,"names_colors_symbols.png").getAbsolutePath());
			writeLegend(title, version, legendNames, legendColors, null, new File(path,"names_colors.png").getAbsolutePath());
			writeLegend(title, version, legendNames, null, null, new File(path,"names.png").getAbsolutePath());
			writeLegend(title, version, legendNames, null, legendSymbols, new File(path,"names_symbols.png").getAbsolutePath());
			String [] listOfFiles =
			{
				new File(path,"names_colors_symbols.png").getAbsolutePath(),
				new File(path,"names_colors.png").getAbsolutePath(),
				new File(path,"names.png").getAbsolutePath(),
				new File(path,"names_symbols.png").getAbsolutePath()
			};
			combineLegends(title + title + title, listOfFiles, new File(path,"combined.png").getAbsolutePath());
		}
		catch (Exception exp)
		{
			exp.printStackTrace(System.out);
		}
	}
	public BufferedImage mBi = null;
	public Graphics2D mG2d = null;

	protected LegendJava(BufferedImage theBi, Graphics2D theG2d)
	{
		mBi = theBi;
		mG2d = theG2d;
	}
	public static Font mFONT = new Font("Arial", Font.PLAIN, 20);

	public static LegendJava getGraphics(int theWidth, int theHeight)
	{
		BufferedImage bi = new BufferedImage(theWidth, theHeight, BufferedImage.TYPE_INT_ARGB);
		Graphics2D ig2 = bi.createGraphics();
		ig2.setFont(mFONT);
		RenderingHints renderHints = new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		renderHints.put(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BICUBIC);
		renderHints.put(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		ig2.setRenderingHints(renderHints);
		return new LegendJava(bi, ig2);
	}

	public static int getTextHeight()
	{
		FontMetrics fontMetrics = getGraphics(1600, 50).mG2d.getFontMetrics();
		int stringHeight = fontMetrics.getAscent();
		return stringHeight;
	}

	public static int getTextWidth(String theText)
	{
		FontMetrics fontMetrics = getGraphics(1600, 50).mG2d.getFontMetrics();
		int stringWidth = fontMetrics.stringWidth(theText);
		return stringWidth;
	}

	public static ArrayList<String> wrapText(String theText, int theWidth)
	{
		ArrayList<String> results = new ArrayList<>();
		int lineHeight = getTextHeight();
		String[] arr = theText.split(" ");
		int nIndex = 0;
		while (nIndex < arr.length)
		{
			String line = arr[nIndex++];
			while ((nIndex < arr.length) && (getTextWidth(line + " " + arr[nIndex]) < theWidth))
			{
				line = line + " " + arr[nIndex];
				nIndex++;
			}
			results.add(line);
		}
		return results;
	}

	public static String getVersion()
	{
		return "LegendJava 2013_05_03_0823";
	}

	public static boolean combineLegends(String theTitle, String [] theListOfFiles, String theFilenamePath)
	{
		boolean success = false;
		try
		{
			System.out.println(getVersion());
			System.out.println("combineLegends theTitle = " + theTitle);
			System.out.println("combineLegends theFilenamePath = " + theFilenamePath);

			ArrayList<BufferedImage> loadedImages = new ArrayList<>();
			for(String file : theListOfFiles)
			{
				loadedImages.add(ImageIO.read(new File(file)));
			}
			// determine largest legend sizes
			int maxWidth = 0;
			int maxHeight = 0;
			for(BufferedImage image : loadedImages)
			{
				if (image.getWidth()>maxWidth)
				{
					maxWidth = image.getWidth();
				}
				if (image.getHeight()>maxHeight)
				{
					maxHeight = image.getHeight();
				}
			}
			// determine layout
			int numRows = 1;
			int numCols = 1;
			if (theListOfFiles.length<=3)
			{
				numRows = 1;
				numCols = theListOfFiles.length;
			}
			else // two rows
			{
				numRows = 2;
				numCols = (int)Math.ceil(theListOfFiles.length/2.0);
			}
			// calculate size
			int paddingLine = 5;
			int paddingImages = 25;
			int imageWidth = numCols * (paddingImages + maxWidth + paddingImages);
			ArrayList<String> titleNames = LegendJava.wrapText(theTitle, imageWidth - paddingLine - paddingLine);
			int lineHeight = LegendJava.getTextHeight();
			int titleHeight = lineHeight + paddingLine + ((lineHeight + paddingLine)*titleNames.size());
			int imageHeight  = titleHeight + (numRows * (paddingImages + maxHeight + paddingImages));
			//// image
			LegendJava gp = LegendJava.getGraphics(imageWidth, imageHeight);
			// write title
			gp.mG2d.setPaint(Color.black);
			int lineHeightLocation = lineHeight;
			for (String name : titleNames)
			{
				gp.mG2d.drawString(name, (imageWidth - LegendJava.getTextWidth(name)) / 2, lineHeightLocation);
				lineHeightLocation = lineHeightLocation + lineHeight + paddingLine;
			}
			//int maxWidth = 0;
			//int maxHeight = 0;
			int countImage = 0;
			// for each row
			for(int rowIndex = 0; rowIndex<numRows; rowIndex++)
			{
				// for each col
				for(int colIndex = 0; (colIndex<numCols)&&(countImage<loadedImages.size()); colIndex++)
				{
					// draw next image
					int myWidth = loadedImages.get(countImage).getWidth();
					int xLoc = colIndex*(paddingImages + maxWidth + paddingImages);
					if (myWidth<maxWidth)
					{
						xLoc = xLoc + ((maxWidth-myWidth)/2);
					}
					int yLoc = titleHeight + (rowIndex*(paddingImages + maxHeight + paddingImages));
					gp.mG2d.drawImage(loadedImages.get(countImage), null, xLoc, yLoc);
					countImage = countImage + 1;
				}
			}

			///////
			System.out.println("combineLegends write");
			ImageIO.write(gp.mBi, "PNG", new File(theFilenamePath));
			System.out.println("combineLegends done");
			success = true;
		}
		catch (Exception exp)
		{
			exp.printStackTrace(System.out);
			exp.printStackTrace(System.err);
		}
		return success;
	}

	public static boolean writeLegend(String theTitle, String theVersion,
			String[] theLegendNames, String[] theLegendColors,
			int[] theRSymbols,
			String theFilenamePath)
	{
		boolean success = false;
		try
		{
			System.out.println(getVersion());
			System.out.println("writeLegendWithSymbols theTitle = " + theTitle);
			System.out.println("writeLegendWithSymbols theVersion = " + theVersion);
			System.out.println("writeLegendWithSymbols theFilenamePath = " + theFilenamePath);
			if ((null!=theLegendColors)&&(0==theLegendColors.length))
			{
				theLegendColors = null;
			}
			if ((null!=theRSymbols)&&(0==theRSymbols.length))
			{
				theRSymbols = null;
			}
			checkForSameSize(theLegendNames, theLegendColors, theRSymbols);
			/*if (null!=theLegendNames)
			{
				for (String name : theLegendNames)
				{
					System.out.println("writeLegend theLegendNames = " + name);
				}
			}
			if (null!=theLegendColors)
			{
				for (String name : theLegendColors)
				{
					System.out.println("writeLegend theLegendColors = " + name);
				}
			}
			if (null!=theRSymbols)
			{
				for (int name : theRSymbols)
				{
					System.out.println("writeLegend theRSymbols = " + name);
				}
			}*/
			int padding = 5;
			int lineHeightWithPad = LegendJava.getTextHeight() + padding;
			int symbolSize = lineHeightWithPad - padding;
			int maxTextWidth = LegendJava.getTextWidth(theVersion);
			for (String name : theLegendNames)
			{
				int tmp = LegendJava.getTextWidth(name);
				if (tmp > maxTextWidth)
				{
					maxTextWidth = tmp;
				}
			}
			int legendWidth = padding + symbolSize + padding + maxTextWidth + padding;
			ArrayList<String> titleNames = LegendJava.wrapText(theTitle, legendWidth - padding - padding);
			
			for (String titleName : titleNames)
			{
				if (LegendJava.getTextWidth(titleName) > legendWidth)
				{
					legendWidth = LegendJava.getTextWidth(titleName);
				}
			}
			
			LegendJava gp = LegendJava.getGraphics(legendWidth,
					((lineHeightWithPad * titleNames.size()) + lineHeightWithPad + (lineHeightWithPad * theLegendNames.length) + lineHeightWithPad + lineHeightWithPad));
			gp.mG2d.setPaint(Color.black);
			int lineHeightLocation = lineHeightWithPad - padding;
			for (String name : titleNames)
			{
				gp.mG2d.drawString(name, (legendWidth - LegendJava.getTextWidth(name)) / 2, lineHeightLocation);
				lineHeightLocation = lineHeightLocation + lineHeightWithPad;
			}
			lineHeightLocation = lineHeightLocation + lineHeightWithPad;
			for (int index = 0; index < theLegendNames.length; index++)
			{
				String name = theLegendNames[index];
				int currentLocation = lineHeightLocation - lineHeightWithPad + padding + (padding / 2);
				gp.mG2d.translate(padding, currentLocation);
				if ((null != theLegendColors)&&(null == theRSymbols))
				{
					gp.mG2d.setPaint(Color.decode(theLegendColors[index]));
					gp.mG2d.fillRect(0, 0, symbolSize, symbolSize);
				}
				else if ((null != theLegendColors)&&(null != theRSymbols))
				{
					//Color myColor = Color.decode(theLegendColors[index]);
					//System.out.println("color = " + theLegendColors[index] + " == " + myColor.toString());
					gp.mG2d.setPaint(Color.decode(theLegendColors[index]));
					ArrayList<Shape> shapes = LegendShapes.getShapeFromSymbol(theRSymbols[index], lineHeightWithPad/2);
					for (Shape shape : shapes)
					{
						gp.mG2d.draw(shape);
					}
				}
				else if ((null == theLegendColors)&&(null != theRSymbols))
				{
					ArrayList<Shape> shapes = LegendShapes.getShapeFromSymbol(theRSymbols[index], lineHeightWithPad/2);
					for (Shape shape : shapes)
					{
						gp.mG2d.draw(shape);
					}
				}
				gp.mG2d.translate(-padding, -currentLocation);
				gp.mG2d.setPaint(Color.black);
				gp.mG2d.drawString(name, padding + symbolSize + padding, lineHeightLocation);
				lineHeightLocation = lineHeightLocation + lineHeightWithPad;
			}
			lineHeightLocation = lineHeightLocation + lineHeightWithPad;
			gp.mG2d.drawString(theVersion, legendWidth - padding - LegendJava.getTextWidth(theVersion), lineHeightLocation);
			System.out.println("writeLegendWithSymbols write");
			ImageIO.write(gp.mBi, "PNG", new File(theFilenamePath));
			System.out.println("writeLegendWithSymbols done");
			success = true;
		}
		catch (Exception exp)
		{
			exp.printStackTrace(System.out);
			exp.printStackTrace(System.err);
		}
		return success;
	}

	public static void checkForSameSize(String[] theLegendNames, String[] theLegendColors, int[] theRSymbols) throws Exception
	{
		// this will throw an error if null - it is required
		int legendSize = theLegendNames.length;
		int colorSize = legendSize;
		if (null != theLegendColors)
		{
			System.out.println("Colors is non-null");
			colorSize = theLegendColors.length;
		}
		int symbolSize = legendSize;
		if (null != theRSymbols)
		{
			System.out.println("Symbols is non-null");
			symbolSize = theRSymbols.length;
		}
		if (legendSize != colorSize)
		{
			throw new Exception("Legend and Color sizes received do not match");
		}
		if (legendSize != symbolSize)
		{
			throw new Exception("Legend and Symbol sizes received do not match");
		}
	}
}
