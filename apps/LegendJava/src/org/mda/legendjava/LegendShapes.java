/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.mda.legendjava;

import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.util.ArrayList;

/**
 *
 * @author tdcasasent
 */
public class LegendShapes
{

	public static ArrayList<Shape> getShapeFromSymbol(int theSymbol, float theSizeFactor)
	{
		ArrayList<Shape> shapeList = new ArrayList<>();
		switch (theSymbol)
		{
			case 0: //	square
			{
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						new Rectangle(Math.round(theSizeFactor), Math.round(theSizeFactor)),
						0, theSizeFactor/4));
			}
			break;
			case 1: //	circle
			{
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						new Ellipse2D.Float(0, 0, theSizeFactor, theSizeFactor),
						0, theSizeFactor/4));
			}
			break;
			case 2: //	triangle-point-up
			{
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						org.jfree.util.ShapeUtilities.createUpTriangle(theSizeFactor / 2),
						theSizeFactor/2, theSizeFactor/2+theSizeFactor/4));
			}
			break;
			case 3: //	plus sign
			{
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						org.jfree.util.ShapeUtilities.createRegularCross(theSizeFactor / 2, 0.1f),
						theSizeFactor/2, theSizeFactor/2+theSizeFactor/4));
			}
			break;
			case 4: //	X
			{
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						org.jfree.util.ShapeUtilities.createDiagonalCross(theSizeFactor / 2, 0.1f),
						theSizeFactor/2, theSizeFactor/2+theSizeFactor/4));
			}
			break;
			case 5: //	Diamond
			{
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						org.jfree.util.ShapeUtilities.createDiamond(theSizeFactor / 2),
						theSizeFactor/2, theSizeFactor/2+theSizeFactor/4));
			}
			break;
			case 6: //	triangle-point-down
			{
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						org.jfree.util.ShapeUtilities.createDownTriangle(theSizeFactor / 2),
						theSizeFactor/2, theSizeFactor/2+theSizeFactor/4));
			}
			break;
			case 7: //	square with X inside
			{
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						new Rectangle(Math.round(theSizeFactor), Math.round(theSizeFactor)),
						0, theSizeFactor/4));
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						org.jfree.util.ShapeUtilities.createDiagonalCross(theSizeFactor / 2, 0.1f),
						theSizeFactor/2, theSizeFactor/2+theSizeFactor/4));
			}
			break;
			case 8: //	star (four lines crossing)
			{
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						org.jfree.util.ShapeUtilities.createRegularCross(theSizeFactor / 2, 0.1f),
						theSizeFactor/2, theSizeFactor/2+theSizeFactor/4));
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						org.jfree.util.ShapeUtilities.createDiagonalCross(theSizeFactor / 2, 0.1f),
						theSizeFactor/2, theSizeFactor/2+theSizeFactor/4));
			}
			break;
			case 9: //	diamond with plus sign inside
			{
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						org.jfree.util.ShapeUtilities.createDiamond(theSizeFactor / 2),
						theSizeFactor/2, theSizeFactor/2+theSizeFactor/4));
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						org.jfree.util.ShapeUtilities.createRegularCross(theSizeFactor / 2, 0.1f),
						theSizeFactor/2, theSizeFactor/2+theSizeFactor/4));
			}
			break;
			case 10: //	circle with plus sign inside
			{
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						new Ellipse2D.Float(0, 0, theSizeFactor, theSizeFactor),
						0, theSizeFactor/4));
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						org.jfree.util.ShapeUtilities.createRegularCross(theSizeFactor / 2, 0.1f),
						theSizeFactor/2, theSizeFactor/2+theSizeFactor/4));
			}
			break;
			case 11: //	Star of David (3 and 6 overlapping)
			{
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						org.jfree.util.ShapeUtilities.createUpTriangle(theSizeFactor / 2),
						theSizeFactor/2, theSizeFactor/2+theSizeFactor/4-theSizeFactor/6));
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						org.jfree.util.ShapeUtilities.createDownTriangle(theSizeFactor / 2),
						theSizeFactor/2, theSizeFactor/2+theSizeFactor/4+theSizeFactor/6));
			}
			break;
			case 12: //	square with plus sign inside
			{
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						new Rectangle(0, 0, Math.round(theSizeFactor), Math.round(theSizeFactor)),
						0, theSizeFactor/4));
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						org.jfree.util.ShapeUtilities.createRegularCross(theSizeFactor / 2, 0.1f),
						theSizeFactor/2, theSizeFactor/2+theSizeFactor/4));
			}
			break;
			case 13: //	circle with X inside
			{
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						new Ellipse2D.Float(0, 0, theSizeFactor, theSizeFactor),
						0, theSizeFactor/4));
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						org.jfree.util.ShapeUtilities.createDiagonalCross(theSizeFactor / 2, 0.1f),
						theSizeFactor/2, theSizeFactor/2+theSizeFactor/4));
			}
			break;
			case 14: //	square with triangle-point-down inside
			{
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						new Rectangle(Math.round(theSizeFactor), Math.round(theSizeFactor)),
						0, theSizeFactor/4));
				shapeList.add(org.jfree.util.ShapeUtilities.createTranslatedShape(
						org.jfree.util.ShapeUtilities.createDownTriangle(theSizeFactor / 2),
						theSizeFactor/2, theSizeFactor/2+theSizeFactor/4));
			}
			break;
		}

		return shapeList;
	}
}
