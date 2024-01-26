# MBatch Copyright (c) 2011-2024 University of Texas MD Anderson Cancer Center
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
# MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>

# Onnes' color distance algorithm as shown in Java by Danny6 Plughoeft
# (which I converted to R) from http://stackoverflow.com/questions/2103368/color-logic-algorithm
# to select as distinct as possible colors.

colorDistance <- function(theColor1, theColor2)
{
	rgb1 <- col2rgb(theColor1)[,1]
	rgb2 <- col2rgb(theColor2)[,1]
	rmean <- (rgb1[["red"]] + rgb2[["red"]]) / 2;
	r <- rgb1[["red"]] - rgb2[["red"]];
	g <- rgb1[["green"]] - rgb2[["green"]];
	b <- rgb1[["blue"]] - rgb2[["blue"]];
	weightR <- 2 + (rmean/256)
	weightG <- 4.0;
	weightB <- 2 + (255-rmean)/256;
	sqrt((weightR*r*r) + (weightG*g*g) + (weightB*b*b))
}

getColorsDistance <- function(theColor, thePossibleColors)
{
	sapply(thePossibleColors, colorDistance, theColor)
}

colorDistEval <- function(theValues)
{
	min(theValues)
}

mycolors <- function()
{
	as.vector(unlist(sapply(setdiff(colors(), c("white", "black")), function(theColor)
	{
		rgb <- col2rgb(theColor)
		if (sum(c(((rgb["red",]>170)||(rgb["red",]<75)),
						((rgb["green",]>170)||(rgb["green",]<75)),
						((rgb["blue",]>170)||(rgb["blue",]<75))))>=2)
		{
			return(NULL)
		}
		else
		{
			return(theColor)
		}
	})))
	# setdiff(colors(), c("white", "black", "floralwhite", "ghostwhite",
	# 										"gray0", "gray1", "gray2", "gray3", "gray4", "gray5", "gray6", "gray7", "gray8", "gray9",
	# 										"gray10", "gray11", "gray12", "gray13", "gray14", "gray15", "gray16", "gray17", "gray18",
	# 										"gray19", "gray20", "gray21", "gray22", "gray23", "gray24", "gray25",
	# 										"gray95", "gray96", "gray97", "gray98", "gray99", "gray100",
	# 										"grey0", "grey1", "grey2", "grey3", "grey4", "grey5", "grey6", "grey7", "grey8", "grey9",
	# 										"grey10", "grey11", "grey12", "grey13", "grey14", "grey15", "grey16", "grey17", "grey18",
	# 										"grey19", "grey20", "grey21", "grey22", "grey23", "grey24", "grey25",
	# 										"grey95", "grey96", "grey97", "grey98", "grey99", "grey100",
	# 										"ivory", "ivory1", "oldlace",	"snow", "snow1", "whitesmoke",
	# 										"lightyellow", "lightyellow1", "mintcream", "honeydew", "honeydew1",
	# 										"azure", "azure1", "lightcyan", "lightcyan1", "cornsilk", "cornsilk1",
	# 										"thistle1", "lemonchiffon", "lemonchiffon1", "lightgoldenrodyellow",
	# 										"chartreuse", "springgreen2", "seashell", "paleturquoise1", "grey", "gray", "gold",
	# 										"pink", "burlywood2", "aliceblue", "beige"
	# ))
}

# get RGB hex values with rgb(t(col2rgb(foo)), maxColorValue=255)
getColorList <- function(theCount, theUsedColors=c("green"))
{
	colorList <- theUsedColors
	while(length(colorList)<theCount)
	{
		possibleColors <- setdiff(mycolors(), colorList)
		distances <- sapply(colorList, getColorsDistance, possibleColors)
		nextColor <- ""
		greatestValue <- 0
		for(newcolor in rownames(distances))
		{
			mymean <- colorDistEval(distances[newcolor,])
			if (mymean>greatestValue)
			{
				greatestValue <- mymean
				nextColor <- newcolor
			}
		}
		if ("" != nextColor)
		{
  		colorList <- c(colorList, nextColor)
		}
		else
		{
		  return(rep_len(colorList, theCount))
		}
	}
	colorList
}

testColors <- function()
{
	knownDiseases <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD",
										 "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP",
										 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV",
										 "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
										 "TGCT", "THCA", "UCEC", "UCS", "UVM", "FPPP", "THYM")
	#  knownDiseases <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "CNTL", "COAD",
	#                     "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP",
	#                     "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV",
	#                     "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
	#                     "TGCT", "THCA", "UCEC", "UCS", "UVM", "FPPP", "THYM")
	tcgaColors <- getColorsForDiseases(knownDiseases)
	x <- 20*(1:length(tcgaColors))
	y <- 20*(1:length(tcgaColors))
	pie(rep(1,length(tcgaColors)), col = tcgaColors, labels=knownDiseases, radius=1)

}

getColorsForDiseases <- function(theDiseaseList)
{
	# for SARC - replaced bisque with burlywood2
	# for TGCT - replaced khaki with darkkhaki
	# for ESCA - replaced lightgray with slategray
	# for KICH - replaced cadetblue1 with cadetblue3
	knownDiseases <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "CNTL", "COAD",
										 "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP",
										 "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV",
										 "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
										 "TGCT", "THCA", "UCEC", "UCS", "UVM", "FPPP", "THYM")
	knownColors <- c( "darkcyan", "green", "blue", "purple", "firebrick", "midnightblue", "brown",
										"red", "slategray", "darkgreen", "magenta", "cadetblue3", "goldenrod", "violet",
										"grey", "olivedrab", "cyan", "gold", "turquoise", "chocolate", "pink",
										"dodgerblue", "mediumvioletred", "forestgreen", "brown1", "burlywood2", "darkgray", "orange",
										"darkkhaki", "seagreen", "tomato", "sienna", "darkorchid", "chartreuse1", "springgreen")
	returnColors <- knownColors
	names(returnColors) <- knownDiseases
	returnColors <- returnColors[which(names(returnColors) %in% theDiseaseList)]
	newDiseases <- setdiff(theDiseaseList, knownDiseases)
	if (length(newDiseases)>0)
	{
		newColors <- getColorList((length(knownColors)+length(newDiseases)), knownColors)[(1+length(knownColors)):(length(knownColors)+length(newDiseases))]
		newColorSet <- c(returnColors, newColors)
		names(newColorSet) <- c(names(returnColors), newDiseases)
		returnColors <- newColorSet
	}
	returnColors <- returnColors[order(names(returnColors))]
	returnColors
}


getColorListAny <- function(theCount, theUsedColors=c("green"))
{
  #message("getColorListAny theCount=", theCount)
  #message("getColorListAny theUsedColors=", theUsedColors)
  colorList <- theUsedColors
  #message("getColorListAny colorList=", colorList)
  #message("getColorListAny length(colorList)=", length(colorList))
	while(length(colorList)<theCount)
	{
		possibleColors <- setdiff(mycolors(), colorList)
		#message("getColorListAny possibleColors=", paste(possibleColors, collapse = ",", sep=","))
		#message("getColorListAny length(possibleColors)=", length(possibleColors))
		distances <- sapply(colorList, getColorsDistance, possibleColors)
		nextColor <- ""
		greatestValue <- 0
		for(newcolor in rownames(distances))
		{
			mymean <- colorDistEval(distances[newcolor,])
			if (mymean>greatestValue)
			{
				greatestValue <- mymean
				nextColor <- newcolor
			}
		}
		if ("" != nextColor)
		{
		  colorList <- c(colorList, nextColor)
		}
		else
		{
		  return(rep_len(colorList, theCount))
		}
		colorList <- c(colorList, nextColor)
	}
	colorList
}
