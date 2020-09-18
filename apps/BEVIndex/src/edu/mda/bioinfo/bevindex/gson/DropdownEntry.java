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

package edu.mda.bioinfo.bevindex.gson;

import java.util.TreeSet;

/**
 *
 * @author Tod-Casasent
 */
public class DropdownEntry implements Comparable<DropdownEntry>
{
	public String entry_label;
	public String dropdown_label;
	public TreeSet<DropdownEntry> dropdown_entries;
	// used in sub dropdowns
	public String diagram_type;
	public String diagram_image; // failover or static image
	public String legend_image; // failover or static image
	// Boxplot diagram_type=box
	public String box_annotations;
	public String box_data;
	public String box_histogram;
	// Correlation Density Plot diagram_type=cdp
	// use failover/static attributes above: public String cdp_diagram_image;
	// Hierarchical Clustering diagram_type=hc
	public String hc_data;
	public String hc_order;
	// Next Gen Clustered Heatmap diagram_type=ngchm
	public String ngchm;
	// PCA-Plus (many-to-many) diagram_type=pcamtm
	// PCA-Plus (one-to-many) diagram_type=pcaotm
	public String pca_annotations;
	public String pca_values;
	// DSC
	public String dsc_values;
	// Supervised Clustering diagram_type=supclu
	// use failover/static attributes above: public String sc_diagram_image;
	// use failover/static attributes above: public String sc_legend_image;

	@Override
	public int compareTo(DropdownEntry o)
	{
		return this.entry_label.compareTo(o.entry_label);
	}
	
	public void addDropdownEntry(DropdownEntry theEntry)
	{
		if (null==this.dropdown_entries)
		{
			this.dropdown_entries = new TreeSet<>();
		}
		this.dropdown_entries.add(theEntry);
	}
}
