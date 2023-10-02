#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright (c) 2011-2022 University of Texas MD Anderson Cancer Center

This program is free software: you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by the Free Software Foundation, either version 2 of
the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.
If not, see <https://www.gnu.org/licenses/>.

MD Anderson Cancer Center Bioinformatics on GitHub <https://github.com/MD-Anderson-Bioinformatics>
MD Anderson Cancer Center Bioinformatics at MDA <https://www.mdanderson.org/research/departments-labs-institutes/departments-divisions/bioinformatics-and-computational-biology.html>
@author: Tod Casasent
"""


from typing import List
import os
from mbatch.test.common import get_sorted_dirs, next_sub_dir_starts_with, get_sorted_files, read_file_to_string


# pylint: disable=too-many-instance-attributes,too-few-public-methods
class MBatchEntry:
    """
    Store information on directory entries - used for making menu entries
    """
    # do not set method variables, as they should be initialized in the init function
    # #####
    # job values
    job_id: str
    job_type: str
    notice: str
    # #####
    # hierarchy level values
    # error log location (empty string if no error)
    error: str
    # Label and dropdown entries for menus
    entry_label: str
    dropdown_label: str
    # used in sub dropdowns
    diagram_type: str
    diagram_image: str  # failover or static image
    legend_image: str  # failover or static image
    # Boxplot diagram_type=box
    box_annotations: str
    box_data: str
    box_histogram: str
    # Correlation Density Plot diagram_type=cdp
    # use failover/static attributes above: public String cdp_diagram_image;
    # Hierarchical Clustering diagram_type=hc
    hc_data: str
    hc_order: str
    # Next Gen Clustered Heatmap diagram_type=ngchm
    ngchm: str
    # PCA-Plus (many-to-many) diagram_type=pca
    # PCA-Plus (one-to-many) diagram_type=pca
    pca_annotations: str
    pca_values: str
    # DSC
    dsc_values: str
    # Supervised Clustering diagram_type=supclu
    # use failover/static attributes above: public String sc_diagram_image;
    # use failover/static attributes above: public String sc_legend_image;
    # umap
    umap_batches: str
    umap_samples: str
    umap_umapdat: str
    # Kruskal Wallis / Dunn's Test diagram_type=discrete
    kwd_kwddata: str
    dropdown_entries: List['MBatchEntry']

    def __init__(self: 'MBatchEntry', the_entry_label: str,
                 the_dropdown_label: str, the_diagram_type: str, the_info_dir: str) -> None:
        """
        init and empty/nan values.
        Members described at class level
        """
        # job values
        self.job_id: str = ""
        self.job_type: str = ""
        self.notice: str = ""
        # hierarchy level values
        self.error = ""
        self.diagram_image = ""  # failover or static image
        self.legend_image = ""  # failover or static image
        self.batch_data = ""
        self.box_annotations = ""
        self.box_data = ""
        self.box_histogram = ""
        self.hc_data = ""
        self.hc_order = ""
        self.ngchm = ""
        self.pca_annotations = ""
        self.pca_values = ""
        self.dsc_values = ""
        self.umap_batches = ""
        self.umap_samples = ""
        self.umap_umapdat = ""
        self.kwd_kwddata = ""
        self.entry_label = the_entry_label
        self.dropdown_label = the_dropdown_label
        self.diagram_type = the_diagram_type
        self.dropdown_entries = []
        # for the title
        self.title = ""
        # call update if needed
        if self.entry_label.startswith("TEST_"):
            self.update_index_from_info(the_info_dir)

    def update_index_from_info(self: 'MBatchEntry', the_info_dir: str) -> None:
        """
        load and set values from the_info_dir/self.entry_label directory
        :param the_info_dir: full path to directory containing TEST_<version> labels
        :return: nothing
        """
        # get job id from job_id.txt
        # job id from BEI or other processing step
        # empty string if file string empty or file does not exist
        job_id_file: str = os.path.join(the_info_dir, self.entry_label, "job_id.txt")
        print(f"update_index_from_info job_id_file={job_id_file}", flush=True)
        job_id: str = read_file_to_string(job_id_file)
        # get version type (such as Original-Analyzed) from version_type.txt
        # type of run being done, may also be BEI-RUN from a BEI run
        # empty string if file string empty or file does not exist
        version_type_file: str = os.path.join(the_info_dir, self.entry_label, "version_type.txt")
        print(f"update_index_from_info version_type_file={version_type_file}", flush=True)
        job_type: str = read_file_to_string(version_type_file)
        # correction notice
        notice: str = ""
        if job_type.startswith("Adjusted") | self.entry_label.endswith("adjusted") :
            notice = "This dataset has been corrected using an automated system without human input. The correction does not imply the presence or absence of batch effects in the original data. The user is solely responsible for assessing batch effects (e.g. by using our assessment tools) and deciding whether or not to use the corrected data, which may or may not have mitigated some useful biological information along with any technical artifacts."
        #
        self.job_id = job_id
        self.job_type = job_type
        self.notice = notice
# pylint: enable=too-many-instance-attributes,too-few-public-methods

# ##################################################################
# ##################################################################


def convert_to_title_file(the_png_filename: str) -> str:
    """
    convert from diagram png file to title txt file.
    Replace _Diagram with _Title
    Replace .PNG (case-insensitive) with .txt by trimming end of string.
    :param the_png_filename: Filname of format xxxxxx_Diagramxxxxx.PNG
    :return: title filename
    """
    the_png_filename = the_png_filename.replace("_Diagram", "_Title")
    # remove last three characters
    the_png_filename = the_png_filename[:-3]
    the_png_filename = the_png_filename + ".txt"
    return the_png_filename


def get_dir_path_from(the_dir: str, the_from: str) -> str:
    """
    Find where the_from string is a directory name in the_dir string,
    and return a new directory string path to that from directory.
    :param the_dir: String with the_from as a directory in it somewhere
    :param the_from: The directory to find
    :return: path down to the_from
    """
    index: int = the_dir.find("/" + the_from + "/")
    sub_path: str = the_dir[index:]
    if not sub_path.endswith("/"):
        sub_path = sub_path + "/"
    return sub_path

# ##################################################################
# ##################################################################


def make_entry_boxplot_diagram(the_dir: str, the_parent: MBatchEntry, the_info_dir: dir) -> None:
    """
    Make the diagram entry.
    :param the_dir: current directory to be investigated
    :param the_parent: parent MBatchEntry object -- add or edit
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: None - populate the_parent and new objects (if any)
    """
    # build Batch Type diagram menu entries
    if os.path.exists(os.path.join(the_dir, "error.log")):
        # error generated
        next_entry: MBatchEntry = MBatchEntry("Diagram", "", "boxplot", the_info_dir)
        dir_path: str = get_dir_path_from(the_dir, "analysis")
        next_entry.error = dir_path + "error.log"
        the_parent.dropdown_entries.append(next_entry)
    else:
        # build results based on available files
        # get list of batch ids from file names
        batch_ids: List[str] = []
        file_list: List[os.DirEntry] = get_sorted_files(the_dir)
        dir_entry: os.DirEntry
        diagram_type: str = ""
        for dir_entry in file_list:
            filename: str = dir_entry.name
            if filename != "BatchData.tsv":
                index: int = filename.rfind("-")
                index += 1
                batch_id: str = filename[index:len(filename)-4]
                batch_ids.append(batch_id)
                if "" == diagram_type:
                    # do not use rfind (reverse find) because Batch Type can have underscore
                    index_a: int = filename.find("_")+1
                    index_b: int = filename.find("_", index_a)
                    diagram_type = filename[index_a:index_b]
        batch_ids = list(set(batch_ids))
        batch_ids.sort()
        ####
        # use batch ids to generate menu entries
        ####
        # analysis is the top level directory (within ZIP-RESULTS)
        # that is added to the archive zip
        dir_path: str = get_dir_path_from(the_dir, "analysis")
        batch_id: str
        for batch_id in batch_ids:
            next_entry: MBatchEntry = MBatchEntry(batch_id, "", "boxplot", the_info_dir)
            next_entry.batch_data = dir_path + "BatchData.tsv"
            next_entry.box_annotations = dir_path + "BoxPlot_" + diagram_type + "_Annotations-" + batch_id + ".tsv"
            next_entry.box_data = dir_path + "BoxPlot_" + diagram_type + "_BoxData-" + batch_id + ".tsv"
            next_entry.box_histogram = dir_path + "BoxPlot_" + diagram_type + "_Histogram-" + batch_id + ".tsv"
            next_entry.diagram_image = dir_path + "BoxPlot_" + diagram_type + "_Diagram-" + batch_id + ".png"
            # read title file for boxplot and set title
            next_entry.title = read_file_to_string(os.path.join(the_dir, "BoxPlot_" + diagram_type + "_Title-" + batch_id + ".txt"))
            the_parent.dropdown_entries.append(next_entry)


def make_entry_boxplot_subdirs(the_dir: str, the_parent: MBatchEntry, the_info_dir: str) -> None:
    """
    Recursive function that drills through directory structure.
    Handles optional DATA and TEST version directories.
    :param the_dir: current directory to be investigated
    :param the_parent: parent MBatchEntry object -- add or edit
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: None - populate the_parent and new objects
    """
    # handles data and test versions, and then the results
    dir_name: str = os.path.basename(the_dir)
    next_entry: MBatchEntry
    do_dirs: bool
    # determine labels for next level of directory/menu
    if next_sub_dir_starts_with(the_dir, "DATA"):
        # add optional data version
        next_entry = MBatchEntry(dir_name, "Data Version", "", the_info_dir)
        do_dirs = True
    elif next_sub_dir_starts_with(the_dir, "TEST"):
        # add optional test version
        next_entry = MBatchEntry(dir_name, "Test Version", "", the_info_dir)
        do_dirs = True
    else:
        # add diagram entry
        next_entry = MBatchEntry(dir_name, "Batch Type", "", the_info_dir)
        do_dirs = False
    # add other subdirectories or diagram info
    if do_dirs:
        # add other subdirectories
        subdir_list: List[os.DirEntry] = get_sorted_dirs(the_dir)
        dir_entry: os.DirEntry
        for dir_entry in subdir_list:
            make_entry_boxplot_subdirs(dir_entry.path, next_entry, the_info_dir)
    else:
        # add diagram data
        make_entry_boxplot_diagram(the_dir, next_entry, the_info_dir)
    the_parent.dropdown_entries.append(next_entry)


def make_entry_boxplot(the_dir: str, the_info_dir: str) -> MBatchEntry:
    """
    Builds MBatchEntry for algorithm directory structure.
    Must match MBatch output.
    :param the_dir: Algorithm directory to be populated
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: populated MBatchEntry object
    """
    my_obj: MBatchEntry = MBatchEntry("Boxplot", "Diagram Type", "", the_info_dir)
    # add Diagram Type entries
    subdir_list: List[os.DirEntry] = get_sorted_dirs(the_dir)
    dir_entry: os.DirEntry
    for dir_entry in subdir_list:
        # second argument, dropdown label, will be Batch Type or Data/Test Version
        make_entry_boxplot_subdirs(dir_entry.path, my_obj, the_info_dir)
    return my_obj

# ##################################################################
# ##################################################################


def make_entry_cdp_subdirs(the_dir: str, the_parent: MBatchEntry, the_info_dir: str) -> None:
    """
    Recursive function that drills through directory structure.
    Handles optional DATA and TEST version directories.
    :param the_dir: current directory to be investigated
    :param the_parent: parent MBatchEntry object -- add or edit
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: None - populate the_parent and new objects
    """
    # handles data and test versions, and then the results
    dir_name: str = os.path.basename(the_dir)
    next_entry: MBatchEntry
    do_dirs: bool
    # determine labels for next level of directory/menu
    if next_sub_dir_starts_with(the_dir, "DATA"):
        # add optional data version
        next_entry = MBatchEntry(dir_name, "Data Version", "", the_info_dir)
        do_dirs = True
    elif next_sub_dir_starts_with(the_dir, "TEST"):
        # add optional test version
        next_entry = MBatchEntry(dir_name, "Test Version", "", the_info_dir)
        do_dirs = True
    else:
        # add diagram entry
        next_entry = MBatchEntry(dir_name, "", "cdp", the_info_dir)
        do_dirs = False
    # add other subdirectories or diagram info
    if do_dirs:
        # add other subdirectories
        subdir_list: List[os.DirEntry] = get_sorted_dirs(the_dir)
        dir_entry: os.DirEntry
        for dir_entry in subdir_list:
            make_entry_cdp_subdirs(dir_entry.path, next_entry, the_info_dir)
    elif os.path.exists(os.path.join(the_dir, "error.log")):
        dir_path: str = get_dir_path_from(the_dir, "analysis")
        next_entry.error = dir_path + "error.log"
    else:
        # add diagram data
        # analysis is the top level directory (within ZIP-RESULTS)
        # that is added to the archive zip
        dir_path: str = get_dir_path_from(the_dir, "analysis")
        next_entry.title = read_file_to_string(os.path.join(the_dir, "CDP_Plot_Data1_Title.txt"))
        next_entry.diagram_image = dir_path + "CDP_Plot_Data1_Diagram.PNG"
    the_parent.dropdown_entries.append(next_entry)


def make_entry_cdp(the_dir: str, the_info_dir: str) -> MBatchEntry:
    """
    Builds MBatchEntry for algorithm directory structure.
    Must match MBatch output.
    :param the_dir: Algorithm directory to be populated
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: populated MBatchEntry object
    """
    # this level the_entry_label is "Correlation Density Plot"
    # this level the_dropdown_label depends
    # eventual the_diagram_type is cdp
    ####
    dropdown_label: str = ""
    diagram_type: str = "cdp"
    if next_sub_dir_starts_with(the_dir, "DATA"):
        dropdown_label = "Data Version"
        diagram_type = ""
    elif next_sub_dir_starts_with(the_dir, "TEST"):
        dropdown_label = "Test Version"
        diagram_type = ""
    my_obj: MBatchEntry = MBatchEntry("Correlation Density Plot", dropdown_label, diagram_type, the_info_dir)
    # add Diagram Type entries
    subdir_list: List[os.DirEntry] = get_sorted_dirs(the_dir)
    dir_entry: os.DirEntry
    for dir_entry in subdir_list:
        # second argument, dropdown label, will be Batch Type or Data/Test Version
        make_entry_cdp_subdirs(dir_entry.path, my_obj, the_info_dir)
    return my_obj

# ##################################################################
# ##################################################################


def make_entry_discrete_subdirs(the_dir: str, the_parent: MBatchEntry, the_info_dir: str) -> None:
    """
    Recursive function that drills through directory structure.
    Handles optional DATA and TEST version directories.
    :param the_dir: current directory to be investigated
    :param the_parent: parent MBatchEntry object -- add or edit
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: None - populate the_parent and new objects
    """
    # handles data and test versions, and then the results
    dir_name: str = os.path.basename(the_dir)
    next_entry: MBatchEntry
    do_dirs: bool
    # determine labels for next level of directory/menu
    if next_sub_dir_starts_with(the_dir, "DATA"):
        # add optional data version
        next_entry = MBatchEntry(dir_name, "Data Version", "", the_info_dir)
        do_dirs = True
    elif next_sub_dir_starts_with(the_dir, "TEST"):
        # add optional test version
        next_entry = MBatchEntry(dir_name, "Test Version", "", the_info_dir)
        do_dirs = True
    else:
        # add diagram entry
        next_entry = MBatchEntry(dir_name, "", "discrete", the_info_dir)
        do_dirs = False
    # add other subdirectories or diagram info
    if do_dirs:
        # add other subdirectories
        subdir_list: List[os.DirEntry] = get_sorted_dirs(the_dir)
        dir_entry: os.DirEntry
        for dir_entry in subdir_list:
            make_entry_discrete_subdirs(dir_entry.path, next_entry, the_info_dir)
    elif os.path.exists(os.path.join(the_dir, "error.log")):
        dir_path: str = get_dir_path_from(the_dir, "analysis")
        next_entry.error = dir_path + "error.log"
    else:
        # add diagram data
        dir_path: str = get_dir_path_from(the_dir, "analysis")
        next_entry.title = read_file_to_string(os.path.join(the_dir, "KW_Dunns_Title.txt"))
        next_entry.diagram_image = dir_path + "KW_Dunns_Diagram.PNG"
        next_entry.kwd_kwddata = dir_path + "KW_Dunns_Diagram.tsv"
    the_parent.dropdown_entries.append(next_entry)


def make_entry_discrete(the_dir: str, the_info_dir: str) -> MBatchEntry:
    """
    Builds MBatchEntry for algorithm directory structure.
    Must match MBatch output.
    :param the_dir: Algorithm directory to be populated
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: populated MBatchEntry object
    """
    my_obj: MBatchEntry = MBatchEntry("Kruskal-Wallis/Dunn's Test", "Batch Type", "", the_info_dir)
    # add Diagram Type entries
    subdir_list: List[os.DirEntry] = get_sorted_dirs(the_dir)
    dir_entry: os.DirEntry
    for dir_entry in subdir_list:
        # second argument, dropdown label, will be Batch Type or Data/Test Version
        make_entry_discrete_subdirs(dir_entry.path, my_obj, the_info_dir)
    return my_obj

# ##################################################################
# ##################################################################


def make_entry_hc_subdirs(the_dir: str, the_parent: MBatchEntry, the_info_dir: str) -> None:
    """
    Recursive function that drills through directory structure.
    Handles optional DATA and TEST version directories.
    :param the_dir: current directory to be investigated
    :param the_parent: parent MBatchEntry object -- add or edit
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: None - populate the_parent and new objects
    """
    # handles data and test versions, and then the results
    dir_name: str = os.path.basename(the_dir)
    next_entry: MBatchEntry
    do_dirs: bool
    # determine labels for next level of directory/menu
    if next_sub_dir_starts_with(the_dir, "DATA"):
        # add optional data version
        next_entry = MBatchEntry(dir_name, "Data Version", "", the_info_dir)
        do_dirs = True
    elif next_sub_dir_starts_with(the_dir, "TEST"):
        # add optional test version
        next_entry = MBatchEntry(dir_name, "Test Version", "", the_info_dir)
        do_dirs = True
    else:
        # add diagram entry
        next_entry = MBatchEntry(dir_name, "", "hc", the_info_dir)
        do_dirs = False
    # add other subdirectories or diagram info
    if do_dirs:
        # add other subdirectories
        subdir_list: List[os.DirEntry] = get_sorted_dirs(the_dir)
        dir_entry: os.DirEntry
        for dir_entry in subdir_list:
            make_entry_hc_subdirs(dir_entry.path, next_entry, the_info_dir)
    elif os.path.exists(os.path.join(the_dir, "error.log")):
        dir_path: str = get_dir_path_from(the_dir, "analysis")
        next_entry.error = dir_path + "error.log"
    else:
        # add diagram data
        # analysis is the top level directory (within ZIP-RESULTS)
        # that is added to the archive zip
        dir_path: str = get_dir_path_from(the_dir, "analysis")
        next_entry.batch_data = dir_path + "BatchData.tsv"
        next_entry.legend_image = dir_path + "HierarchicalClustering_Diagram.png"
        next_entry.diagram_image = dir_path + "HierarchicalClustering_Legend-ALL.png"
        next_entry.title = read_file_to_string(os.path.join(the_dir, "HierarchicalClustering_Title.txt"))
        next_entry.hc_data = dir_path + "HCData.tsv"
        next_entry.hc_order = dir_path + "HCOrder.tsv"
    the_parent.dropdown_entries.append(next_entry)


def make_entry_hc(the_dir: str, the_info_dir: str) -> MBatchEntry:
    """
    Builds MBatchEntry for algorithm directory structure.
    Must match MBatch output.
    :param the_dir: Algorithm directory to be populated
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: populated MBatchEntry object
    """
    # this level the_dropdown_label depends
    # eventual the_diagram_type is hc
    ####
    dropdown_label: str = ""
    diagram_type: str = "hc"
    if next_sub_dir_starts_with(the_dir, "DATA"):
        dropdown_label = "Data Version"
        diagram_type = ""
    elif next_sub_dir_starts_with(the_dir, "TEST"):
        dropdown_label = "Test Version"
        diagram_type = ""
    my_obj: MBatchEntry = MBatchEntry("Hierarchical Clustering", dropdown_label, diagram_type, the_info_dir)
    # add Diagram Type entries
    subdir_list: List[os.DirEntry] = get_sorted_dirs(the_dir)
    dir_entry: os.DirEntry
    for dir_entry in subdir_list:
        # second argument
        make_entry_hc_subdirs(dir_entry.path, my_obj, the_info_dir)
    return my_obj

# ##################################################################
# ##################################################################


def make_entry_ngchm_subdirs(the_dir: str, the_parent: MBatchEntry, the_info_dir: str) -> None:
    """
    Recursive function that drills through directory structure.
    Handles optional DATA and TEST version directories.
    :param the_dir: current directory to be investigated
    :param the_parent: parent MBatchEntry object -- add or edit
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: None - populate the_parent and new objects
    """
    # handles data and test versions, and then the results
    dir_name: str = os.path.basename(the_dir)
    next_entry: MBatchEntry
    do_dirs: bool
    # determine labels for next level of directory/menu
    if next_sub_dir_starts_with(the_dir, "DATA"):
        # add optional data version
        next_entry = MBatchEntry(dir_name, "Data Version", "", the_info_dir)
        do_dirs = True
    elif next_sub_dir_starts_with(the_dir, "TEST"):
        # add optional test version
        next_entry = MBatchEntry(dir_name, "Test Version", "", the_info_dir)
        do_dirs = True
    else:
        # add diagram entry
        next_entry = MBatchEntry(dir_name, "", "ngchm", the_info_dir)
        do_dirs = False
    # add other subdirectories or diagram info
    if do_dirs:
        # add other subdirectories
        subdir_list: List[os.DirEntry] = get_sorted_dirs(the_dir)
        dir_entry: os.DirEntry
        for dir_entry in subdir_list:
            make_entry_ngchm_subdirs(dir_entry.path, next_entry, the_info_dir)
    elif os.path.exists(os.path.join(the_dir, "error.log")) and not os.path.exists(os.path.join(the_dir, "All_ngchm.ngchm.html")):
        dir_path: str = get_dir_path_from(the_dir, "analysis")
        next_entry.error = dir_path + "error.log"
    else:
        # add diagram data
        # analysis is the top level directory (within ZIP-RESULTS)
        # that is added to the archive zip
        dir_path: str = get_dir_path_from(the_dir, "analysis")
        next_entry.title = read_file_to_string(os.path.join(the_dir, "All_ngchm.ngchm_Title.txt"))
        # use .html file -- visualizer looks for .ngchm first
        next_entry.ngchm = dir_path + "All_ngchm.ngchm.html"
    the_parent.dropdown_entries.append(next_entry)


def make_entry_ngchm(the_dir: str, the_info_dir: str) -> MBatchEntry:
    """
    Builds MBatchEntry for algorithm directory structure.
    Must match MBatch output.
    :param the_dir: Algorithm directory to be populated
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: populated MBatchEntry object
    """
    # eventual the_diagram_type is ngchm
    ####
    dropdown_label: str = ""
    diagram_type: str = "ngchm"
    if next_sub_dir_starts_with(the_dir, "DATA"):
        dropdown_label = "Data Version"
        diagram_type = ""
    elif next_sub_dir_starts_with(the_dir, "TEST"):
        dropdown_label = "Test Version"
        diagram_type = ""
    my_obj: MBatchEntry = MBatchEntry("NGCHM", dropdown_label, diagram_type, the_info_dir)
    # add Diagram Type entries
    subdir_list: List[os.DirEntry] = get_sorted_dirs(the_dir)
    dir_entry: os.DirEntry
    for dir_entry in subdir_list:
        make_entry_ngchm_subdirs(dir_entry.path, my_obj, the_info_dir)
    return my_obj

# ##################################################################
# ##################################################################


def make_entry_pca_diagram(the_dir: str, the_parent: MBatchEntry) -> None:
    """
    Make the diagram entry.
    :param the_dir: current directory to be investigated
    :param the_parent: parent MBatchEntry object -- add or edit
    :return: None - populate the_parent and new objects (if any)
    """
    # build Batch Type diagram menu entries
    if os.path.exists(os.path.join(the_dir, "error.log")):
        # error generated
        dir_path: str = get_dir_path_from(the_dir, "analysis")
        the_parent.error = dir_path + "error.log"
    else:
        # build results based on available files
        dir_path: str = get_dir_path_from(the_dir, "analysis")
        the_parent.batch_data = dir_path + "BatchData.tsv"
        the_parent.title = read_file_to_string(os.path.join(the_dir, "PCA_Title.txt"))
        the_parent.diagram_image = dir_path + "PCA-Plus/ALL_Comp1_Comp2_Diagram.png"
        the_parent.legend_image = dir_path + "PCA-Plus/ALL_Comp1_Comp2_Legend-ALL.png"
        the_parent.pca_annotations = dir_path + "PCAAnnotations.tsv"
        the_parent.pca_values = dir_path + "PCAValues.tsv"


def make_entry_pca_subdirs(the_dir: str, the_parent: MBatchEntry, the_info_dir: str) -> None:
    """
    Recursive function that drills through directory structure.
    Handles optional DATA and TEST version directories.
    :param the_dir: current directory to be investigated
    :param the_parent: parent MBatchEntry object -- add or edit
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: None - populate the_parent and new objects
    """
    # handles data and test versions, and then the results
    dir_name: str = os.path.basename(the_dir)
    next_entry: MBatchEntry
    do_dirs: bool
    # determine labels for next level of directory/menu
    if (next_sub_dir_starts_with(the_dir, "ManyToMany")) | (next_sub_dir_starts_with(the_dir, "OneToMany")):
        # add optional data version
        next_entry = MBatchEntry(dir_name, "Diagram Type", "", the_info_dir)
        do_dirs = True
    elif next_sub_dir_starts_with(the_dir, "DATA"):
        # add optional data version
        next_entry = MBatchEntry(dir_name, "Data Version", "", the_info_dir)
        do_dirs = True
    elif next_sub_dir_starts_with(the_dir, "TEST"):
        # add optional test version
        next_entry = MBatchEntry(dir_name, "Test Version", "", the_info_dir)
        do_dirs = True
    else:
        # next_sub_dir_starts_with(the_dir, "PCA-Plus"):
        # add diagram entry
        next_entry = MBatchEntry(dir_name, "", "pca", the_info_dir)
        do_dirs = False
    # add other subdirectories or diagram info
    if do_dirs:
        # add other subdirectories
        subdir_list: List[os.DirEntry] = get_sorted_dirs(the_dir)
        dir_entry: os.DirEntry
        for dir_entry in subdir_list:
            make_entry_pca_subdirs(dir_entry.path, next_entry, the_info_dir)
    else:
        # add diagram data
        make_entry_pca_diagram(the_dir, next_entry)
    the_parent.dropdown_entries.append(next_entry)


def make_entry_pca(the_dir: str, the_info_dir: str) -> MBatchEntry:
    """
    Builds MBatchEntry for algorithm directory structure.
    Must match MBatch output.
    :param the_dir: Algorithm directory to be populated
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: populated MBatchEntry object
    """
    my_obj: MBatchEntry = MBatchEntry("PCA+", "Batch Type", "", the_info_dir)
    # add Diagram Type entries
    subdir_list: List[os.DirEntry] = get_sorted_dirs(the_dir)
    dir_entry: os.DirEntry
    for dir_entry in subdir_list:
        # second argument, dropdown label, will be Diagram Type
        # will either have diagram or Data/Test Version
        make_entry_pca_subdirs(dir_entry.path, my_obj, the_info_dir)
    return my_obj

# ##################################################################
# ##################################################################

# ##################################################################
# ##################################################################


def make_entry_dsc_subdirs(the_dir: str, the_parent: MBatchEntry, the_info_dir: str) -> None:
    """
    Recursive function that drills through directory structure.
    Handles optional DATA and TEST version directories.
    :param the_dir: current directory to be investigated
    :param the_parent: parent MBatchEntry object -- add or edit
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: None - populate the_parent and new objects
    """
    # handles data and test versions, and then the results
    dir_name: str = os.path.basename(the_dir)
    next_entry: MBatchEntry
    do_dirs: bool
    # determine labels for next level of directory/menu
    if next_sub_dir_starts_with(the_dir, "DATA"):
        # add optional data version
        next_entry = MBatchEntry(dir_name, "Data Version", "", the_info_dir)
        do_dirs = True
    elif next_sub_dir_starts_with(the_dir, "TEST"):
        # add optional test version
        next_entry = MBatchEntry(dir_name, "Test Version", "", the_info_dir)
        do_dirs = True
    else:
        # add diagram entry
        next_entry = MBatchEntry(dir_name, "", "dsc", the_info_dir)
        do_dirs = False
    # add other subdirectories or diagram info
    if do_dirs:
        # add other subdirectories
        subdir_list: List[os.DirEntry] = get_sorted_dirs(the_dir)
        dir_entry: os.DirEntry
        for dir_entry in subdir_list:
            make_entry_dsc_subdirs(dir_entry.path, next_entry, the_info_dir)
    elif os.path.exists(os.path.join(the_dir, "error.log")):
        dir_path: str = get_dir_path_from(the_dir, "analysis")
        next_entry.error = dir_path + "error.log"
    else:
        # add diagram data
        # analysis is the top level directory (within ZIP-RESULTS)
        # that is added to the archive zip
        dir_path: str = get_dir_path_from(the_dir, "analysis")
        next_entry.dsc_values = dir_path + "DSCOverview.tsv"
    the_parent.dropdown_entries.append(next_entry)


def make_entry_dsc(the_dir: str, the_info_dir: str) -> MBatchEntry:
    """
    Builds MBatchEntry for algorithm directory structure.
    Must match MBatch output.
    :param the_dir: Algorithm directory to be populated
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: populated MBatchEntry object
    """
    # this level the_entry_label is "Dispersion Seperability Criteria"
    # this level the_dropdown_label depends
    # eventual the_diagram_type is dsc
    ####
    dropdown_label: str = ""
    diagram_type: str = "dsc"
    if next_sub_dir_starts_with(the_dir, "DATA"):
        dropdown_label = "Data Version"
        diagram_type = ""
    elif next_sub_dir_starts_with(the_dir, "TEST"):
        dropdown_label = "Test Version"
        diagram_type = ""
    my_obj: MBatchEntry = MBatchEntry("Dispersion Separability Criteria", dropdown_label, diagram_type, the_info_dir)
    # add Diagram Type entries
    subdir_list: List[os.DirEntry] = get_sorted_dirs(the_dir)
    dir_entry: os.DirEntry
    for dir_entry in subdir_list:
        # second argument, dropdown label, will be Batch Type or Data/Test Version
        make_entry_dsc_subdirs(dir_entry.path, my_obj, the_info_dir)
    return my_obj

# ##################################################################
# ##################################################################


def make_entry_sc_subdirs(the_batch_type: str, the_dir: str, the_parent: MBatchEntry, the_info_dir: str) -> None:
    """
    Recursive function that drills through directory structure.
    Handles optional DATA and TEST version directories.
    :param the_batch_type: Batch Type Name for NGCHM file name
    :param the_dir: current directory to be investigated
    :param the_parent: parent MBatchEntry object -- add or edit
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: None - populate the_parent and new objects
    """
    # handles data and test versions, and then the results
    print(f"make_entry_sc_subdirs the_dir={the_dir}", flush=True)
    dir_name: str = os.path.basename(the_dir)
    next_entry: MBatchEntry
    do_dirs: bool
    # determine labels for next level of directory/menu
    if next_sub_dir_starts_with(the_dir, "DATA"):
        # add optional data version
        next_entry = MBatchEntry(dir_name, "Data Version", "", the_info_dir)
        do_dirs = True
    elif next_sub_dir_starts_with(the_dir, "TEST"):
        # add optional test version
        next_entry = MBatchEntry(dir_name, "Test Version", "", the_info_dir)
        do_dirs = True
    else:
        # add diagram entry
        next_entry = MBatchEntry(dir_name, "", "ngchm", the_info_dir)
        do_dirs = False
    # add other subdirectories or diagram info
    if do_dirs:
        # add other subdirectories
        subdir_list: List[os.DirEntry] = get_sorted_dirs(the_dir)
        dir_entry: os.DirEntry
        for dir_entry in subdir_list:
            print(f"make_entry_sc_subdirs dir_entry.path={dir_entry.path}", flush=True)
            make_entry_sc_subdirs(the_batch_type, dir_entry.path, next_entry, the_info_dir)
    elif os.path.exists(os.path.join(the_dir, "error.log")) and not os.path.exists(os.path.join(the_dir, f"{the_batch_type}_ngchm.ngchm.html")):
        dir_path: str = get_dir_path_from(the_dir, "analysis")
        next_entry.error = dir_path + "error.log"
    else:
        # add diagram data
        # analysis is the top level directory (within ZIP-RESULTS)
        # that is added to the archive zip
        dir_path: str = get_dir_path_from(the_dir, "analysis")
        # use .html file -- visualizer looks for .ngchm first
        next_entry.title = read_file_to_string(os.path.join(the_dir, f"{the_batch_type}_ngchm.ngchm_Title.txt"))
        next_entry.ngchm = dir_path + f"{the_batch_type}_ngchm.ngchm.html"
    the_parent.dropdown_entries.append(next_entry)


def make_entry_sc(the_dir: str, the_info_dir: str) -> MBatchEntry:
    """
    Builds MBatchEntry for algorithm directory structure.
    Must match MBatch output.
    :param the_dir: Algorithm directory to be populated
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: populated MBatchEntry object
    """
    # eventual the_diagram_type is ngchm
    ####
    print(f"make_entry_sc the_dir={the_dir}", flush=True)
    my_obj: MBatchEntry = MBatchEntry("Supervised Clustering", "Batch Type", "", the_info_dir)
    # add Diagram Type entries
    subdir_list: List[os.DirEntry] = get_sorted_dirs(the_dir)
    dir_entry: os.DirEntry
    for dir_entry in subdir_list:
        print(f"make_entry_sc dir_entry.path={dir_entry.path}", flush=True)
        batch_type: str = os.path.basename(dir_entry.path)
        make_entry_sc_subdirs(batch_type, dir_entry.path, my_obj, the_info_dir)
    return my_obj

# ##################################################################
# ##################################################################


def make_entry_umap_subdirs(the_dir: str, the_parent: MBatchEntry, the_info_dir: str) -> None:
    """
    Recursive function that drills through directory structure.
    Handles optional DATA and TEST version directories.
    :param the_dir: current directory to be investigated
    :param the_parent: parent MBatchEntry object -- add or edit
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: None - populate the_parent and new objects
    """
    # handles data and test versions, and then the results
    dir_name: str = os.path.basename(the_dir)
    next_entry: MBatchEntry
    do_dirs: bool
    # determine labels for next level of directory/menu
    if next_sub_dir_starts_with(the_dir, "DATA"):
        # add optional data version
        next_entry = MBatchEntry(dir_name, "Data Version", "", the_info_dir)
        do_dirs = True
    elif next_sub_dir_starts_with(the_dir, "TEST"):
        # add optional test version
        next_entry = MBatchEntry(dir_name, "Test Version", "", the_info_dir)
        do_dirs = True
    else:
        # add diagram entry
        next_entry = MBatchEntry(dir_name, "", "umap", the_info_dir)
        do_dirs = False
    # add other subdirectories or diagram info
    if do_dirs:
        # add other subdirectories
        subdir_list: List[os.DirEntry] = get_sorted_dirs(the_dir)
        dir_entry: os.DirEntry
        for dir_entry in subdir_list:
            make_entry_umap_subdirs(dir_entry.path, next_entry, the_info_dir)
    else:
        if os.path.exists(os.path.join(the_dir, "error.log")):
            dir_path: str = get_dir_path_from(the_dir, "analysis")
            next_entry.error = dir_path + "error.log"
        else:
            # add diagram data
            dir_path: str = get_dir_path_from(the_dir, "analysis")
            next_entry.diagram_image = dir_path + "UMAP_Diagram.png"
            next_entry.title = read_file_to_string(os.path.join(the_dir, "UMAP_Title.txt"))
            next_entry.legend_image = dir_path + "UMAP_Legend.png"
            next_entry.umap_batches = dir_path + "UMAP_Data-batc.tsv"
            next_entry.umap_samples = dir_path + "UMAP_Data-meta.tsv"
            next_entry.umap_umapdat = dir_path + "UMAP_Data-umap.tsv"
    the_parent.dropdown_entries.append(next_entry)


def make_entry_umap(the_dir: str, the_info_dir: str) -> MBatchEntry:
    """
    Builds MBatchEntry for algorithm directory structure.
    Must match MBatch output.
    :param the_dir: Algorithm directory to be populated
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: populated MBatchEntry object
    """
    my_obj: MBatchEntry = MBatchEntry("UMAP", "Batch Type", "", the_info_dir)
    # add Diagram Type entries
    subdir_list: List[os.DirEntry] = get_sorted_dirs(the_dir)
    dir_entry: os.DirEntry
    for dir_entry in subdir_list:
        # second argument, dropdown label, will be Batch Type or Data/Test Version
        make_entry_umap_subdirs(dir_entry.path, my_obj, the_info_dir)
    return my_obj

# ##################################################################
# ##################################################################
