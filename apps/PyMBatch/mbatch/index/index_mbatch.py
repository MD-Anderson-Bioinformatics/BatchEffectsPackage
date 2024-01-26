#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright (c) 2011-2024 University of Texas MD Anderson Cancer Center

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


import io
import os
from typing import List, Optional
import jsonpickle

from mbatch.index.index_entry import make_entry_volcanoplot
from mbatch.index.index_original_data import OriginalData
from mbatch.index.index_entry import MBatchEntry
from mbatch.index.index_entry import make_entry_boxplot, make_entry_cdp, make_entry_discrete
from mbatch.index.index_entry import make_entry_hc, make_entry_ngchm, make_entry_pca
from mbatch.index.index_entry import make_entry_sc, make_entry_umap, make_entry_dsc
from mbatch.test.common import get_sorted_dirs


# pylint: disable=too-many-instance-attributes,too-many-arguments,too-few-public-methods
class MBatchIndex:
    """
    dataset_id:   dataset id from conversion step
    source:       value from original_data.json, used to populate dropdown menus for dataset
    variant:      value from original_data.json, used to populate dropdown menus for dataset
    project:      value from original_data.json, used to populate dropdown menus for dataset
    category:     value from original_data.json, used to populate dropdown menus for dataset
    platform:     value from original_data.json, used to populate dropdown menus for dataset
    data:         value from original_data.json, used to populate dropdown menus for dataset
    algorithm:    value from original_data.json, used to populate dropdown menus for dataset
    details:      value from original_data.json, used to populate dropdown menus for dataset
    mbatch:       object representing results for algorithm menu
    """
    # do not set method variables, as they should be initialized in the init function
    dataset_id: str
    source: str
    program: str
    project: str
    category: str
    platform: str
    data: str
    details: str
    mbatch: Optional[MBatchEntry]

    def __init__(self: 'MBatchIndex', the_dataset_id: str, the_orig_data: OriginalData) -> None:
        """
        Values with which to initialize this object
        :param the_dataset_id:   dataset id from conversion step
        :param the_orig_data:    original_data.json data (source, variant, etc)
        """
        # init and empty/nan values.
        # Members described at class level
        self.dataset_id = the_dataset_id
        self.source = the_orig_data.source
        self.program = the_orig_data.program
        self.project = the_orig_data.project
        self.category = the_orig_data.category
        self.platform = the_orig_data.platform
        self.data = the_orig_data.data
        self.details = the_orig_data.details
        self.mbatch = None

    def write_to_json(self: 'MBatchIndex', the_filepath: str) -> None:
        """
        Write the object to a JSON file for reading by visualization component
        :param the_filepath: full path plus filename to which to write
        :return: Nothing
        """
        out_file: io.TextIOWrapper
        with open(the_filepath, 'w', encoding='utf-8') as out_file:
            out_file.write(jsonpickle.encode(self, indent=4, unpicklable=False))
# pylint: enable=too-many-instance-attributes,too-many-arguments,too-few-public-methods


# pylint: disable=too-many-arguments
def make_mbatch_index(the_original_data: OriginalData, the_dataset_id: str, the_analysis_dir: str, the_info_dir: str) -> MBatchIndex:
    """
    Completely populated MBatchIndex object, ready to be written to json index.
    :param the_original_data: original_data.json data (source, variant, etc)
    :param the_dataset_id:    dataset id from conversion step
    :param the_analysis_dir:  Directory named 'analysis' with results of MBatch run
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: MBatchIndex object completely populated
    """
    # build base object
    my_obj: MBatchIndex = MBatchIndex(the_dataset_id, the_original_data)
    # TODO handle corrected data
    # build menu options from available results on disk
    my_obj.mbatch = make_top_mbatch_entry(the_analysis_dir, the_info_dir)
    return my_obj
# pylint: enable=too-many-arguments


def make_top_mbatch_entry(the_analysis_dir: str, the_info_dir: str) -> MBatchEntry:
    """
    Make top-lebel MBatchEntry instance representing menu structure
    for all the algorithms in this directory.
    Directory structure must match MBatch generated output.
    Directory must be named 'analysis'
    Will throw an exception if unknown algorithm is found.
    :param the_analysis_dir: Directory named 'analysis' with results of MBatch run
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :return: MBatchEntry instance representing MBatch results.
    """
    # build base object
    my_obj: MBatchEntry = MBatchEntry("", "Algorithm", "", the_info_dir)
    # build menu options from available results on disk
    subdir_list: List[os.DirEntry] = get_sorted_dirs(the_analysis_dir)
    dir_entry: os.DirEntry
    for dir_entry in subdir_list:
        if dir_entry.name == "BoxPlot":
            my_obj.dropdown_entries.append(make_entry_boxplot(dir_entry.path, the_info_dir))
        elif dir_entry.name == "VolcanoPlot":
            my_obj.dropdown_entries.append(make_entry_volcanoplot(dir_entry.path, the_info_dir))
        elif dir_entry.name == "CDP":
            my_obj.dropdown_entries.append(make_entry_cdp(dir_entry.path, the_info_dir))
        elif dir_entry.name == "Discrete":
            my_obj.dropdown_entries.append(make_entry_discrete(dir_entry.path, the_info_dir))
        elif dir_entry.name == "HierarchicalClustering":
            my_obj.dropdown_entries.append(make_entry_hc(dir_entry.path, the_info_dir))
        elif dir_entry.name == "NGCHM":
            my_obj.dropdown_entries.append(make_entry_ngchm(dir_entry.path, the_info_dir))
        elif dir_entry.name == "PCA":
            my_obj.dropdown_entries.append(make_entry_pca(dir_entry.path, the_info_dir))
        elif dir_entry.name == "DSC":
            my_obj.dropdown_entries.append(make_entry_dsc(dir_entry.path, the_info_dir))
        elif dir_entry.name == "SupervisedClustering":
            my_obj.dropdown_entries.append(make_entry_sc(dir_entry.path, the_info_dir))
        elif dir_entry.name == "UMAP":
            my_obj.dropdown_entries.append(make_entry_umap(dir_entry.path, the_info_dir))
        else:
            # signal error, unknown algorithm
            raise ValueError(f"The directory {dir_entry.path} does not represent a known algorithm")
    return my_obj
