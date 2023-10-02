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

# TODO: make this more object oriented with shared code
from typing import List, Dict
import io
import os
from mbatch.test.common import read_headers, add_warnings, get_newest_file, read_file_to_string
from mbatch.index.index_original_data import OriginalData, read_json_original_data


# pylint: disable=too-many-instance-attributes,too-few-public-methods
class VisualIndexElementDsc:
    """
    Class to hold index data for visualizer datasets
    MEMBER VALUES
    """
    # do not set method variables, as they should be initialized in the init function
    # camel-case done to match headers in index file
    m_path_results: str
    m_path_data: str
    m_id: str
    m_source: str
    m_program: str
    m_project: str
    m_category: str
    m_platform: str
    m_data: str
    m_details: str
    m_data_version: str
    m_test_version: str
    # specific to this index
    m_analysis_path: str
    m_overall_dsc: str
    m_overall_dsc_pvalue: str
    m_dsc_1_2: str
    m_dsc_1_2_pvalue: str
    m_dsc_1_3: str
    m_dsc_1_3_pvalue: str
    m_dsc_1_4: str
    m_dsc_1_4_pvalue: str
    m_dsc_2_3: str
    m_dsc_2_3_pvalue: str
    m_dsc_2_4: str
    m_dsc_2_4_pvalue: str
    m_dsc_3_4: str
    m_dsc_3_4_pvalue: str
    #
    job_type: str

    # pylint: disable=too-many-arguments,too-many-locals
    def __init__(self: 'VisualIndexElementDsc', the_path_results: str, the_path_data: str, the_id: str,
                 the_source: str, the_program: str, the_project: str, the_category: str,
                 the_platform: str, the_data: str, the_details: str, the_data_version: str, the_test_version: str,
                 the_analysis_path: str, the_overall_dsc: str, the_overall_dsc_pvalue: str,
                 the_dsc_1_2: str, the_dsc_1_2_pvalue: str, the_dsc_1_3: str, the_dsc_1_3_pvalue: str,
                 the_dsc_1_4: str, the_dsc_1_4_pvalue: str, the_dsc_2_3: str, the_dsc_2_3_pvalue: str,
                 the_dsc_2_4: str, the_dsc_2_4_pvalue: str, the_dsc_3_4: str, the_dsc_3_4_pvalue: str,
                 the_job_type: str) -> None:
        """
        Members described at class level
        :param the_path_results:
        :param the_path_data:
        :param the_id:
        :param the_source:
        :param the_program:
        :param the_project:
        :param the_category:
        :param the_platform:
        :param the_data:
        :param the_details:
        :param the_data_version:
        :param the_test_version:
        :param the_job_type:
        """
        self.m_path_results = the_path_results
        self.m_path_data = the_path_data
        self.m_id = the_id
        self.m_source = the_source
        self.m_program = the_program
        self.m_project = the_project
        self.m_category = the_category
        self.m_platform = the_platform
        self.m_data = the_data
        self.m_details = the_details
        self.m_data_version = the_data_version
        self.m_test_version = the_test_version
        # specific to this index
        self.m_analysis_path = the_analysis_path
        self.m_overall_dsc = the_overall_dsc
        self.m_overall_dsc_pvalue = the_overall_dsc_pvalue
        self.m_dsc_1_2 = the_dsc_1_2
        self.m_dsc_1_2_pvalue = the_dsc_1_2_pvalue
        self.m_dsc_1_3 = the_dsc_1_3
        self.m_dsc_1_3_pvalue = the_dsc_1_3_pvalue
        self.m_dsc_1_4 = the_dsc_1_4
        self.m_dsc_1_4_pvalue = the_dsc_1_4_pvalue
        self.m_dsc_2_3 = the_dsc_2_3
        self.m_dsc_2_3_pvalue = the_dsc_2_3_pvalue
        self.m_dsc_2_4 = the_dsc_2_4
        self.m_dsc_2_4_pvalue = the_dsc_2_4_pvalue
        self.m_dsc_3_4 = the_dsc_3_4
        self.m_dsc_3_4_pvalue = the_dsc_3_4_pvalue
        self.job_type = the_job_type
    # pylint: enable=too-many-arguments,too-many-locals

    def get_key(self: 'VisualIndexElementDsc') -> str:
        """
        create a key based on all members
        :return: key formed from member variables
        """
        return f"{self.m_source}--{self.m_program}--{self.m_project}--{self.m_category}--{self.m_platform}--{self.m_data}--{self.m_details}--{self.m_data_version}--{self.m_test_version}--{self.m_analysis_path}--{self.job_type}"

    # pylint: disable=too-many-statements
    # noinspection DuplicatedCode
    def write_element(self: 'VisualIndexElementDsc', the_out: io.TextIOWrapper) -> None:
        """
        write self element to given stream
        :param the_out: output stream for element
        :return: nothing
        """
        the_out.write(self.m_path_results)
        the_out.write("\t")
        the_out.write(self.m_path_data)
        the_out.write("\t")
        the_out.write(self.m_id)
        the_out.write("\t")
        the_out.write(self.m_source)
        the_out.write("\t")
        the_out.write(self.m_program)
        the_out.write("\t")
        the_out.write(self.m_project)
        the_out.write("\t")
        the_out.write(self.m_category)
        the_out.write("\t")
        the_out.write(self.m_platform)
        the_out.write("\t")
        the_out.write(self.m_data)
        the_out.write("\t")
        the_out.write(self.m_details)
        the_out.write("\t")
        the_out.write(self.m_data_version)
        the_out.write("\t")
        the_out.write(self.m_test_version)
        the_out.write("\t")
        # specific to this index
        the_out.write(self.m_analysis_path)
        the_out.write("\t")
        the_out.write(self.m_overall_dsc)
        the_out.write("\t")
        the_out.write(self.m_overall_dsc_pvalue)
        the_out.write("\t")
        the_out.write(self.m_dsc_1_2)
        the_out.write("\t")
        the_out.write(self.m_dsc_1_2_pvalue)
        the_out.write("\t")
        the_out.write(self.m_dsc_1_3)
        the_out.write("\t")
        the_out.write(self.m_dsc_1_3_pvalue)
        the_out.write("\t")
        the_out.write(self.m_dsc_1_4)
        the_out.write("\t")
        the_out.write(self.m_dsc_1_4_pvalue)
        the_out.write("\t")
        the_out.write(self.m_dsc_2_3)
        the_out.write("\t")
        the_out.write(self.m_dsc_2_3_pvalue)
        the_out.write("\t")
        the_out.write(self.m_dsc_2_4)
        the_out.write("\t")
        the_out.write(self.m_dsc_2_4_pvalue)
        the_out.write("\t")
        the_out.write(self.m_dsc_3_4)
        the_out.write("\t")
        the_out.write(self.m_dsc_3_4_pvalue)
        the_out.write("\t")
        the_out.write(self.job_type)
        the_out.write("\n")
    # pylint: enable=too-many-statements
# pylint: enable=too-many-instance-attributes,too-few-public-methods


# pylint: disable=too-few-public-methods
class VisualIndexDsc:
    """
    Class to hold list of visualizer datasets
    MEMBER VALUES
    """
    # do not set method variables, as they should be initialized in the init function
    # camel-case done to match headers in index file
    m_ele_dict: Dict[str, VisualIndexElementDsc]
    m_filename: str
    m_base_dir: str

    def __init__(self: 'VisualIndexDsc', the_filename: str, the_base_dir: str) -> None:
        """
        initialize base values
        :param the_filename: filename for this element
        :param the_base_dir: base directory for paths/etc
        """
        self.m_filename = the_filename
        self.m_ele_dict = {}
        self.m_base_dir = the_base_dir

    def populate_index(self: 'VisualIndexDsc') -> None:
        """
        build index from this node of the data
        :return: nothing
        """
        print(f"VisualIndexDsc::populate_index m_filename={self.m_filename}", flush=True)
        if os.path.exists(self.m_filename):
            print("VisualIndexDsc::populate_index populate index from file", flush=True)
            in_file: io.TextIOWrapper
            line: str
            with open(self.m_filename, 'r', encoding='utf-8') as in_file:
                keys: List[str] = read_headers(in_file)
                line: str = in_file.readline().rstrip('\n')
                while '' != line:
                    values: List[str] = line.split("\t")
                    tsv_dict: Dict[str, str] = dict(zip(keys, values))
                    newval: VisualIndexElementDsc = VisualIndexElementDsc(tsv_dict["PathResults"], tsv_dict["PathData"],
                                                                          tsv_dict["ID"], tsv_dict["Source"],
                                                                          tsv_dict["Program"], tsv_dict["Project"],
                                                                          tsv_dict["Category"], tsv_dict["Platform"],
                                                                          tsv_dict["Data"], tsv_dict["Details"],
                                                                          tsv_dict["DataVersion"], tsv_dict["TestVersion"],
                                                                          tsv_dict["AnalysisPath"],
                                                                          tsv_dict["Overall-DSC"], tsv_dict["Overall-DSC-pvalue"],
                                                                          tsv_dict["DSC(1,2)"], tsv_dict["DSC-pvalue(1,2)"],
                                                                          tsv_dict["DSC(1,3)"], tsv_dict["DSC-pvalue(1,3)"],
                                                                          tsv_dict["DSC(1,4)"], tsv_dict["DSC-pvalue(1,4)"],
                                                                          tsv_dict["DSC(2,3)"], tsv_dict["DSC-pvalue(2,3)"],
                                                                          tsv_dict["DSC(2,4)"], tsv_dict["DSC-pvalue(2,4)"],
                                                                          tsv_dict["DSC(3,4)"], tsv_dict["DSC-pvalue(3,4)"],
                                                                          tsv_dict["job_type"])
                    self.m_ele_dict[newval.get_key()] = newval
                    line = in_file.readline().rstrip('\n')
        else:
            print("VisualIndexDsc::populate_index use empty index, since file does not exist", flush=True)

    def write_index_file(self: 'VisualIndexDsc') -> None:
        """
        write index file from this element
        :return: nothing
        """
        out_file: io.TextIOWrapper
        write_list: List[VisualIndexElementDsc] = list(self.m_ele_dict.values())
        write_list.sort(key=lambda my_element: my_element.m_path_data, reverse=False)
        with open(self.m_filename, 'w', encoding='utf-8') as out_file:
            out_file.write("PathResults\tPathData\tID\tSource\tProgram\tProject\tCategory\tPlatform\tData\tDetails\tDataVersion\tTestVersion" +
                           "\tAnalysisPath\tOverall-DSC\tOverall-DSC-pvalue\tDSC(1,2)\tDSC-pvalue(1,2)\tDSC(1,3)\tDSC-pvalue(1,3)" +
                           "\tDSC(1,4)\tDSC-pvalue(1,4)\tDSC(2,3)\tDSC-pvalue(2,3)\tDSC(2,4)\tDSC-pvalue(2,4)\tDSC(3,4)\tDSC-pvalue(3,4)\tjob_type\n")
            my_entry: VisualIndexElementDsc
            for my_entry in write_list:
                my_entry.write_element(out_file)

    # pylint: disable=too-many-arguments,too-many-locals
    def add_entry(self: 'VisualIndexDsc', the_path_results: str, the_path_data: str, the_id: str,
                  the_source: str, the_program: str, the_project: str, the_category: str,
                  the_platform: str, the_data: str, the_details: str, the_data_version: str, the_test_version: str,
                  the_analysis_path: str, the_overall_dsc: str, the_overall_dsc_pvalue: str,
                  the_dsc_1_2: str, the_dsc_1_2_pvalue: str, the_dsc_1_3: str, the_dsc_1_3_pvalue: str,
                  the_dsc_1_4: str, the_dsc_1_4_pvalue: str, the_dsc_2_3: str, the_dsc_2_3_pvalue: str,
                  the_dsc_2_4: str, the_dsc_2_4_pvalue: str, the_dsc_3_4: str, the_dsc_3_4_pvalue: str, the_job_type: str) -> None:
        """
        Add a node to the index/data list
        :param the_path_results:
        :param the_path_data:
        :param the_id:
        :param the_source:
        :param the_program:
        :param the_project:
        :param the_category:
        :param the_platform:
        :param the_data:
        :param the_details:
        :param the_data_version:
        :param the_test_version:
        :param the_analysis_path:
        :param the_overall_dsc:
        :param the_overall_dsc_pvalue:
        :param the_dsc_1_2:
        :param the_dsc_1_2_pvalue:
        :param the_dsc_1_3:
        :param the_dsc_1_3_pvalue:
        :param the_dsc_1_4:
        :param the_dsc_1_4_pvalue:
        :param the_dsc_2_3:
        :param the_dsc_2_3_pvalue:
        :param the_dsc_2_4:
        :param the_dsc_2_4_pvalue:
        :param the_dsc_3_4:
        :param the_dsc_3_4_pvalue:
        :param the_job_type:
        :return:
        """
        newval: VisualIndexElementDsc = VisualIndexElementDsc(the_path_results, the_path_data, the_id,
                                                              the_source, the_program, the_project, the_category,
                                                              the_platform, the_data, the_details, the_data_version,
                                                              the_test_version, the_analysis_path, the_overall_dsc,
                                                              the_overall_dsc_pvalue,
                                                              the_dsc_1_2, the_dsc_1_2_pvalue, the_dsc_1_3, the_dsc_1_3_pvalue,
                                                              the_dsc_1_4, the_dsc_1_4_pvalue, the_dsc_2_3, the_dsc_2_3_pvalue,
                                                              the_dsc_2_4, the_dsc_2_4_pvalue, the_dsc_3_4, the_dsc_3_4_pvalue,
                                                              the_job_type)
        key: str = newval.get_key()
        if key not in self.m_ele_dict:
            print(f"VisualIndexDsc Adding new entry to index {key}", flush=True)
            self.m_ele_dict[key] = newval
        else:
            add_warnings(f"VisualIndexDsc *should not happen* Found duplicate index {key}")
    # pylint: enable=too-many-arguments,too-many-locals

    # pylint: disable=too-many-locals
    # noinspection DuplicatedCode
    def find_and_add_entries(self: 'VisualIndexDsc', the_output_dir: str, the_results_zip: str, the_data_zip: str, the_result_dir: str) -> None:
        """
        find completed data analysis and add it to the entries
        :param the_output_dir:
        :param the_results_zip:
        :param the_data_zip:
        :param the_result_dir:
        :param the_test_version:
        :param the_data_path:
        :return:
        """
        data_file: str = get_newest_file(os.path.join(the_result_dir, "analysis"), "DSCOverview.tsv")
        if data_file is not None:
            print(f"VisualIndexDsc::find_and_add_entries populate from file {data_file}", flush=True)
            # for the_results_zip and the_data_zip, replace the_output_dir with /DAPI/DATA, to get PathResults and PathData
            val_path_results: str = the_results_zip.replace(the_output_dir, self.m_base_dir)
            val_path_data: str = the_data_zip.replace(the_output_dir, self.m_base_dir)
            # read from the_result_dir/source_id.txt for ID
            val_id: str = ""
            in_file: io.TextIOWrapper
            with open(os.path.join(the_result_dir, "source_id.txt"), 'r', encoding="utf-8") as in_file:
                val_id = in_file.read().replace('\n', '')
            # get the job type (Original vs Corrected*)
            version_type_file: str = os.path.join(the_result_dir, "version_type.txt")
            print(f"populate_values_from_files version_type_file={version_type_file}", flush=True)
            job_type: str = read_file_to_string(version_type_file)
            # read the_result_dir/original_data.json for Source Program Project Category Platform Data Details Version
            original_data: OriginalData = read_json_original_data(os.path.join(the_result_dir, "original_data.json"))
            line: str
            with open(data_file, 'r', encoding='utf-8') as in_file:
                keys: List[str] = read_headers(in_file)
                line: str = in_file.readline().rstrip('\n')
                while '' != line:
                    values: List[str] = line.split("\t")
                    tsv_dict: Dict[str, str] = dict(zip(keys, values))
                    test_version: str = tsv_dict["dataset"]
                    test_version = test_version[test_version.index("TEST_") + 5:]
                    self.add_entry(val_path_results, val_path_data, val_id,
                                   original_data.source, original_data.program, original_data.project,
                                   original_data.category, original_data.platform, original_data.data,
                                   original_data.details, original_data.version, test_version,
                                   tsv_dict["dataset"], tsv_dict["Overall-DSC"], tsv_dict["Overall-DSC-pvalue"],
                                   tsv_dict["DSC(1,2)"], tsv_dict["DSC-pvalue(1,2)"], tsv_dict["DSC(1,3)"], tsv_dict["DSC-pvalue(1,3)"],
                                   tsv_dict["DSC(1,4)"], tsv_dict["DSC-pvalue(1,4)"], tsv_dict["DSC(2,3)"], tsv_dict["DSC-pvalue(2,3)"],
                                   tsv_dict["DSC(2,4)"], tsv_dict["DSC-pvalue(2,4)"], tsv_dict["DSC(3,4)"], tsv_dict["DSC-pvalue(3,4)"],
                                   job_type)
                    line = in_file.readline().rstrip('\n')
    # pylint: enable=too-many-locals
# pylint: enable=too-few-public-methods
