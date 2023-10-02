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
from typing import List, Dict, Set
import io
import os
import typing
import csv
import pandas
from mbatch.test.common import read_headers, add_warnings, get_filenames
from mbatch.test.common import index_string_to_dict_str_float, index_string_to_dict_str_int, convert_str_to_int
from mbatch.test.common import convert_int_to_str, add_error, get_newest_file, read_file_to_string
from mbatch.index.index_original_data import OriginalData, read_json_original_data


# pylint: disable=too-many-instance-attributes,too-few-public-methods
class VisualIndexElementBase:
    """
    Class to hold index data for visualizer datasets
    MEMBER VALUES
    """
    # do not set method variables, as they should be initialized in the init function
    # camel-case done to match headers in index file
    m_path_results: str
    m_path_data: str
    m_id: str
    m_files: str
    m_source: str
    m_program: str
    m_project: str
    m_category: str
    m_platform: str
    m_data: str
    m_details: str
    m_data_version: str
    m_test_version: str
    # number of samples
    samples_matrix: int
    samples_mutations: int
    # number of features
    features_matrix: int
    features_mutations: int
    # percent of each batch type which is unknown
    batch_type_unk_pct: Dict[str, float]
    # single batch counts
    single_batch_count: Dict[str, int]
    # number of unique batches
    unique_batch_count: Dict[str, int]
    # strongly (>=50) batch types
    corr_batch_types: Dict[str, int]
    # batch type count
    batch_type_count: int
    # job type, such as Original or Corrected
    job_type: str
    #
    pipeline_status: str

    # pylint: disable=too-many-arguments,too-many-locals
    def __init__(self: 'VisualIndexElementBase',
                 the_path_results: str, the_path_data: str,
                 the_id: str, the_files: str,
                 the_source: str, the_program: str,
                 the_project: str, the_category: str,
                 the_platform: str, the_data: str,
                 the_details: str, the_data_version: str,
                 the_test_version: str, the_samples_matrix: int,
                 the_samples_mutations: int, the_features_matrix: int,
                 the_features_mutations: int, the_batch_type_unk_pct: Dict[str, float],
                 the_single_batch: Dict[str, int], the_unique_batch_count: Dict[str, int],
                 the_corr_batch_types: Dict[str, int], the_batch_type_count: int,
                 the_pipeline_status: str, the_job_type: str) -> None:
        """
        Members described at class level
        :param the_path_results:
        :param the_path_data:
        :param the_id:
        :param the_files:
        :param the_source:
        :param the_program:
        :param the_project:
        :param the_category:
        :param the_platform:
        :param the_data:
        :param the_details:
        :param the_data_version:
        :param the_test_version:
        :param the_samples_matrix:
        :param the_samples_mutations:
        :param the_features_matrix:
        :param the_features_mutations:
        :param the_batch_type_unk_pct:
        :param the_single_batch:
        :param the_unique_batch_count:
        :param the_corr_batch_types:
        :param the_batch_type_count:
        :param the_pipeline_status:
        :param the_job_type
        """
        self.m_path_results = the_path_results
        self.m_path_data = the_path_data
        self.m_id = the_id
        self.m_files = the_files
        self.m_source = the_source
        self.m_program = the_program
        self.m_project = the_project
        self.m_category = the_category
        self.m_platform = the_platform
        self.m_data = the_data
        self.m_details = the_details
        self.m_data_version = the_data_version
        self.m_test_version = the_test_version
        self.samples_matrix = the_samples_matrix

        self.samples_mutations = the_samples_mutations
        self.features_matrix = the_features_matrix

        self.features_mutations = the_features_mutations
        self.batch_type_unk_pct = the_batch_type_unk_pct

        self.single_batch_count = the_single_batch
        self.unique_batch_count = the_unique_batch_count

        self.corr_batch_types = the_corr_batch_types
        self.batch_type_count = the_batch_type_count

        self.job_type = the_job_type

        self.pipeline_status = the_pipeline_status
    # pylint: enable=too-many-arguments,too-many-locals

    def get_key(self: 'VisualIndexElementBase') -> str:
        """
        Make key for this portion of index
        :return: concatenate a bunch of member values
        """
        return f"{self.m_source}--{self.m_program}--{self.m_project}--{self.m_category}--{self.m_platform}--{self.m_data}--{self.m_details}--{self.m_data_version}--{self.m_test_version}--{self.job_type}"

    def get_str_dict_batch_type_unk_pct(self: 'VisualIndexElementBase') -> str:
        """
        convert the batch_type_unk_pct: Dict[str, float] into a parseable string
        for writing to the index for display in the GUi
        :return:  parseable string for dict
        """
        result_list: List[str] = []
        key: str
        value: float
        for key, value in self.batch_type_unk_pct.items():
            if value > 0.0:
                # export as percent values
                # that is 20.5 NOT .205
                result_list.append('"' + key + ":" + str(value) + '"')
        return ", ".join(result_list)

    def get_str_dict_single_batch(self: 'VisualIndexElementBase') -> str:
        """
        convert the single_batch_count: Dict[str, int] into a parseable string
        for writing to the index for display in the GUi
        :return:  parseable string for dict
        """
        result_list: List[str] = []
        key: str
        value: int
        for key, value in self.single_batch_count.items():
            result_list.append('"' + key + ":" + str(value) + '"')
        return ", ".join(result_list)

    def get_str_dict_unique_batch_counts(self: 'VisualIndexElementBase') -> str:
        """
        parseable list of batch types
        :return:  parseable string of batch types
        """
        result_list: List[str] = []
        key: str
        value: int
        for key, value in self.unique_batch_count.items():
            result_list.append('"' + key + ":" + str(value) + '"')
        return ",".join(result_list)

    def get_str_dict_corr_batch_types(self: 'VisualIndexElementBase') -> str:
        """
        convert the corr_batch_types: Dict[str, int] into a parseable string
        for writing to the index for display in the GUi
        :return:  parseable string for dict
        """
        result_list: List[str] = []
        key: str
        value: float
        for key, value in self.corr_batch_types.items():
            result_list.append('"' + key + ":" + str(value) + '"')
        return ", ".join(result_list)

    # noinspection DuplicatedCode
    def write_element(self: 'VisualIndexElementBase', the_out: io.TextIOWrapper) -> None:
        """
        Write the self element to the_out
        :param the_out: write to this stream
        :return: nothing
        """
        # PathResults
        # PathData
        the_out.write(self.m_path_results)
        the_out.write("\t")
        the_out.write(self.m_path_data)
        the_out.write("\t")
        # ID
        # Files
        the_out.write(self.m_id)
        the_out.write("\t")
        the_out.write(self.m_files)
        the_out.write("\t")
        # Source
        # Program
        the_out.write(self.m_source)
        the_out.write("\t")
        the_out.write(self.m_program)
        the_out.write("\t")
        # Project
        # Category
        the_out.write(self.m_project)
        the_out.write("\t")
        the_out.write(self.m_category)
        the_out.write("\t")
        # Platform
        # Data
        the_out.write(self.m_platform)
        the_out.write("\t")
        the_out.write(self.m_data)
        the_out.write("\t")
        # Details
        # DataVersion
        the_out.write(self.m_details)
        the_out.write("\t")
        the_out.write(self.m_data_version)
        the_out.write("\t")
        # TestVersion
        # samples_matrix
        the_out.write(self.m_test_version)
        the_out.write("\t")
        the_out.write(convert_int_to_str(self.samples_matrix))
        the_out.write("\t")
        # samples_mutations
        # features_matrix
        the_out.write(convert_int_to_str(self.samples_mutations))
        the_out.write("\t")
        the_out.write(convert_int_to_str(self.features_matrix))
        the_out.write("\t")
        # features_mutations
        # unknown_batches
        the_out.write(convert_int_to_str(self.features_mutations))
        the_out.write("\t")
        the_out.write(self.get_str_dict_batch_type_unk_pct())
        the_out.write("\t")
        # single_batch
        # batch_unique_cnt
        the_out.write(self.get_str_dict_single_batch())
        the_out.write("\t")
        the_out.write(self.get_str_dict_unique_batch_counts())
        the_out.write("\t")
        # correlated_batch_types
        # batch_type_count
        the_out.write(self.get_str_dict_corr_batch_types())
        the_out.write("\t")
        the_out.write(convert_int_to_str(self.batch_type_count))
        the_out.write("\t")
        # pipeline_Status
        the_out.write(self.pipeline_status)
        the_out.write("\t")
        # pipeline_Status
        the_out.write(self.job_type)
        the_out.write("\n")
    # pylint: enable=too-many-instance-attributes,too-few-public-methods


# pylint: disable=too-few-public-methods
class VisualIndexBase:
    """
    Class to hold list of visualizer datasets
    MEMBER VALUES
    """
    # do not set method variables, as they should be initialized in the init function
    # camel-case done to match headers in index file
    m_ele_dict: Dict[str, VisualIndexElementBase]
    m_filename: str
    m_base_dir: str

    def __init__(self: 'VisualIndexBase', the_filename: str, the_base_dir: str) -> None:
        """
        initialize base values
        :param the_filename: filename for this element
        :param the_base_dir: base directory for paths/etc
        """
        self.m_filename = the_filename
        self.m_ele_dict = {}
        self.m_base_dir = the_base_dir

    def populate_index(self: 'VisualIndexBase') -> None:
        """
        build index from this node of the data
        :return: nothing
        """
        print(f"VisualIndexBase::populate_index m_filename={self.m_filename}", flush=True)
        if os.path.exists(self.m_filename):
            print("populate index from file", flush=True)
            in_file: io.TextIOWrapper
            line: str
            with open(self.m_filename, 'r', encoding='utf-8') as in_file:
                keys: List[str] = read_headers(in_file)
                line: str = in_file.readline().rstrip('\n')
                while '' != line:
                    values: List[str] = line.split("\t")
                    tsv_dict: Dict[str, str] = dict(zip(keys, values))
                    #  the_unique_batch_count: Dict[str, int], the_corr_batch_types: Dict[str, int]
                    newval: VisualIndexElementBase = VisualIndexElementBase(tsv_dict["PathResults"], tsv_dict["PathData"],
                                                                            tsv_dict["ID"], tsv_dict["Files"],
                                                                            tsv_dict["Source"], tsv_dict["Program"],
                                                                            tsv_dict["Project"], tsv_dict["Category"],
                                                                            tsv_dict["Platform"], tsv_dict["Data"],
                                                                            tsv_dict["Details"], tsv_dict["DataVersion"],
                                                                            tsv_dict["TestVersion"],
                                                                            convert_str_to_int(tsv_dict["samples_matrix"]),
                                                                            convert_str_to_int(tsv_dict["samples_mutations"]),
                                                                            convert_str_to_int(tsv_dict["features_matrix"]),
                                                                            convert_str_to_int(tsv_dict["features_mutations"]),
                                                                            index_string_to_dict_str_float(tsv_dict["unknown_batches"]),
                                                                            index_string_to_dict_str_int(tsv_dict["single_batch"]),
                                                                            index_string_to_dict_str_int(tsv_dict["batch_unique_cnt"]),
                                                                            index_string_to_dict_str_int(tsv_dict["correlated_batch_types"]),
                                                                            convert_str_to_int(tsv_dict["batch_type_count"]),
                                                                            tsv_dict["pipeline_status"],
                                                                            tsv_dict["job_type"])
                    self.m_ele_dict[newval.get_key()] = newval
                    line = in_file.readline().rstrip('\n')
        else:
            print("use empty index, since file does not exist", flush=True)

    def write_index_file(self: 'VisualIndexBase') -> None:
        """
        write index file from this element
        :return: nothing
        """
        out_file: io.TextIOWrapper
        write_list: List[VisualIndexElementBase] = list(self.m_ele_dict.values())
        write_list.sort(key=lambda my_element: my_element.m_path_data, reverse=False)
        with open(self.m_filename, 'w', encoding='utf-8') as out_file:
            out_file.write("PathResults\tPathData\t")
            out_file.write("ID\tFiles\t")
            out_file.write("Source\tProgram\t")
            out_file.write("Project\tCategory\t")
            out_file.write("Platform\tData\t")
            out_file.write("Details\tDataVersion\t")
            out_file.write("TestVersion\tsamples_matrix\t")
            out_file.write("samples_mutations\tfeatures_matrix\t")
            out_file.write("features_mutations\tunknown_batches\t")
            out_file.write("single_batch\tbatch_unique_cnt\t")
            out_file.write("correlated_batch_types\tbatch_type_count\t")
            out_file.write("pipeline_status\t")
            out_file.write("job_type\n")
            # \tdata_format\tpca\tboxplot\tcdp\th_c\tdsc\tdiscrete\tngchm\tsuperclust\tumap
            my_entry: VisualIndexElementBase
            for my_entry in write_list:
                my_entry.write_element(out_file)

    # pylint: disable=too-many-arguments,too-many-locals
    def add_entry(self: 'VisualIndexBase',
                  the_path_results: str, the_path_data: str,
                  the_id: str, the_files: str,
                  the_source: str, the_program: str,
                  the_project: str, the_category: str,
                  the_platform: str, the_data: str,
                  the_details: str, the_data_version: str,
                  the_test_version: str, the_samples_matrix: int,
                  the_samples_mutations: int, the_features_matrix: int,
                  the_features_mutations: int, the_batch_type_unk_pct: Dict[str, float],
                  the_single_batch: Dict[str, int], the_unique_batch_count: Dict[str, int],
                  the_corr_batch_types: Dict[str, int], the_batch_type_count: int,
                  the_pipeline_status: str, the_job_type: str) -> None:
        """
        Add a node to the index/data list
        :param the_path_results:
        :param the_path_data:
        :param the_id:
        :param the_files:
        :param the_source:
        :param the_program:
        :param the_project:
        :param the_category:
        :param the_platform:
        :param the_data:
        :param the_details:
        :param the_data_version:
        :param the_test_version:
        :param the_samples_matrix:
        :param the_samples_mutations:
        :param the_features_matrix:
        :param the_features_mutations:
        :param the_batch_type_unk_pct:
        :param the_single_batch:
        :param the_unique_batch_count:
        :param the_corr_batch_types:
        :param the_batch_type_count:
        :param the_pipeline_status:
        :param the_job_type:
        :return:
        """
        newval: VisualIndexElementBase = VisualIndexElementBase(the_path_results, the_path_data,
                                                                the_id, the_files,
                                                                the_source, the_program,
                                                                the_project, the_category,
                                                                the_platform, the_data,
                                                                the_details, the_data_version,
                                                                the_test_version, the_samples_matrix,
                                                                the_samples_mutations, the_features_matrix,
                                                                the_features_mutations, the_batch_type_unk_pct,
                                                                the_single_batch,  the_unique_batch_count,
                                                                the_corr_batch_types, the_batch_type_count,
                                                                the_pipeline_status, the_job_type)
        key: str = newval.get_key()
        if key not in self.m_ele_dict:
            print(f"VisualIndexBase Adding new entry to index {key}", flush=True)
            self.m_ele_dict[key] = newval
        else:
            add_warnings(f"VisualIndexBase *should not happen* Found duplicate index {key}")
    # pylint: enable=too-many-arguments,too-many-locals

    # pylint: disable=too-many-arguments,too-many-locals
    # noinspection DuplicatedCode
    def find_and_add_entries(self: 'VisualIndexBase', the_output_dir: str, the_results_zip: str,
                             the_data_zip: str, the_result_dir: str, the_test_version: str,
                             the_data_path: str) -> None:
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
        print(f"VisualIndexDsc::find_and_add_entries populate for {the_output_dir}", flush=True)
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
        # collect all filenames from the_result_dir/analysis for Files
        file_list: str = get_filenames(os.path.join(the_result_dir, "analysis"))
        success_file: str = os.path.join(the_result_dir, "analysis", "MBATCH_SUCCESS.txt")
        # read the_result_dir/original_data.json for Source Program Project Category Platform Data Details Version
        original_data: OriginalData = read_json_original_data(os.path.join(the_result_dir, "original_data.json"))
        # original data dir
        orig_data_dir: str = os.path.join(os.path.dirname(os.path.dirname(the_data_path)), "original")
        # count number of samples
        samples_matrix: int = get_dataset_samples(the_data_path, 'matrix.tsv')
        samples_mutations: int = get_mut_samples(orig_data_dir, 'mutations.tsv')
        # count number of features
        features_matrix: int = get_dataset_features(the_data_path, 'matrix.tsv')
        features_mutations: int = get_mut_features(orig_data_dir, 'mutations.tsv')
        # percent unknown
        batch_type_unk_pct: Dict[str, float] = get_batch_type_unk_pct(the_data_path, 'batches.tsv')
        # unique batches
        single_batch: Dict[str, int] = get_single_batch_type(the_data_path, 'batches.tsv')
        unique_batch_count: Dict[str, int] = get_batch_type_unique_cnt(the_data_path, 'batches.tsv')
        # strongly correlated batch types
        corr_batch_types: Dict[str, int] = get_batch_type_corr_pct(the_data_path, 'batches.tsv')
        batch_type_count = len(unique_batch_count.items())
        pipeline_status = "FAILURE"
        if os.path.exists(success_file):
            pipeline_status = "SUCCESS"
        self.add_entry(val_path_results, val_path_data, val_id,
                       file_list, original_data.source, original_data.program, original_data.project,
                       original_data.category, original_data.platform, original_data.data,
                       original_data.details, original_data.version,
                       the_test_version, samples_matrix,
                       samples_mutations, features_matrix,
                       features_mutations, batch_type_unk_pct,
                       single_batch, unique_batch_count,
                       corr_batch_types, batch_type_count,
                       pipeline_status, job_type)
    # pylint: enable=too-many-arguments
# pylint: enable=too-few-public-methods,too-many-locals


# pylint: disable=too-many-nested-blocks
def get_dataset_samples(the_data_path: str, the_file: str) -> int:
    """
    Count samples from matrix file
    :param the_data_path: full path to the data
    :param the_file: filename
    :return: integer count of number of samples from the matrix_data.tsv
    """
    result: int = 0
    try:
        read_file: typing.Optional[str] = get_newest_file(the_data_path, the_file)
        if read_file is not None:
            in_file: io.TextIOWrapper
            with open(read_file, 'r', encoding="utf-8") as in_file:
                first_line: str = in_file.readline()
                if len(first_line) > 1:
                    samples: List[str] = first_line.split("\t", -1)
                    result = len(samples)
    except Exception:
        add_error(f"get_dataset_samples error handling {the_file} for {the_data_path}")
        raise
    return result
# pylint: enable=too-many-nested-blocks


# pylint: disable=too-many-nested-blocks
def get_mut_samples(the_data_path: str, the_file: str) -> int:
    """
    Count samples in a mutations file
    :param the_data_path: full path to the data
    :param the_file: filename
    :return: integer count of number of samples from the mutations.tsv file
    """
    result: int = -1
    try:
        read_file: typing.Optional[str] = get_newest_file(the_data_path, the_file)
        if read_file is not None:
            in_file: io.TextIOWrapper
            with open(read_file, 'r', encoding="utf-8") as in_file:
                my_line: str
                sample_ids: Set[str] = set()
                for my_line in in_file:
                    if not my_line.startswith("Tumor_Sample_Barcode"):
                        split_line: List[str] = my_line.split("\t", -1)
                        sample_ids.add(split_line[0])
                result = len(sample_ids)
    except Exception:
        add_error(f"get_mut_samples error handling {the_file} for {the_data_path}")
        raise
    return result
# pylint: enable=too-many-nested-blocks


# pylint: disable=too-many-nested-blocks
def get_dataset_features(the_data_path: str, the_file: str) -> int:
    """
    Count features in matrix file
    :param the_data_path: full path to the data
    :param the_file: filename
    :return: integer count of number of features from the matrix_data.tsv
    """
    result: int = -1
    try:
        read_file: typing.Optional[str] = get_newest_file(the_data_path, the_file)
        if read_file is not None:
            in_file: io.TextIOWrapper
            with open(read_file, 'r', encoding="utf-8") as in_file:
                result = len(list(in_file))
                result -= 1
    except Exception:
        add_error(f"get_dataset_features error handling {the_file} for {the_data_path}")
        raise
    return result
# pylint: enable=too-many-nested-blocks


# pylint: disable=too-many-nested-blocks
def get_mut_features(the_data_path: str, the_file: str) -> int:
    """
    Count features in mutation file
    :param the_data_path: full path to the data
    :param the_file: filename
    :return: integer count of number of features from the mutations.tsv file
    """
    result: int = -1
    try:
        read_file: typing.Optional[str] = get_newest_file(the_data_path, the_file)
        if read_file is not None:
            in_file: io.TextIOWrapper
            with open(read_file, 'r', encoding="utf-8") as in_file:
                my_line: str
                feature_ids: Set[str] = set()
                for my_line in in_file:
                    if not my_line.startswith("Tumor_Sample_Barcode"):
                        split_line: List[str] = my_line.split("\t", -1)
                        # assumes second entry is always the Gene entry
                        # (should be true for std data)
                        feature_ids.add(split_line[1])
                result = len(feature_ids)
    except Exception:
        add_error(f"get_mut_features error handling {the_file} for {the_data_path}")
        raise
    return result
# pylint: enable=too-many-nested-blocks


# pylint: disable=too-many-nested-blocks
def get_batch_type_unique_cnt(the_data_path: str, the_file: str) -> Dict[str, int]:
    """
    count number of unique batches in each batch type
    :param the_data_path: full path to the data
    :param the_file: filename
    :return: dict of number of unique batches for each batch type
    """
    result: Dict[str, int] = {}
    try:
        read_file: typing.Optional[str] = get_newest_file(the_data_path, the_file)
        if read_file is not None:
            # read as string, so string compare works and values may not be numbers
            # need quote none to avoid mis-interpreting number of columns
            data_frame: pandas.DataFrame = pandas.read_csv(read_file, sep="\t", dtype=str, quoting=csv.QUOTE_NONE)
            col: str
            for col in data_frame.columns:
                if not "Sample" == col:
                    result[col] = len(data_frame[col].unique())
    except Exception:
        add_error(f"get_batch_type_unique_cnt error handling {the_file} for {the_data_path}")
        raise
    return result
# pylint: enable=too-many-nested-blocks


# pylint: disable=too-many-nested-blocks
def get_single_batch_type(the_data_path: str, the_file: str) -> Dict[str, int]:
    """
    count number of single batch types
    :param the_data_path: full path to the data
    :param the_file: filename
    :return: dict of number of unique batches for each batch type
    """
    result: Dict[str, int] = {}
    try:
        read_file: typing.Optional[str] = get_newest_file(the_data_path, the_file)
        if read_file is not None:
            # read as string, so string compare works and values may not be numbers
            # need quote none to avoid mis-interpreting number of columns
            data_frame: pandas.DataFrame = pandas.read_csv(read_file, sep="\t", dtype=str, quoting=csv.QUOTE_NONE)
            col: str
            for col in data_frame.columns:
                if not "Sample" == col:
                    val: int = len(data_frame[col].unique())
                    if val<=1:
                        result[col] = val
    except Exception:
        add_error(f"get_batch_type_unique_cnt error handling {the_file} for {the_data_path}")
        raise
    return result
# pylint: enable=too-many-nested-blocks


# pylint: disable=too-many-nested-blocks
def get_batch_type_unk_pct(the_data_path: str, the_file: str) -> Dict[str, float]:
    """
    percent of each batch type that is unknown
    :param the_data_path: full path to the data
    :param the_file: filename
    :return: dict of percent for each batch type that is unknown
    """
    result: Dict[str, float] = {}
    try:
        read_file: typing.Optional[str] = get_newest_file(the_data_path, the_file)
        if read_file is not None:
            # read as string, so string compare works and values may not be numbers
            # need quote none to avoid mis-interpreting number of columns
            data_frame: pandas.DataFrame = pandas.read_csv(read_file, sep="\t", dtype=str,
                                                           quoting=csv.QUOTE_NONE)
            col: str
            for col in data_frame.columns:
                if not "Sample" == col:
                    # percent unknown
                    if "Unknown" in data_frame[col].unique():
                        unk: int = data_frame[col].value_counts()["Unknown"]
                        cnt: int = len(data_frame[col])
                        result[col] = round(((unk / cnt)*100.0))
                    else:
                        result[col] = 0
    except Exception:
        add_error(f"get_batch_type_unk_pct error handling {the_file} for {the_data_path}")
        raise
    return result
# pylint: enable=too-many-nested-blocks


# pylint: disable=too-many-nested-blocks
def get_batch_type_corr_pct(the_data_path: str, the_file: str) -> Dict[str, int]:
    """
    do pairwise correlation among batch types
    :param the_data_path: full path to the data
    :param the_file: filename
    :return: dict of correlation for batch types with >= 50%
    """
    result: Dict[str, int] = {}
    try:
        read_file: typing.Optional[str] = get_newest_file(the_data_path, the_file)
        if read_file is not None:
            # read as string, so string compare works and values may not be numbers
            # need quote none to avoid mis-interpreting number of columns
            data_frame: pandas.DataFrame = pandas.read_csv(read_file, sep="\t",
                                                           dtype=str,
                                                           quoting=csv.QUOTE_NONE)
            df_corr: pandas.DataFrame = data_frame.apply(lambda x:
                                                         pandas.factorize(x)[0]).corr(
                method='pearson', min_periods=1)
            col_name: str
            df_series: pandas.Series
            for col_name, df_series in df_corr.iterrows():
                row_name: str
                cor_value: float
                for row_name, cor_value in df_series.items():
                    if row_name != col_name:
                        if cor_value >= 0.5:
                            ind_list: List[str] = sorted([row_name, col_name])
                            ind_str: str = "+".join(ind_list)
                            result[ind_str] = round(cor_value * 100.0)
    except Exception:
        add_error(f"get_batch_type_corr_pct error handling {the_file} for {the_data_path}")
        raise
    return result
# pylint: enable=too-many-nested-blocks
