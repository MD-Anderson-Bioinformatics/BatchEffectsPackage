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


import os
import shutil
import typing
from typing import Tuple, List
from mbatch.index.index_original_data import OriginalData, read_json_original_data
from mbatch.index.index_mbatch import MBatchIndex, make_mbatch_index
from mbatch.index.index_repos import populate_from_existing_repo
from mbatch.test.common import read_file_to_string, get_newest_dir
from mbatch.test.common import timestamp_file_if_exists
from mbatch.pipeline.std_data import StandardizedData


M_AUTOCORRECT_NOTICE: str = "This dataset has been corrected using an automated system without human input. " \
                            "The correction does not imply the presence or absence of batch effects in the original data. " \
                            "The user is solely responsible for assessing batch effects " \
                            "(e.g. by using our assessment tools) and deciding whether or not to use the corrected data, " \
                            "which may or may not have mitigated some useful biological information along with any technical artifacts."


def populate_values_from_files(the_directory: str) -> Tuple[OriginalData, str, str, str, str, str, str]:
    """
    Pulls values from various files and returns them in a tuple.
    original_data.json data (source, etc)
    dataset id from conversion step
    job id from BEI or other processing step
    get data version from data_stamp.txt
    get test version from test_stamp.txt
    get version type (such as Original-Analyzed) from version_type.txt
    autocorrection notice string (if needed)
    :param the_directory: fill path to directory from which to read files
    :return: Tuple of above data
    """
    # get info from original_data.json
    # expect this to have JSON data for source, etc
    # if file string variable is empty or file does not exist, object will be populated with empty strings
    original_data_file: str = os.path.join(the_directory, "original_data.json")
    print(f"populate_values_from_files original_data_file={original_data_file}", flush=True)
    original_data: OriginalData = read_json_original_data(original_data_file)
    # get source id from source_id.txt
    # dataset id from conversion step
    # empty string if file string empty or file does not exist
    print(f"pre populate_values_from_files the_directory={the_directory}", flush=True)
    source_id_file: str = os.path.join(the_directory, "source_id.txt")
    print(f"populate_values_from_files source_id_file={source_id_file}", flush=True)
    source_id: str = read_file_to_string(source_id_file)
    # get job id from job_id.txt
    # job id from BEI or other processing step
    # empty string if file string empty or file does not exist
    job_id_file: str = os.path.join(the_directory, "job_id.txt")
    print(f"populate_values_from_files job_id_file={job_id_file}", flush=True)
    job_id: str = read_file_to_string(job_id_file)
    # get data version from data_stamp.txt
    # This is "DATA_" version from new MBatch/MBatchUtils
    # May also be a random timestamp or string
    data_version_file: str = os.path.join(the_directory, "data_stamp.txt")
    print(f"populate_values_from_files data_version_file={data_version_file}", flush=True)
    data_version: str = read_file_to_string(data_version_file)
    # get test version from test_stamp.txt
    # This is "TEST_" version from new MBatch/MBatchUtils
    # May also be a random timestamp from a BEI run
    test_version_file: str = os.path.join(the_directory, "test_stamp.txt")
    print(f"populate_values_from_files test_version_file={test_version_file}", flush=True)
    test_version: str = read_file_to_string(test_version_file)
    # get version type (such as Original-Analyzed) from version_type.txt
    # type of run being done, may also be BEI-RUN from a BEI run
    # empty string if file string empty or file does not exist
    version_type_file: str = os.path.join(the_directory, "version_type.txt")
    print(f"populate_values_from_files version_type_file={version_type_file}", flush=True)
    job_type: str = read_file_to_string(version_type_file)
    # determine if job needs autocorrection notice
    # non-corrected data will have version_types that start with Original- or Aggregate-
    ac_notice: str = ""
    if (not job_type.startswith("Original")) & (not job_type.startswith("Aggregate")) & (not job_type.startswith("BEI-DATA")) & (not job_type.startswith("MBA-DATA")):
        ac_notice = M_AUTOCORRECT_NOTICE
    # TODO handle corrected data
    return original_data, source_id, data_version, test_version, job_id, job_type, ac_notice


# pylint: disable=too-many-locals,too-many-arguments
def create_index_archive(the_results_dir: str, the_data_dir: str, the_zip_dir: str, the_info_dir: str,
                         the_update_only_flag: bool,
                         the_new_data: typing.Optional[StandardizedData],
                         the_std_list: typing.Optional[List[StandardizedData]]) -> Tuple[str, str]:
    """
    Build index and create zip archive
    :param the_results_dir: directory with MBatch results
    :param the_data_dir: directory with actual data
    :param the_zip_dir: directory in which to place ZIP file
    :param the_info_dir: full path to directory containing TEST_<version> labels
    :param the_update_only_flag: if True, don't check for job id, as this is a update run, for sets that already exist
    :param the_new_data: new data object for latest analysis
    :param the_std_list: dictionary of data objects
    :return: full pathname for ZIP file
    """
    print(f"create_index_archive the_results_dir={the_results_dir}", flush=True)
    print(f"create_index_archive the_data_dir={the_data_dir}", flush=True)
    print(f"create_index_archive the_zip_dir={the_zip_dir}", flush=True)
    original_data: OriginalData
    dataset_id: str
    data_version: str
    test_version: str
    job_id: str
    job_type: str
    ac_notice: str
    newest_dir: str = get_newest_dir(os.path.join(the_results_dir, "info"))
    original_data, dataset_id, data_version, test_version, job_id, job_type, ac_notice = \
        populate_values_from_files(newest_dir)
    print(f"create_index_archive original_data={original_data}", flush=True)
    print(f"create_index_archive dataset_id={dataset_id}", flush=True)
    print(f"create_index_archive data_version={data_version}", flush=True)
    print(f"create_index_archive test_version={test_version}", flush=True)
    print(f"create_index_archive job_id={job_id}", flush=True)
    print(f"create_index_archive job_type={job_type}", flush=True)
    print(f"create_index_archive ac_notice={ac_notice}", flush=True)
    analysis_dir: str = os.path.join(the_results_dir, "analysis")
    # use dataset_id to look for existing ZIPs
    dataset_file_id: str = dataset_id
    if the_new_data is not None:
        if the_std_list is not None:
            print("create_index_archive check for existing repo", flush=True)
            reuse_old_dataset_id: str = populate_from_existing_repo(the_data_dir, the_results_dir, the_new_data, the_std_list, the_update_only_flag)
            if reuse_old_dataset_id != '':
                dataset_file_id = reuse_old_dataset_id
                print(f"create_index_archive reuse_old_dataset_id={dataset_file_id}", flush=True)
            else:
                print(f"create_index_archive new dataset_file_id={dataset_file_id}", flush=True)
    print("create_index_archive call make_mbatch_index", flush=True)
    mbatch_index: MBatchIndex = make_mbatch_index(original_data, dataset_id, analysis_dir, the_info_dir)
    print(f"create_index_archive call after 1 make_mbatch_index {os.path.join(newest_dir, 'index.json')}", flush=True)
    mbatch_index.write_to_json(os.path.join(newest_dir, "index.json"))
    print(f"create_index_archive call after 2 make_mbatch_index {os.path.join(the_results_dir, 'index.json')}", flush=True)
    mbatch_index.write_to_json(os.path.join(the_results_dir, "index.json"))
    # create results zip
    zip_file_results: str = os.path.join(the_zip_dir, dataset_file_id + "-results")
    print(f"create_index_archive zip_file_results={zip_file_results}", flush=True)
    timestamp_file_if_exists(f"{zip_file_results + '.zip'}", test_version)
    shutil.make_archive(zip_file_results, "zip", root_dir=the_results_dir)
    print(f"create_index_archive zip={zip_file_results + '.zip'}", flush=True)
    # create data zip
    zip_file_data: str = os.path.join(the_zip_dir, dataset_file_id + "-data")
    print(f"create_index_archive zip_file_data={zip_file_data}", flush=True)
    print(f"create_index_archive the_data_dir={the_data_dir}", flush=True)
    timestamp_file_if_exists(f"{zip_file_data + '.zip'}", test_version)
    shutil.make_archive(zip_file_data, "zip", root_dir=the_data_dir)
    print(f"create_index_archive zip={zip_file_data + '.zip'}", flush=True)
    print(f"create_index_archive completed 1 {zip_file_results + '.zip'}", flush=True)
    print(f"create_index_archive completed 2 {zip_file_data + '.zip'}", flush=True)
    return f"{zip_file_results}.zip", f"{zip_file_data}.zip"
# pylint: enable=too-many-locals,too-many-arguments
