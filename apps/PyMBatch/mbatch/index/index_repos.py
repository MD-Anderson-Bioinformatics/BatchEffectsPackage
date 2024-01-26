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
import zipfile
import typing
from typing import List
from mbatch.pipeline.std_data import StandardizedData


def extract_data(the_zip_path: str, the_out_dir: str) -> None:
    """
    Extract updatable results from ZIP to the given directory
    :param the_zip_path: ZIP file to extract
    :param the_out_dir: Directory to which to extract
    :return: Nothing
    """
    zip_file: zipfile.ZipFile
    file_name: str
    # extract versions directory
    with zipfile.ZipFile(the_zip_path, 'r') as zip_file:
        for file_name in zip_file.namelist():
            if file_name.startswith('versions/'):
                zip_file.extract(file_name, the_out_dir)


def extract_results(the_zip_path: str, the_out_dir: str) -> None:
    """
    Extract updatable results from ZIP to the given directory
    :param the_zip_path: ZIP file to extract
    :param the_out_dir: Directory to which to extract
    :return: Nothing
    """
    zip_file: zipfile.ZipFile
    file_name: str
    # extract analysis and info directories
    with zipfile.ZipFile(the_zip_path, 'r') as zip_file:
        for file_name in zip_file.namelist():
            if file_name.startswith('analysis/'):
                zip_file.extract(file_name, the_out_dir)
            if file_name.startswith('info/'):
                zip_file.extract(file_name, the_out_dir)


def populate_from_existing_repo(the_data_directory: str, the_result_directory: str, the_new_data: StandardizedData,
                                the_std_list: List[StandardizedData], the_update_only_flag: bool) -> str:
    """
    Find existing repo, extract existing results and data into new results and data directories,
    then return the dataset id used in the original repo, for use in the new one.
    :param the_data_directory: new data directory in BEI/MBA analysis
    :param the_result_directory: new results directory in BEI/MBA analysis
    :param the_new_data: new/updated StandardizedData object
    :param the_std_list: dictionary of all StandardizedData object entries
    :param the_update_only_flag: if True, don't check for job id, as this is a update run, for sets that already exist
    :return:
    """
    # std_archive	version	data_archive	result_archive
    # search through the_std_list for older versions
    old_data: typing.Optional[StandardizedData] = None
    tmp_data: StandardizedData
    for tmp_data in the_std_list:
        if tmp_data.std_archive == the_new_data.std_archive:
            # if tmp_data if for the same dataset as the new data...
            if old_data is None:
                if (tmp_data.job_id != the_new_data.job_id) | the_update_only_flag:
                    # if no old_data yet, and tmp_data job id is not same as new data, use tmp as old data
                    old_data = tmp_data
    # do not need to search more because path to results and data ZIP is always the same
    dataset_file_id: str = ''
    if old_data is not None:
        dataset_file_id = os.path.basename(old_data.data_archive).split("-")[0]
        extract_data(old_data.data_archive, the_data_directory)
        extract_results(old_data.result_archive, the_result_directory)
    return dataset_file_id
