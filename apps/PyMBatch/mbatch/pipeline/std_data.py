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

import io
import os
import csv
import zipfile
import hashlib
import json
from typing import Optional, Dict, List
import pandas
import requests
from mbatch.test.common import write_tsv_list, read_headers


# must be same name as attributes for StandardizedData class
STD_HEADERS: List[str] = [
    'std_archive', 'version', 'data_archive', 'result_archive', 'job_id', 'job_status'
]


# pylint: disable=too-many-instance-attributes,too-few-public-methods,too-many-arguments
class StandardizedData:
    """
    StandardizedData is a class where each instance represents a particular version of a dataset.
    self.std_archive: str is the sub-path from the the_input_dir down to the archive ZIP
    self.version: str is the version number within the ZIP of the data being referenced
    self.data_archive: str is the full output of directory with ZIP file.
        Path is results dir, run source, path, std_archive, dataset_id + "-data.zip"
    self.result_archive: str is the full output of directory with ZIP file.
        Path is results dir, run source, path, std_archive, dataset_id + "-results.zip"
    self.job_id: str is a string giving the job id created by BEI
    self.job_status: str is the current job status.
        Options are:
            created - job had been created by BEI
            queued - job is queued to run in BEI
            running - job is currently running
            succeeded - job completed successfully
            failed - job completed unsuccessfully
    """
    # declare but do not set member attributes
    std_archive: str
    version: str
    data_archive: str
    result_archive: str
    job_id: str
    job_status: str

    @classmethod
    def from_dict(cls, the_dict: Optional[Dict[str, str]]) -> 'StandardizedData':
        """
        Create a StandardizedData object from a Dictionary
        :param the_dict: dictionary with keys matching attributes of StandardizedData class
        :return: instance of StandardizedData
        """
        return cls(the_dict['std_archive'], the_dict['version'],
                   the_dict['data_archive'], the_dict['result_archive'],
                   the_dict['job_id'], the_dict['job_status'])

    # pylint: disable=too-many-arguments
    def __init__(self: 'StandardizedData', the_std_archive: str, the_version: str,
                 the_data_archive: str, the_result_archive: str,
                 the_job_id: str, the_job_status: str) -> None:
        """
        Initialized values for a StandardizedData object
        :param the_std_archive: str is the sub-path from the the_input_dir down to the archive ZIP
        :param the_version: str is the version number within the ZIP of the data being referenced
        :param the_data_archive: tr is the full output of directory with ZIP file
        :param the_result_archive: str is the full output of directory with ZIP file
        :param the_job_id: str is a string giving the job id created by BEI
        :param the_job_status: str is the current job status
        """
        super().__init__()
        self.std_archive: str = the_std_archive
        self.version: str = the_version
        self.data_archive: str = the_data_archive
        self.result_archive: str = the_result_archive
        self.job_id: str = the_job_id
        # created, queued, running, succeeded, failed
        self.job_status: str = the_job_status
    # pylint: enable=too-many-arguments

    # pylint: disable=too-many-branches
    def needs_log_transform(self: 'StandardizedData', the_run_source: str) -> bool:
        """
        Tells if dataset needs to be log transformed
        :param the_run_source: string source, GDC or MWB
        :return: True if needs log transforming
        """
        log_trans: bool = False
        if 'MWB' == the_run_source:
            log_trans = True
        else:
            if '/ASCAT2/' in self.std_archive:
                log_trans = True
            elif '/Aliquot Ensemble Somatic Variant Merging and Masking/' in self.std_archive:
                log_trans = False
            elif '/STAR - Counts/' in self.std_archive:
                log_trans = False
            elif '/SeSAMe Methylation Beta Estimation/' in self.std_archive:
                log_trans = False
            elif '/Masked Somatic Mutation/' in self.std_archive:
                log_trans = True
            elif '/Gene Expression Quantification/' in self.std_archive:
                log_trans = True
            elif '/BCGSC miRNA Profiling/' in self.std_archive:
                log_trans = True
            elif '/Masked Copy Number Segment/' in self.std_archive:
                log_trans = False
            elif '/Copy Number Segment/' in self.std_archive:
                log_trans = False
            elif '/miRNA-Seq/' in self.std_archive:
                log_trans = True
            elif '/Protein Expression Quantification/' in self.std_archive:
                log_trans = False
            elif '/Gene Level Copy Number/' in self.std_archive:
                log_trans = False
            elif ('/scRNA-Seq/' in self.std_archive) & ('/Differential Gene Expression/' in self.std_archive) & ('/Seurat - 10x Chromium/' in self.std_archive):
                log_trans = False
            else:
                print(f"Unknown datatype {self.std_archive}")
        return log_trans
        # pylint: enable=too-many-branches

    # pylint: disable=too-many-branches
    def needs_replace_na(self: 'StandardizedData', the_run_source: str) -> bool:
        """
        Tells if dataset needs NAs to be replaced with zeros
        :param the_run_source: string source, GDC or MWB
        :return: True if needs NAs to be replaced with zeros
        """
        replace_nas: bool = False
        if 'MWB' == the_run_source:
            replace_nas = False
        else:
            if '/ASCAT2/' in self.std_archive:
                replace_nas = False
            elif '/Aliquot Ensemble Somatic Variant Merging and Masking/' in self.std_archive:
                replace_nas = True
            elif '/STAR - Counts/' in self.std_archive:
                replace_nas = False
            elif '/SeSAMe Methylation Beta Estimation/' in self.std_archive:
                replace_nas = False
            elif '/Masked Somatic Mutation/' in self.std_archive:
                replace_nas = True
            elif '/Gene Expression Quantification/' in self.std_archive:
                replace_nas = False
            elif '/BCGSC miRNA Profiling/' in self.std_archive:
                replace_nas = False
            elif '/Masked Copy Number Segment/' in self.std_archive:
                replace_nas = False
            elif '/Copy Number Segment/' in self.std_archive:
                replace_nas = False
            elif '/miRNA-Seq/' in self.std_archive:
                replace_nas = False
            elif '/Protein Expression Quantification/' in self.std_archive:
                replace_nas = False
            elif '/Gene Level Copy Number/' in self.std_archive:
                replace_nas = False
            elif ('/scRNA-Seq/' in self.std_archive) & ('/Differential Gene Expression/' in self.std_archive) & ('/Seurat - 10x Chromium/' in self.std_archive):
                replace_nas = False
            else:
                print(f"Unknown datatype {self.std_archive}")
        return replace_nas
        # pylint: enable=too-many-branches

    def get_title(self: 'StandardizedData', the_run_version: str) -> str:
        """
        Build title from std_archive string and versions
        :return: title for results
        """
        title: str = self.std_archive[:self.std_archive.rfind("/")]
        title = title.replace('/', ' / ')
        title = f"{title} / DATA_{self.version} / TEST_{the_run_version}"
        return title

    def write_index_row(self: 'StandardizedData', the_out_file: io.TextIOWrapper) -> None:
        """
        Write out this entry to a file,
        Writes a newline.
        :param the_out_file: stream to which to write strings
        :return: None
        """
        data_list: List[str] = []
        for attr in STD_HEADERS:
            data_list.append(getattr(self, attr))
        write_tsv_list(the_out_file, data_list, True, False)

    def has_valid_data_p(self: 'StandardizedData', the_base_zip_dir: str) -> bool:
        """
        Check inside the ZIP archive at the_base_zip_dir/std_archive if version has a matrix.tsv file
        Check there is more than one Sample
        Check there is more than one feature
        :param the_base_zip_dir: base directory for paths to ZIP files
        :return: True if there is a matrix file and it is larger than one in both directions
        """
        is_good_data_p: bool = False
        std_archive_full: str = os.path.join(the_base_zip_dir, self.std_archive)
        # use manual / as it is ZIP file, not OS file
        my_filename: str = f"versions/DATA_{self.version}/matrix.tsv"
        zip_file: zipfile.ZipFile
        read_zip_file: io.TextIOWrapper
        with zipfile.ZipFile(std_archive_full, 'r') as zip_file:
            if my_filename in zip_file.namelist():
                with zip_file.open(my_filename, mode="r") as read_zip_file:
                    # set no index column with index_col=False so my_batches['foo'] returns a column of data
                    my_matrix: pandas.DataFrame = pandas.read_csv(read_zip_file, sep='\t', quoting=csv.QUOTE_NONE,
                                                                  encoding='utf-8', index_col=False, dtype=str)
                    col: int = my_matrix.shape[1]
                    if col > 1:
                        row: int = my_matrix.shape[0]
                        if row > 1:
                            is_good_data_p = True
        return is_good_data_p

    def has_batch_info_p(self: 'StandardizedData', the_base_zip_dir: str) -> bool:
        """
        Check inside the ZIP archive at the_base_zip_dir/std_archive if version has a batches.tsv file
        :param the_base_zip_dir: base directory for paths to ZIP files
        :return: True is archive version has batches.tsv
        """
        has_batch_info: bool = False
        std_archive_full: str = os.path.join(the_base_zip_dir, self.std_archive)
        # use manual / as it is ZIP file, not OS file
        my_filename: str = f"versions/DATA_{self.version}/batches.tsv"
        zip_file: zipfile.ZipFile
        with zipfile.ZipFile(std_archive_full, 'r') as zip_file:
            if my_filename in zip_file.namelist():
                has_batch_info = True
        return has_batch_info

    def read_batch_info(self: 'StandardizedData', the_base_zip_dir: str) -> pandas.DataFrame:
        """
        Read batches.tsv for this dataset version.
        Returns empty DataFrame if not batch information.
        Return DataFrame of batch information. One column is samples (aliquot_barcode) rest are potential batch types
        :param the_base_zip_dir: base directory for paths to ZIP files
        :return: DataFrame of batch information. One column is samples (aliquot_barcode) rest are potential batch types
        """
        my_batches: pandas.DataFrame
        std_archive_full: str = os.path.join(the_base_zip_dir, self.std_archive)
        # use manual / as it is ZIP file, not OS file
        my_filename: str = f"versions/DATA_{self.version}/batches.tsv"
        # Do not type, since pandas.read_csv is not typed properly
        # in that TextIOWrapper works, but the type checker doesn't know it.
        # in_file: io.TextIOWrapper
        zip_file: zipfile.ZipFile
        with zipfile.ZipFile(std_archive_full, 'r') as zip_file:
            with zip_file.open(my_filename, mode="r") as in_file:
                # set no index column with index_col=False so my_batches['foo'] returns a column of data
                my_batches = pandas.read_csv(in_file, sep='\t', quoting=csv.QUOTE_NONE, encoding='utf-8', index_col=False, dtype=str)
                my_batches.replace(r'[\t\n\r\"]', '', regex=True, inplace=True)
        return my_batches

    def create_new_job(self: 'StandardizedData', the_url: str) -> str:
        """
        Create a new job and return the job id.
        :param the_url: BEI URL for make job requests
        :return: job id from BEI
        """
        response: requests = requests.get(the_url + "JOBnew", timeout=60, allow_redirects=True)
        print(f"create_new_job response={response}", flush=True)
        # check response code
        if response.ok:
            self.job_id = response.text
            self.job_status = 'created'
        else:
            raise f"create_new_job error {response.status_code}"
        self.job_status = 'created'
        return self.job_id

    def queue_job(self: 'StandardizedData', the_url: str) -> None:
        """
        Queue the job for this dataset for running by BEI.
        :param the_url: BEI URL for make job requests
        :return: Nothing
        """
        # make HTTP call to run job and change job_status to queued
        response: requests = requests.get(the_url + "JOBRunMBatch?jobId=" + self.job_id, timeout=60, allow_redirects=True)
        print(f"queue_job response={response}", flush=True)
        # check response code
        if response.ok:
            self.job_status = 'queued'
        else:
            raise f"queue_job error {response.status_code}"

    def status_job(self: 'StandardizedData', the_url: str) -> bool:
        """
        Get the current status of the job for this dataset.
        Set status in instance, return True if status changed.
        :param the_url: BEI URL for make job requests
        :return: Set status in instance, return True if status changed.
        """
        changed: bool = False
        if ('succeeded' != self.job_status) & ('failed' != self.job_status):
            response: requests = requests.get(the_url + "JOBstatus?jobId=" + self.job_id, timeout=60, allow_redirects=True)
            print(f"status_job response={response} for job_id={self.job_id}", flush=True)
            # check response code
            if response.ok:
                # returned JSON. Turn into dict and look at "status"
                my_status_str: str = response.text
                my_status_dict: dict = json.loads(my_status_str)
                my_status: str = my_status_dict['status']
                # status_job -> created, running, succeeded, failed
                if 'MBATCHRUN_RUNNING_WAIT' == my_status:
                    if 'running' != self.job_status:
                        self.job_status = 'running'
                        changed = True
                elif 'MBATCHRUN_END_SUCCESS' == my_status:
                    if 'succeeded' != self.job_status:
                        self.job_status = 'succeeded'
                        changed = True
                elif 'MBATCHRUN_END_COMPLETED' == my_status:
                    if 'succeeded' != self.job_status:
                        self.job_status = 'succeeded'
                        changed = True
                elif 'MBATCHRUN_END_FAILURE' == my_status:
                    if 'failed' != self.job_status:
                        self.job_status = 'failed'
                        changed = True
            else:
                raise f"status_job error {response.status_code} for job_id={self.job_id}"
        return changed

    def make_fake_batch_info(self: 'StandardizedData', the_base_zip_dir: str, the_batch_file: str,
                             the_sample_column_name: str) -> None:
        """
        Create a fake batch info file for this dataset.
        :param the_base_zip_dir: Directory base for ZIP files (needs std_archive)
        :param the_batch_file: Full path and filename to which to write created batch information
        :param the_sample_column_name: column string for samples in batch (different for GDC vs MWB)
        :return: Nothing
        """
        # write fake batch info for processing data without batches
        # read first line of matrix and get sample ids
        list_of_samples: List[str] = []
        # no slash at front
        my_filename: str = f"versions/DATA_{self.version}/matrix.tsv"
        my_arch: str = os.path.join(the_base_zip_dir, self.std_archive)
        print(f'make_fake_batch_info my_arch={my_arch}', flush=True)
        print(f'make_fake_batch_info my_filename={my_filename}', flush=True)
        zip_file: zipfile.ZipFile
        with zipfile.ZipFile(my_arch, 'r') as zip_file:
            read_zip_file: io.TextIOWrapper
            # read as bytes for read_headers
            with zip_file.open(my_filename, mode="r") as read_zip_file:
                keys: List[str] = read_headers(read_zip_file)
                tmp_key: str
                for tmp_key in keys:
                    if '' != tmp_key:
                        list_of_samples.append(tmp_key)
        # write batches.tsv with headers aliquot_barcode and "example"
        list_of_samples.sort()
        out_file: io.TextIOWrapper
        with open(the_batch_file, 'w', encoding='utf-8') as out_file:
            out_file.write(f"{the_sample_column_name}\texample\n")
            barcode: str
            for barcode in list_of_samples:
                # write "Unknown" for each batch
                out_file.write(f"{barcode}\tUnknown\n")

    # pylint: disable=too-many-arguments
    def calc_md5(self: 'StandardizedData') -> str:
        """
        Calculate the MD5 sum for a Versioned Standardized Dataset
        :return: MD5 sum from sub-path to ZIP and version
        """
        my_hash: hashlib.md5 = hashlib.md5()
        my_hash.update(self.std_archive.encode('utf-8'))
        my_hash.update(self.version.encode('utf-8'))
        return my_hash.hexdigest()
    # pylint: enable=too-many-arguments

    def build_results_path(self: 'StandardizedData', the_results_dir: str, the_run_source: str) -> str:
        """
        Build path to where results and data archives should be written
        :param the_results_dir: top directory for data/results archives
        :param the_run_source: Source for run (GDC or MWB)
        :return: String for full path, the_results_dir plus the_run_source, plus std_archive with filename removed
        """
        out_path: str = os.path.join(the_results_dir, the_run_source)
        out_path = os.path.join(out_path, self.std_archive)
        out_path = os.path.dirname(out_path)
        return out_path
# pylint: enable=too-many-instance-attributes,too-few-public-methods,too-many-arguments


# ########################################################
# writing files
# ########################################################

def write_std_pipeline_index(the_index_file: str, the_entries: List[StandardizedData]) -> None:
    """
    Write pipeline index file.
    :param the_index_file: Full path to index file, including filename.
    :param the_entries: List of StandardizedData objects
    :return: Nothing
    """
    print(f"write_std_pipeline_index the_index_file={the_index_file}", flush=True)
    the_entries.sort(key=lambda my_std_data: (my_std_data.std_archive, my_std_data.version, my_std_data.job_id), reverse=False)
    out_file: io.TextIOWrapper
    with open(the_index_file, 'w', encoding='utf-8') as out_file:
        write_tsv_list(out_file, STD_HEADERS, True, False)
        my_entry: StandardizedData
        for my_entry in the_entries:
            my_entry.write_index_row(out_file)


# ########################################################
# reading files
# ########################################################

def read_std_pipeline_index(the_index_file: str) -> List[StandardizedData]:
    """
    Read index file and return Dictionary of StandardizedData objects
    :param the_index_file: Full path including filename to pipeline index file.
    :return: Dictionary with keys being tuples of std_archive and version, and values being StandardizedData
    """
    print(f"read_std_pipeline_index the_index_file={the_index_file}", flush=True)
    index_file: io.TextIOWrapper
    # std_list use entry.std_archive, entry.version, entry.job_id
    my_entries: List[StandardizedData] = []
    with open(the_index_file, 'r', encoding='utf-8') as index_file:
        keys: List[str] = read_headers(index_file)
        line: str = index_file.readline().rstrip('\n')
        index: int = 1
        while '' != line:
            if 0 == index % 1000:
                print(f"read_std_pipeline_index line={index}", flush=True)
            index += 1
            values: List[str] = line.split("\t")
            tsv_dict: Dict[str, str] = dict(zip(keys, values))
            entry: StandardizedData = StandardizedData.from_dict(tsv_dict)
            my_entries.append(entry)
            line = index_file.readline().rstrip('\n')
    return my_entries


# ########################################################
# archive data
# ########################################################

def get_archive_versions(the_archive_path: str) -> List[str]:
    """
    Get list of versions available for dataset
    :param the_archive_path: Full path, including filename to ZIP file archive.
    :return: List of version timestamp for dataset.
    """
    version_list: List[str] = []
    zip_file: zipfile.ZipFile
    with zipfile.ZipFile(the_archive_path, 'r') as zip_file:
        my_file: str
        for my_file in zip_file.namelist():
            if my_file.startswith("versions/DATA_"):
                if my_file.endswith("/"):
                    # remove first 14 characters, and last one
                    my_file = my_file[14:][:-1]
                    version_list.append(my_file)
    return version_list


# ########################################################
# building index dictionary
# ########################################################


def find_historical_standardized_data(the_std_list: List[StandardizedData], the_version: str,
                                      the_data_archive: str, the_result_archive: str) -> Optional[StandardizedData]:
    """
    Find the historical StandardizedData object that has the version, data, and result paths.
    That way we can get the std_archive
    :param the_std_list: list of StandardizedData objects
    :param the_version: version string
    :param the_data_archive: data archive zip path
    :param the_result_archive: result archive zip path
    :return: StandardizedData object matching version, and data/result archive zip paths
    """
    old_std: Optional[StandardizedData] = None
    my_std: StandardizedData
    for my_std in the_std_list:
        if the_version == my_std.version:
            if the_data_archive == my_std.data_archive:
                if the_result_archive == my_std.result_archive:
                    old_std = my_std
    return old_std


def std_list_contains(the_std_list: List[StandardizedData], the_std_archive: str, the_my_version: str) -> bool:
    """
    Check if archive and version are already in list.
    :param the_std_list: List of StanardizedData (pipline.tsv entries) to search
    :param the_std_archive: archive path to check for
    :param the_my_version: data version to check for
    :return: bool, return True if archive and version are in list
    """
    found: bool = False
    my_std: StandardizedData
    for my_std in the_std_list:
        if my_std.std_archive == the_std_archive:
            if my_std.version == the_my_version:
                found = True
    return found


def build_std_pipeline_index(the_input_dir: str, the_index_file: str, the_write_flag: bool) -> List[StandardizedData]:
    """
    Build and save pipeline index file from current file and available datasets.
    :param the_input_dir: full directory path containing standardized data ZIP files (uses Standardized Data std_archive path)
    :param the_index_file: full path including file name to index file for MBatch index file
    :param the_write_flag: if True, skip writing index. Used for checking new status of pipeline
    :return: Dictionary with keys being tuples of std_archive and version, and values being StandardizedData
    """
    std_list: List[StandardizedData] = []
    # read index file if it exists
    if os.path.exists(the_index_file):
        print(f"read index {the_index_file}", flush=True)
        std_list = read_std_pipeline_index(the_index_file)
    # check for unprocessed standardized data archives, add to std_list, write to pipeline index
    found: bool = False
    dir_path: str
    # dir_names: List[str] replaced with _ to avoid unused-variable warning
    filenames: List[str]
    my_file: str
    # pylint: disable=too-many-nested-blocks
    for dir_path, _, filenames in os.walk(the_input_dir):
        for my_file in filenames:
            if my_file.lower().endswith('.zip'):
                std_archive_full: str = os.path.join(dir_path, my_file)
                std_archive: str = std_archive_full.replace(the_input_dir, '')
                if std_archive.startswith(os.sep):
                    std_archive = std_archive[1:]
                version_list: List[str] = get_archive_versions(std_archive_full)
                if version_list is not None:
                    my_version: str
                    for my_version in version_list:
                        print(f"Check version {my_version} and path {dir_path}", flush=True)
                        if not std_list_contains(std_list, std_archive, my_version):
                            found = True
                            # the_std_archive, the_version, the_data_archive, the_result_archive
                            print(f"New archive {my_version} and path {dir_path}", flush=True)
                            std_list.append(StandardizedData(std_archive, my_version, '', '', '', ''))
    # pylint: enable=too-many-nested-blocks
    if the_write_flag:
        if found:
            print(f"write index {the_index_file}", flush=True)
            write_std_pipeline_index(the_index_file, std_list)
    return std_list
