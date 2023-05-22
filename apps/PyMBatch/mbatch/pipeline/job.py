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

import os
import io
import zipfile
import json
from typing import List
import jsonpickle
import pandas
from mbatch.index.index_original_data import object_decoder_from_convert, OriginalData
from mbatch.pipeline.std_data import StandardizedData
from mbatch.gdcapi.standardized_data import write_converted_dataframe
from mbatch.test.common import delete_directory_contents, copy_archive_file_to_regular, get_current_timestamp


def copy_json_original_files(the_zip_file: str, the_source: str, the_version: str, the_out_dir: str) -> None:
    """
    Copy files to original dir
    :param the_zip_file: ZIP file to copy from
    :param the_source: source string
    :param the_version: version string
    :param the_out_dir: output directory
    :return:
    """
    # copy index.json from ZIP file to original_data.json
    zip_file: zipfile.ZipFile
    read_zip_file: io.TextIOWrapper
    with zipfile.ZipFile(the_zip_file, 'r') as zip_file:
        with zip_file.open("index.json", mode="r") as read_zip_file:
            # read file (does not have source attribute)
            json_data: str = read_zip_file.read()
            # make object
            my_obj: OriginalData = json.loads(json_data, object_hook=object_decoder_from_convert)
            # set source and version attributes
            my_obj.source = the_source
            my_obj.version = the_version
            # no b in open call, since not copying from ZIP (which requires binary)
            # writing a string/json which is just 'w'
            out_file: io.TextIOWrapper
            with open(os.path.join(the_out_dir, "original_data.json"), 'w', encoding='utf-8') as out_file:
                # write out file
                out_file.write(jsonpickle.encode(my_obj, indent=4, unpicklable=False))
                # json.dump(jsonpickle.encode(my_obj, indent=4, unpicklable=False), out_file)
                # json.dump(json.dumps(my_obj), out_file)


def write_zip_files(the_std_data: StandardizedData, the_job_zip_dir: str, the_run_version: str,
                    the_input_dir: str, the_source: str) -> None:
    """
    Create normal files for data and job information in job directories to be ZIPped
    :param the_std_data: object representing a standardized data set
    :param the_job_zip_dir: directory in job which will be ZIPped into result or data ZIP
    :param the_run_version: timestamp string for MBatch run version
    :param the_input_dir: full directory path containing standardized data ZIP files (uses Standardized Data std_archive path)
    :param the_source: string giving source (GDC or MWB)
    :return: Nothing
    """
    print(f'write_zip_files the_job_zip_dir={the_job_zip_dir}', flush=True)
    print(f'write_zip_files the_run_version={the_run_version}', flush=True)
    print(f'write_zip_files the_input_dir={the_input_dir}', flush=True)
    zip_arch: str = os.path.join(the_input_dir, the_std_data.std_archive)
    print(f'write_zip_files zip_arch={zip_arch}', flush=True)
    out_file: io.TextIOWrapper
    with open(os.path.join(the_job_zip_dir, "source_id.txt"), 'w', encoding='utf-8') as out_file:
        out_file.write(f"{the_std_data.calc_md5()}\n")
    with open(os.path.join(the_job_zip_dir, "version_stamp.txt"), 'w', encoding='utf-8') as out_file:
        out_file.write(f"{the_run_version}\n")
    # TODO: update to handle Original and Corrected*
    with open(os.path.join(the_job_zip_dir, "version_type.txt"), 'w', encoding='utf-8') as out_file:
        out_file.write("Original\n")
    copy_json_original_files(zip_arch, the_source, the_std_data.version, the_job_zip_dir)


def check_batch_type(the_batch_info: pandas.DataFrame, the_sample_count: int, the_greater_than: float, the_less_than: float) -> List[str]:
    """
    Check the batch information using rule-of-thumb guidelines to pick "interesting" batch types.
    Calling function will reduce level of interest until at least one batch type is selected.
    :param the_batch_info: DataFrame of batch information. One column is samples (aliquot_barcode) rest are potential batch types
    :param the_sample_count: number of samples in this data
    :param the_greater_than: fraction to use in comparison for interesting batch type: number of batches is greater than (.9 * number of samples)
    :param the_less_than: fraction to use in comparison for interesting batch type: number of batches is less than (.6 * number of samples)
    :return: list of sample column (aliquot_barcode) and interesting batch types
    """
    valid_columns: List[str] = []
    batch_type: str
    ignore_columns: List[str] = ['aliquot_uuid', 'patient_barcode', 'patient_uuid']
    # pylint: disable=too-many-nested-blocks
    for batch_type in the_batch_info:
        if 'aliquot_barcode' == batch_type:
            valid_columns.append(batch_type)
        elif 'Sample' == batch_type:
            valid_columns.append(batch_type)
        elif batch_type not in ignore_columns:
            my_series: pandas.Series = the_batch_info[batch_type]
            # Evaluate Batch Types Based On Number of Batches
            # Remove Non-Batches
            # 1. Get total number of batches.
            # 2. Remove counts for any batches that are "-" or "".
            # If one or zero Batches are left, do not use this Batch Type
            uniq_vals: pandas.Series = my_series.value_counts()
            batch_name: str
            batch_count: int = 0
            for batch_name in uniq_vals.keys():
                if ('-' != batch_name) & ('' != batch_name):
                    batch_count += 1
            if batch_count > 1:
                print(f"count is {batch_count} and sample count is {the_sample_count} batch type is {batch_type}", flush=True)
                # IF number of batches is less than number of samples
                # AND the number of batches is less than (.9 * number of samples)
                # THEN use this Batch Type (for now)
                if (batch_count < the_sample_count) & (batch_count < (the_greater_than * the_sample_count)):
                    # Evaluate Batch Types Based on Samples with Batch Values
                    # Count number of samples with Batch values that are not "-" or "".
                    # IF the number of batches is greater than (.6 * number of samples)
                    # THEN do use this Batch Type
                    if batch_count > (the_less_than * the_sample_count):
                        valid_columns.append(batch_type)
    # pylint: enable=too-many-nested-blocks
    return valid_columns


def easier_check_batch_type(the_batch_info: pandas.DataFrame, the_min_count: int, the_max_count: int) -> List[str]:
    """
    Check the batch information using rule-of-thumb guidelines to pick "interesting" batch types.
    Calling function will reduce level of interest until at least one batch type is selected.
    :param the_batch_info: DataFrame of batch information. One column is samples (aliquot_barcode) rest are potential batch types
    :param the_min_count: minimum number of batches to accept
    :param the_max_count: maximum number of batches to accept
    :return: list of sample column (aliquot_barcode) and interesting batch types
    """
    # TODO: test this with metabolomics data
    valid_columns: List[str] = []
    batch_type: str
    ignore_columns: List[str] = ['aliquot_uuid', 'patient_barcode', 'patient_uuid']
    # pylint: disable=too-many-nested-blocks
    for batch_type in the_batch_info:
        if 'aliquot_barcode' == batch_type:
            valid_columns.append(batch_type)
        elif 'Sample' == batch_type:
            valid_columns.append(batch_type)
        elif batch_type not in ignore_columns:
            my_series: pandas.Series = the_batch_info[batch_type]
            # Evaluate Batch Types Based On Number of Batches
            # Remove Non-Batches
            # 1. Get total number of batches.
            # 2. Remove counts for any batches that are "-" or "".
            # If one or zero Batches are left, do not use this Batch Type
            uniq_vals: pandas.Series = my_series.value_counts()
            batch_name: str
            batch_count: int = 0
            for batch_name in uniq_vals.keys():
                if ('-' != batch_name) & ('' != batch_name):
                    batch_count += 1
            if batch_count > 1:
                if batch_count < the_max_count:
                    if batch_count > the_min_count:
                        print(f"count is {batch_count} batch type is {batch_type}", flush=True)
                        valid_columns.append(batch_type)
    # pylint: enable=too-many-nested-blocks
    return valid_columns


# pylint: disable=too-many-nested-blocks
def filter_usable_batch_types(the_batch_info: pandas.DataFrame, the_tcga_flag: bool) -> pandas.DataFrame:
    """
    Filter the batch data for interesting batch types
    :param the_batch_info: DataFrame of batch information. One column is samples (aliquot_barcode) rest are potential batch types
    :param the_tcga_flag: if True, then use normal batch types for TCGA data
    :return: DataFrame of batch information with columns reduced to samples (aliquot_barcode) and interesting batches
    """
    if the_tcga_flag:
        valid_columns: List[str] = ['aliquot_barcode', 'Sample', 'ship_date', 'tissue_source_site', 'batch_id', 'source_center', 'sample_type_name']
        the_batch_info = the_batch_info[valid_columns]
    else:
        # keep aliquot_barcode
        greater_than: float = 0.9
        less_than: float = 0.5
        valid_columns: List[str] = check_batch_type(the_batch_info, the_batch_info.shape[0], greater_than, less_than)
        # reduce criteria to find some sort of batch type if available
        max_count: int = the_batch_info.shape[0] - 1
        min_count = the_batch_info.shape[0] - 1
        while (len(valid_columns) < 2) & (min_count > 0):
            min_count -= 1
            valid_columns = easier_check_batch_type(the_batch_info, min_count, max_count)
        the_batch_info = the_batch_info[valid_columns]
    return the_batch_info
# pylint: enable=too-many-nested-blocks


# pylint: disable=too-many-arguments,too-many-locals,too-many-statements
def create_job(the_std_data: StandardizedData, the_bei_url: str, the_bei_dir: str, the_input_dir: str,
               the_run_version: str, the_has_batch_info_flag: bool, the_tcga_flag: bool,
               the_util_dir: str, the_run_source: str, the_sample_column_name: str) -> None:
    """
    Create the job and populate its job directories.
    :param the_std_data: object representing a standardized data set
    :param the_bei_url: URL to submit BEI jobs
    :param the_bei_dir: full path to directory containing job directories
    :param the_input_dir: full directory path containing standardized data ZIP files (uses Standardized Data std_archive path)
    :param the_run_version: timestamp string for MBatch run version
    :param the_has_batch_info_flag: is True if data has batch information
    :param the_tcga_flag: is True if data is TCGA data
    :param the_util_dir: full directory path to util dir
    :param the_run_source: string source, GDC or MWB
    :param the_sample_column_name: column string for samples in batch (different for GDC vs MWB)
    :return: Nothing
    """
    print(f'create_job start {the_std_data.version} {the_std_data.std_archive}', flush=True)
    # also sets the job_id attribute of the_std_data
    job_id: str = the_std_data.create_new_job(the_bei_url)
    job_dir: str = os.path.join(the_bei_dir, job_id)
    # clear job directory
    delete_directory_contents(job_dir)
    # make Results and Data directories
    result_dir: str = os.path.join(job_dir, "ZIP-RESULTS")
    os.makedirs(result_dir, exist_ok=True)
    data_dir: str = os.path.join(job_dir, "ZIP-DATA")
    os.makedirs(data_dir, exist_ok=True)
    ########################################################
    # Populate Results directory
    write_zip_files(the_std_data, result_dir, the_run_version, the_input_dir, the_run_source)
    ########################################################
    # create and populate Data/original directory
    original_dir: str = os.path.join(data_dir, "original")
    os.makedirs(original_dir, exist_ok=True)
    # Populate Data directory
    write_zip_files(the_std_data, data_dir, the_run_version, the_input_dir, the_run_source)
    # copy changeable files to Data/original directory
    # only copies if they exist
    sample_id_batch_type: str = the_sample_column_name
    batch_types: List[str]
    if the_has_batch_info_flag:
        batch_info: pandas.DataFrame = the_std_data.read_batch_info(the_input_dir)
        # reduce batches to usable batch information
        batch_info = filter_usable_batch_types(batch_info, the_tcga_flag)
        if len(batch_info.columns) > 1:
            write_converted_dataframe(batch_info, original_dir, "batches.tsv", the_sample_column_name)
            # batch_types
            batch_types = list(batch_info.columns)
            batch_types.remove(sample_id_batch_type)
        else:
            the_std_data.make_fake_batch_info(the_input_dir, os.path.join(original_dir, "batches.tsv"), the_sample_column_name)
            batch_types = ['example']
    else:
        # create fake batch info
        # aliquot_barcode example
        the_std_data.make_fake_batch_info(the_input_dir, os.path.join(original_dir, "batches.tsv"), the_sample_column_name)
        batch_types = ['example']
    copy_archive_file_to_regular(os.path.join(the_input_dir, the_std_data.std_archive), f"versions/DATA_{the_std_data.version}/matrix.tsv", os.path.join(original_dir, "matrix_data.tsv"))
    copy_archive_file_to_regular(os.path.join(the_input_dir, the_std_data.std_archive), f"versions/DATA_{the_std_data.version}/clinical.tsv", os.path.join(original_dir, "clinical.tsv"))
    copy_archive_file_to_regular(os.path.join(the_input_dir, the_std_data.std_archive), f"versions/DATA_{the_std_data.version}/mutations.tsv", os.path.join(original_dir, "mutations.tsv"))
    # below two files only exist for Metabolomics Workbench Data
    copy_archive_file_to_regular(os.path.join(the_input_dir, the_std_data.std_archive), f"versions/DATA_{the_std_data.version}/ngchm_link_map.tsv", os.path.join(original_dir, "ngchm_link_map.tsv"))
    copy_archive_file_to_regular(os.path.join(the_input_dir, the_std_data.std_archive), f"versions/DATA_{the_std_data.version}/row_col_types.tsv", os.path.join(original_dir, "row_col_types.tsv"))
    # copy historical static version of files to Data/original directory
    # only copies if they exist
    if the_has_batch_info_flag:
        copy_archive_file_to_regular(os.path.join(the_input_dir, the_std_data.std_archive), f"versions/DATA_{the_std_data.version}/batches.tsv", os.path.join(original_dir, "original_batches.tsv"))
    else:
        # create fake batch info
        the_std_data.make_fake_batch_info(the_input_dir, os.path.join(original_dir, "original_batches.tsv"), the_sample_column_name)
    copy_archive_file_to_regular(os.path.join(the_input_dir, the_std_data.std_archive), f"versions/DATA_{the_std_data.version}/matrix.tsv", os.path.join(original_dir, "original_matrix_data.tsv"))
    copy_archive_file_to_regular(os.path.join(the_input_dir, the_std_data.std_archive), f"versions/DATA_{the_std_data.version}/clinical.tsv", os.path.join(original_dir, "original_clinical.tsv"))
    copy_archive_file_to_regular(os.path.join(the_input_dir, the_std_data.std_archive), f"versions/DATA_{the_std_data.version}/mutations.tsv", os.path.join(original_dir, "original_mutations.tsv"))
    # below two files only exist for Metabolomics Workbench Data
    copy_archive_file_to_regular(os.path.join(the_input_dir, the_std_data.std_archive), f"versions/DATA_{the_std_data.version}/ngchm_link_map.tsv", os.path.join(original_dir, "original_ngchm_link_map.tsv"))
    copy_archive_file_to_regular(os.path.join(the_input_dir, the_std_data.std_archive), f"versions/DATA_{the_std_data.version}/row_col_types.tsv", os.path.join(original_dir, "original_row_col_types.tsv"))
    # setup config files
    # default template lives in BatchEffects_clean/BatchEffectsPackage/data/testing_static/PyMBatch/mbatch
    # but is in util directory during pipeline run
    # batch_id,plate_id,sample_type_name,sex,ship_date,source_center,tissue_source_site,sample_type_name
    # TODO: copy/update MBatchConfig file
    config_template: str = os.path.join(the_util_dir, "MBatchConfig_template.tsv")
    config_file: io.TextIOWrapper
    file_data: str
    with open(config_template, 'r', encoding='utf-8') as config_file:
        file_data = config_file.read()
    file_data = file_data.replace('<jobId>', the_std_data.job_id)
    file_data = file_data.replace('<configDesc>', f"MBatch run on {get_current_timestamp()}")
    file_data = file_data.replace('<sampleIdBatchType>', sample_id_batch_type)
    file_data = file_data.replace('<batchTypes>', ",".join(batch_types))
    file_data = file_data.replace('<filterLogTransformFlag>', str(the_std_data.needs_log_transform(the_run_source)).upper())
    file_data = file_data.replace('<title>', the_std_data.get_title(the_run_version))
    file_data = file_data.replace('<DataVersion>', f"DATA_{the_std_data.version}")
    file_data = file_data.replace('<TestVersion>', f"TEST_{the_run_version}")
    file_data = file_data.replace('<replaceNAs>', f"{the_std_data.needs_replace_na(the_run_source)}".upper())
    config_file_path: str = os.path.join(result_dir, "MBatchConfig.tsv")
    with open(config_file_path, 'w', encoding='utf-8') as config_file:
        config_file.write(file_data)
    print(f'create_job done {job_id} {the_std_data.version} {the_std_data.std_archive}', flush=True)
# pylint: enable=too-many-arguments,too-many-locals,too-many-statements


def queue_job(the_std_data: StandardizedData, the_bei_url: str) -> None:
    """
    Tell BEI to start the job when a compute image is available
    :param the_std_data: object representing a standardized data set
    :param the_bei_url: URL for BEI, used to trigger the run of the dataset
    :return: Nothing
    """
    # make HTTP call to run job and change job_status to queued
    the_std_data.queue_job(the_bei_url)
