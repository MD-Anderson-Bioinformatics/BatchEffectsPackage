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

import glob
import os
import shutil
from typing import List
from mbatch.pipeline.std_data import StandardizedData, build_std_pipeline_index
from mbatch.pipeline.std_data import write_std_pipeline_index, read_std_pipeline_index
from mbatch.pipeline.job import create_job, queue_job
from mbatch.test.common import add_error, delete_from_dirs, delete_directory_contents, extract_zip_to_dir
from mbatch.test.test_index import create_index_archive
from mbatch.visualindex.visual_index_base import VisualIndexBase, VisualIndexElementBase
from mbatch.visualindex.visual_index_base import get_dataset_features, get_dataset_samples
from mbatch.visualindex.visual_index_base import get_mut_samples, get_mut_features
from mbatch.visualindex.visual_index_dsc import VisualIndexDsc
from mbatch.visualindex.visual_index_kwd import VisualIndexKwd


# noinspection DuplicatedCode
# pylint: disable=too-many-arguments,too-many-locals,too-many-statements,too-many-branches
def execute_pipeline(the_input_dir: str, the_output_dir: str, the_index_file: str, the_util_dir: str,
                     the_bei_url: str, the_bei_dir: str, the_run_version: str, the_run_source: str,
                     the_base_dir: str, the_sample_column_name: str) -> None:
    """
    Execute the pipeline, running MBatch on any unprocessed datasets
    :param the_input_dir: full directory path containing standardized data ZIP files (uses Standardized Data std_archive path)
    :param the_output_dir: full directory path to which to write MBatch results and data ZIP files
    :param the_index_file: full path including file name to index file for MBatch index file
    :param the_util_dir: full path to util directory
    :param the_bei_url: URL to submit BEI jobs
    :param the_bei_dir: full path to directory containing job directories
    :param the_run_version: timestamp string for MBatch run version
    :param the_run_source: string source, GDC or MWB
    :param the_base_dir: string detailing Docker internal path like /DAPI_MQA/DATA
    :param the_sample_column_name: column string for samples in batch (different for GDC vs MWB)
    :return: nothing
    """
    # read index files
    data_index_file: str = os.path.join(os.path.dirname(the_index_file), f"data_index_{the_run_source.lower()}.tsv")
    dsc_index_file: str = os.path.join(os.path.dirname(the_index_file), f"dsc_index_{the_run_source.lower()}.tsv")
    kwd_index_file: str = os.path.join(os.path.dirname(the_index_file), f"kwd_index_{the_run_source.lower()}.tsv")
    vib: VisualIndexBase = VisualIndexBase(data_index_file, the_base_dir)
    vib.populate_index()
    vid: VisualIndexDsc = VisualIndexDsc(dsc_index_file, the_base_dir)
    vid.populate_index()
    vik: VisualIndexKwd = VisualIndexKwd(kwd_index_file, the_base_dir)
    vik.populate_index()
    # process pipeline
    std_list: List[StandardizedData] = build_std_pipeline_index(the_input_dir, the_index_file, True)
    my_std: StandardizedData
    for my_std in std_list:
        print(f'Standardized Data {my_std.std_archive}', flush=True)
        if my_std.job_id == '':
            if my_std.has_valid_data_p(the_input_dir):
                # Handle Start New Job
                print(f'Start New Job {my_std.version} {my_std.std_archive}', flush=True)
                create_job(my_std, the_bei_url, the_bei_dir, the_input_dir, the_run_version,
                           my_std.has_batch_info_p(the_input_dir), 'TCGA' in my_std.std_archive,
                           the_util_dir, the_run_source, the_sample_column_name, "null", "null", "null")
            else:
                print(f'No-Data for Job {my_std.version} {my_std.std_archive}', flush=True)
                my_std.job_id = "no-processing"
                my_std.job_status = "no-processing"
            # write index file update
            print(f"write updated index {the_index_file}", flush=True)
            write_std_pipeline_index(the_index_file, std_list)
        else:
            if my_std.job_id == 'no-processing':
                print(f'No-Processing for {my_std.version} {my_std.std_archive}', flush=True)
            else:
                # status_job -> created, queued, running, succeeded, failed, no-processing
                # call returns and *updates* instance status
                changed: bool = my_std.status_job(the_bei_url)
                job_status: str = my_std.job_status
                if changed & ('succeeded' == job_status):
                    # Handle succeeded Job
                    print(f'Succeeded {my_std.version} {my_std.std_archive}', flush=True)
                    zip_dir: str = my_std.build_results_path(the_output_dir, the_run_source)
                    job_dir: str = os.path.join(the_bei_dir, my_std.job_id)
                    result_dir: str = os.path.join(job_dir, "ZIP-RESULTS")
                    data_dir: str = os.path.join(job_dir, "ZIP-DATA")
                    results_zip: str
                    data_zip: str
                    # clean up Rplots.pdf files left by other packages
                    delete_from_dirs(job_dir, "Rplots.pdf")
                    results_zip, data_zip = create_index_archive(result_dir, data_dir, zip_dir, os.path.join(result_dir, "info"), my_std, std_list)
                    my_std.result_archive = results_zip
                    my_std.data_archive = data_zip
                    # add dataset information to index and write index
                    data_path: str = os.path.join(job_dir, "ZIP-DATA", "current", "pipeline")
                    vib.find_and_add_entries(the_output_dir, results_zip, data_zip, result_dir, the_run_version, data_path)
                    vib.write_index_file()
                    # add DSC information to index (also for corrections) and write index
                    vid.find_and_add_entries(the_output_dir, results_zip, data_zip, result_dir)
                    vid.write_index_file()
                    # add KWD information to index (also for corrections?) and write index
                    vik.find_and_add_entries(the_output_dir, results_zip, data_zip, result_dir)
                    vik.write_index_file()
                elif changed & ('failed' == job_status):
                    # Handle Failed Job - basically by ignoring
                    print(f'Failure {my_std.version} {my_std.std_archive}', flush=True)
                    add_error(f'Failure Job {my_std.job_id} Version {my_std.version} Archive {my_std.std_archive}')
                elif changed & ('running' == job_status):
                    # Handle Running Job
                    print(f'Running {my_std.version} {my_std.std_archive}', flush=True)
                elif changed & ('queued' == job_status):
                    # Handle queued Job
                    print(f'Queued {my_std.version} {my_std.std_archive}', flush=True)
                if changed:
                    # write index file update
                    print(f"write updated index {the_index_file}", flush=True)
                    write_std_pipeline_index(the_index_file, std_list)
        if my_std.job_status == 'created':
            # make HTTP call to run job and change job_status to running
            queue_job(my_std, the_bei_url)
            # write index file update
            print(f"write updated index {the_index_file}", flush=True)
            write_std_pipeline_index(the_index_file, std_list)
    # corrections done elsewhere, separately
# pylint: enable=too-many-arguments,too-many-locals,too-many-statements,too-many-branches


# noinspection DuplicatedCode
# pylint: disable=too-many-arguments,too-many-locals,too-many-statements
def check_pipeline(the_input_dir: str, the_output_dir: str, the_index_file: str, the_util_dir: str,
                     the_bei_url: str, the_bei_dir: str, the_run_version: str, the_run_source: str,
                     the_base_dir: str, the_sample_column_name: str) -> None:
    """
    Execute the pipeline, running MBatch on any unprocessed datasets
    :param the_input_dir: full directory path containing standardized data ZIP files (uses Standardized Data std_archive path)
    :param the_output_dir: full directory path to which to write MBatch results and data ZIP files
    :param the_index_file: full path including file name to index file for MBatch index file
    :param the_util_dir: full path to util directory
    :param the_bei_url: URL to submit BEI jobs
    :param the_bei_dir: full path to directory containing job directories
    :param the_run_version: timestamp string for MBatch run version
    :param the_run_source: string source, GDC or MWB
    :param the_base_dir: string detailing Docker internal path like /DAPI_MQA/DATA
    :param the_sample_column_name: column string for samples in batch (different for GDC vs MWB)
    :return: nothing
    """
    # read index files
    data_index_file: str = os.path.join(os.path.dirname(the_index_file), f"data_index_{the_run_source.lower()}.tsv")
    dsc_index_file: str = os.path.join(os.path.dirname(the_index_file), f"dsc_index_{the_run_source.lower()}.tsv")
    kwd_index_file: str = os.path.join(os.path.dirname(the_index_file), f"kwd_index_{the_run_source.lower()}.tsv")
    vib: VisualIndexBase = VisualIndexBase(data_index_file, the_base_dir)
    vib.populate_index()
    vid: VisualIndexDsc = VisualIndexDsc(dsc_index_file, the_base_dir)
    vid.populate_index()
    vik: VisualIndexKwd = VisualIndexKwd(kwd_index_file, the_base_dir)
    vik.populate_index()
    # process pipeline
    std_list: List[StandardizedData] = build_std_pipeline_index(the_input_dir, the_index_file, False)
    my_std: StandardizedData
    print('*************************************************', flush=True)
    for my_std in std_list:
        if my_std.job_id == '':
            if my_std.has_valid_data_p(the_input_dir):
                # Handle Start New Job
                print(f'Start New Job {my_std.version} {my_std.std_archive}', flush=True)
            else:
                print(f'No-Data for Job {my_std.version} {my_std.std_archive}', flush=True)
        else:
            if my_std.job_id == 'no-processing':
                print(f'No-Processing for {my_std.version} {my_std.std_archive}', flush=True)
    print('*************************************************', flush=True)
# pylint: enable=too-many-arguments,too-many-locals,too-many-statements


# pylint: disable=too-many-arguments,too-many-locals,too-many-statements
def update_pipeline_data_index(the_input_dir: str, the_output_dir: str, the_index_file: str, the_util_dir: str,
                     the_bei_url: str, the_bei_dir: str, the_run_version: str, the_run_source: str,
                     the_base_dir: str, the_sample_column_name: str) -> None:
    """
    Execute the pipeline, running MBatch on any unprocessed datasets
    :param the_input_dir: full directory path containing standardized data ZIP files (uses Standardized Data std_archive path)
    :param the_output_dir: full directory path to which to write MBatch results and data ZIP files
    :param the_index_file: full path including file name to index file for MBatch index file
    :param the_util_dir: full path to util directory
    :param the_bei_url: URL to submit BEI jobs
    :param the_bei_dir: full path to directory containing job directories
    :param the_run_version: timestamp string for MBatch run version
    :param the_run_source: string source, GDC or MWB
    :param the_base_dir: string detailing Docker internal path like /DAPI_MQA/DATA
    :param the_sample_column_name: column string for samples in batch (different for GDC vs MWB)
    :return: nothing
    """
    # read pipeline index file
    std_list: List[StandardizedData] = read_std_pipeline_index(the_index_file)
    # read data/advanced index file
    data_index_file: str = os.path.join(os.path.dirname(the_index_file), f"data_index_{the_run_source.lower()}.tsv")
    vib: VisualIndexBase = VisualIndexBase(data_index_file, the_base_dir)
    vib.populate_index()
    data_entry: StandardizedData
    for data_entry in std_list:
        if "succeeded" == data_entry.job_status:
            element: VisualIndexElementBase
            for element in vib.m_ele_dict.values():
                zip_path: str = element.m_path_data
                # replace the_base_dir with the_output_dir
                zip_path = zip_path.replace(the_base_dir, the_output_dir)
                if zip_path == data_entry.data_archive:
                    print(f'zip_path = {zip_path}', flush=True)
                    # find job output directory
                    file_dir_path: str = os.path.join(the_bei_dir, data_entry.job_id, "ZIP-DATA", "original")
                    # count number of samples
                    samples_matrix: int = get_dataset_samples(file_dir_path, 'matrix_data.tsv')
                    print(f'samples_matrix = {samples_matrix}', flush=True)
                    samples_mutations: int = get_mut_samples(file_dir_path, 'mutations.tsv')
                    print(f'samples_mutations = {samples_mutations}', flush=True)
                    # count number of features
                    features_matrix: int = get_dataset_features(file_dir_path, 'matrix_data.tsv')
                    print(f'features_matrix = {features_matrix}', flush=True)
                    features_mutations: int = get_mut_features(file_dir_path, 'mutations.tsv')
                    print(f'features_mutations = {features_mutations}', flush=True)
                    element.samples_matrix = samples_matrix
                    element.samples_mutations = samples_mutations
                    element.features_matrix = features_matrix
                    element.features_mutations = features_mutations
    vib.write_index_file()
# pylint: enable=too-many-arguments,too-many-locals,too-many-statements


# pylint: disable=too-many-arguments,too-many-locals,too-many-statements
def update_pipeline_rebuild_results_archives(the_input_dir: str, the_output_dir: str, the_index_file: str, the_util_dir: str,
                     the_bei_url: str, the_bei_dir: str, the_run_version: str, the_run_source: str,
                     the_base_dir: str, the_sample_column_name: str) -> None:
    """
    Execute the pipeline, running MBatch on any unprocessed datasets
    :param the_input_dir: full directory path containing standardized data ZIP files (uses Standardized Data std_archive path)
    :param the_output_dir: full directory path to which to write MBatch results and data ZIP files
    :param the_index_file: full path including file name to index file for MBatch index file
    :param the_util_dir: full path to util directory
    :param the_bei_url: URL to submit BEI jobs
    :param the_bei_dir: full path to directory containing job directories
    :param the_run_version: timestamp string for MBatch run version
    :param the_run_source: string source, GDC or MWB
    :param the_base_dir: string detailing Docker internal path like /DAPI_MQA/DATA
    :param the_sample_column_name: column string for samples in batch (different for GDC vs MWB)
    :return: nothing
    """
    # read pipeline index file
    std_list: List[StandardizedData] = read_std_pipeline_index(the_index_file)
    my_std: StandardizedData
    for my_std in std_list:
        print(f'Standardized Data {my_std.std_archive}', flush=True)
        if "succeeded" == my_std.job_status:
            # rebuild results ZIP
            print(f'my_std.std_archive = {my_std.std_archive}', flush=True)
            zip_dir: str = my_std.build_results_path(the_output_dir, the_run_source)
            job_dir: str = os.path.join(the_bei_dir, my_std.job_id)
            result_dir: str = os.path.join(job_dir, "ZIP-RESULTS")
            data_dir: str = os.path.join(job_dir, "ZIP-DATA")
            results_zip: str
            data_zip: str
            # should already have old data and results added, but double check
            results_zip, data_zip = create_index_archive(result_dir, data_dir, zip_dir, os.path.join(result_dir, "info"), my_std, std_list)
            print(f'results_zip = {results_zip}', flush=True)
            print(f'data_zip = {data_zip}', flush=True)
# pylint: enable=too-many-arguments,too-many-locals,too-many-statements


# noinspection DuplicatedCode
# pylint: disable=too-many-arguments,too-many-locals,too-many-statements
def reindex_pipeline(the_dir_result_data: str, the_temp_dir: str, the_old_tag: str) -> None:
    """
    Reindex pipeline -
    :param the_dir_result_data: full directory path to top level with results/data zips
    :param the_temp_dir: full directory path to temporary directory for working with zips
    :param the_old_tag: tag to add to end of old ZIP file
    :return: nothing
    """
    # find files that end in -results.zip
    start_dir: str = os.path.join(the_dir_result_data, "**", "*-results.zip")
    result_zips: List[str] = glob.glob(start_dir, recursive=True)
    result_zips.sort()
    results_filename: str
    print('*************************************************', flush=True)
    count: int = 0
    max_count: int = 0
    for results_filename in result_zips:
        print(f'{count} -- {results_filename}', flush=True)
        count += 1
    max_count = count
    count = 0
    print('*************************************************', flush=True)
    for results_filename in result_zips:
        print('----------------------------------------------------------------------', flush=True)
        print(f'{count}/{max_count} -- {results_filename}', flush=True)
        data_filename: str = results_filename.replace('-results.', '-data.')
        results_dir: str = os.path.join(the_temp_dir, "ZIP-RESULTS")
        data_dir: str = os.path.join(the_temp_dir, "ZIP-DATA")
        print(f'results_filename={results_filename}', flush=True)
        print(f'data_filename={data_filename}', flush=True)
        print(f'results_dir={results_dir}', flush=True)
        print(f'data_dir={data_dir}', flush=True)
        # clean contents of the_temp_dir
        delete_directory_contents(the_temp_dir)
        os.mkdir(results_dir)
        os.mkdir(data_dir)
        # extract results_filename into the_temp_dir/ZIP-RESULTS
        extract_zip_to_dir(results_filename, results_dir)
        # extract data_filename into the_temp_dir/ZIP-DATA
        extract_zip_to_dir(data_filename, data_dir)
        # rebuild index file
        new_results_zip: str
        new_data_zip: str
        zip_dir: str = the_temp_dir
        new_results_zip, new_data_zip = create_index_archive(results_dir, data_dir, zip_dir, os.path.join(results_dir, "info"), None, None)
        print(f'new_results_zip={new_results_zip}', flush=True)
        print(f'new_data_zip={new_data_zip}', flush=True)
        # rename old ZIP files
        rename_results: str = results_filename + the_old_tag
        rename_data: str = data_filename + the_old_tag
        print(f'rename_results={rename_results}', flush=True)
        os.rename(results_filename, rename_results)
        print(f'rename_data={rename_data}', flush=True)
        os.rename(data_filename, rename_data)
        # copy new ZIP files into place
        print(f'{new_results_zip} to {results_filename}', flush=True)
        shutil.copyfile(new_results_zip, results_filename)
        print(f'{new_data_zip} to {data_filename}', flush=True)
        shutil.copyfile(new_data_zip, data_filename)
        count += 1
    print('*************************************************', flush=True)
# pylint: enable=too-many-arguments,too-many-locals,too-many-statements
