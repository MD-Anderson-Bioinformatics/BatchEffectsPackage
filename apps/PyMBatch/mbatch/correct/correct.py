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
from typing import Dict, List
from mbatch.pipeline.std_data import StandardizedData, build_std_pipeline_index
from mbatch.pipeline.std_data import find_historical_standardized_data, write_std_pipeline_index
from mbatch.visualindex.visual_index_base import VisualIndexBase, VisualIndexElementBase
from mbatch.visualindex.visual_index_dsc import VisualIndexDsc, VisualIndexElementDsc
from mbatch.visualindex.visual_index_kwd import VisualIndexKwd, VisualIndexElementKwd
from mbatch.pipeline.job import create_job, queue_job
from mbatch.test.common import add_error, delete_from_dirs
from mbatch.test.test_index import create_index_archive


def extract_batch(the_analysis_path: str) -> str:
    """
    Get the batch type from the analysis path.
    Analysis path for PCA+ is /dir/batchtype/dir/dir, so split and get batch type, which is [2].
    :param the_analysis_path: path to split
    :return: batch type string to return
    """
    splitted: List[str] = the_analysis_path.split("/")
    return splitted[2]


# noinspection DuplicatedCode
# pylint: disable=too-many-arguments,too-many-locals,too-many-statements,too-many-nested-blocks,too-many-branches
def continue_correct(the_input_dir: str, the_output_dir: str, the_index_file: str, the_util_dir: str,
                     the_bei_url: str, the_bei_dir: str, the_run_version: str, the_run_source: str,
                     the_base_dir: str, the_sample_column_name: str) -> None:
    """
    Check the pipeline for correctable datasets, running MBatch on any unprocessed datasets
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
    # TODO: code is partially duplicated in check_correct
    # TODO: code is partially duplicated in pipeline.py
    # read index files
    data_index_file: str = os.path.join(os.path.dirname(the_index_file), f"data_index_{the_run_source.lower()}.tsv")
    dsc_index_file: str = os.path.join(os.path.dirname(the_index_file), f"dsc_index_{the_run_source.lower()}.tsv")
    # kwd_index_file: str = os.path.join(os.path.dirname(the_index_file), f"kwd_index_{the_run_source.lower()}.tsv")
    vib: VisualIndexBase = VisualIndexBase(data_index_file, the_base_dir)
    vib.populate_index()
    vid: VisualIndexDsc = VisualIndexDsc(dsc_index_file, the_base_dir)
    vid.populate_index()
    # vik: VisualIndexKwd = VisualIndexKwd(kwd_index_file, the_base_dir)
    # vik.populate_index()
    # process pipeline
    std_list: List[StandardizedData] = build_std_pipeline_index(the_input_dir, the_index_file, False)
    print('*************************************************', flush=True)
    my_std: StandardizedData
    for my_std in std_list:
        if my_std.job_id == '':
            add_error(f'(continue_correct) Should not have StandardizedData without job id here {my_std.version} {my_std.std_archive}')
        else:
            if my_std.job_id != 'no-processing':
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
                    # vik.find_and_add_entries(the_output_dir, results_zip, data_zip, result_dir)
                    # vik.write_index_file()
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
    print('*************************************************', flush=True)
# pylint: enable=too-many-arguments,too-many-locals,too-many-statements,too-many-nested-blocks,too-many-branches


# noinspection DuplicatedCode
# pylint: disable=too-many-arguments,too-many-locals,too-many-statements,too-many-nested-blocks,too-many-branches
def initiate_correct(the_input_dir: str, the_output_dir: str, the_index_file: str, the_util_dir: str,
                     the_bei_url: str, the_bei_dir: str, the_run_version: str, the_run_source: str,
                     the_base_dir: str, the_sample_column_name: str) -> None:
    """
    Check the pipeline for correctable datasets, running MBatch on any unprocessed datasets
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
    # TODO: code is partially duplicated in check_correct
    # TODO: code is partially duplicated in pipeline.py
    # read index files
    data_index_file: str = os.path.join(os.path.dirname(the_index_file), f"data_index_{the_run_source.lower()}.tsv")
    dsc_index_file: str = os.path.join(os.path.dirname(the_index_file), f"dsc_index_{the_run_source.lower()}.tsv")
    # kwd_index_file: str = os.path.join(os.path.dirname(the_index_file), f"kwd_index_{the_run_source.lower()}.tsv")
    vib: VisualIndexBase = VisualIndexBase(data_index_file, the_base_dir)
    vib.populate_index()
    vid: VisualIndexDsc = VisualIndexDsc(dsc_index_file, the_base_dir)
    vid.populate_index()
    # vik: VisualIndexKwd = VisualIndexKwd(kwd_index_file, the_base_dir)
    # vik.populate_index()
    # process pipeline
    std_list: List[StandardizedData] = build_std_pipeline_index(the_input_dir, the_index_file, False)
    ##########################################################################
    # collect DSC correctable datasets
    ##########################################################################
    # new corrections to start
    new_correction_jobs_dsc: Dict[str, VisualIndexElementDsc] = {}
    print('*************************************************', flush=True)
    my_dsc: VisualIndexElementDsc
    for my_dsc in vid.m_ele_dict.values():
        if float(my_dsc.m_overall_dsc) > 0.30:
            if '<0.0005' == my_dsc.m_overall_dsc_pvalue:
                # check if already corrected
                corrected: bool = False
                my_data_index_entry: VisualIndexElementBase
                for my_data_index_entry in vib.m_ele_dict.values():
                    if my_data_index_entry.m_path_results == my_dsc.m_path_results:
                        if my_data_index_entry.m_data_version == my_dsc.m_data_version:
                            if "Adjusted-EB_withPara-DSC" == my_data_index_entry.job_type:
                                corrected = True
                if not corrected:
                    # add to list after filtering for largest DSC
                    if my_dsc.m_id not in new_correction_jobs_dsc:
                        new_correction_jobs_dsc[my_dsc.m_id] = my_dsc
                    else:
                        old_dsc: VisualIndexElementDsc = new_correction_jobs_dsc[my_dsc.m_id]
                        if float(my_dsc.m_overall_dsc) > float(old_dsc.m_overall_dsc):
                            new_correction_jobs_dsc[my_dsc.m_id] = my_dsc
    ##########################################################################
    # collect KWD correctable datasets
    ##########################################################################
    # new corrections to start
    # new_correction_jobs_kwd: Dict[str, VisualIndexElementKwd] = {}
    # print('*************************************************', flush=True)
    # my_kwd: VisualIndexElementKwd
    # for my_kwd in vik.m_ele_dict.values():
    #     if float(my_kwd.m_neg_log10_pvalue) >= float(my_kwd.m_neg_log10_cutoff):
    #         # second elements of analysis path is batch type
    #         # my_kwd.m_batches_called is (batch-id-1, batch-id-2)
    #         # check if already corrected
    #         corrected: bool = False
    #         my_data_index_entry: VisualIndexElementBase
    #         for my_data_index_entry in vib.m_ele_dict.values():
    #             if my_data_index_entry.m_path_results == my_kwd.m_path_results:
    #                 if my_data_index_entry.m_data_version == my_kwd.m_data_version:
    #                     if "Adjusted-EB_withPara-KWD" == my_data_index_entry.job_type:
    #                         corrected = True
    #         if not corrected:
    #             # add to list after filtering for largest DSC
    #             if my_kwd.m_id not in new_correction_jobs_kwd:
    #                 new_correction_jobs_kwd[my_kwd.m_id] = my_kwd
    #             else:
    #                 old_kwd: VisualIndexElementKwd = new_correction_jobs_kwd[my_kwd.m_id]
    #                 if float(my_kwd.m_neg_log10_pvalue) > float(old_kwd.m_neg_log10_pvalue):
    #                     new_correction_jobs_kwd[my_kwd.m_id] = my_kwd
    ##########################################################################
    # run DSC corrections
    ##########################################################################
    for my_dsc in new_correction_jobs_dsc.values():
        old_std: StandardizedData = find_historical_standardized_data(std_list, my_dsc.m_data_version, my_dsc.m_path_data, my_dsc.m_path_results)
        if old_std is not None:
            my_std: StandardizedData = StandardizedData(old_std.std_archive, my_dsc.m_data_version, '', '', '', '')
            # Handle Start New CORRECTION Job
            print(f'Start New CORRECTION Job {round(float(my_dsc.m_overall_dsc),2)} {my_dsc.m_overall_dsc_pvalue} {my_dsc.m_data_version} {my_dsc.m_path_data} {extract_batch(my_dsc.m_analysis_path)}', flush=True)
            create_job(my_std, the_bei_url, the_bei_dir, the_input_dir,
                       the_run_version, my_std.has_batch_info_p(the_input_dir), 'TCGA' in my_std.std_archive,
                       the_util_dir, the_run_source, the_sample_column_name,
                       extract_batch(my_dsc.m_analysis_path), "EB_withPara", "DSC")
            # write pipeline index std_list
            if my_std.job_status == 'created':
                # NOTE: std_list cannot directly track back to vib unfortunately
                # make HTTP call to run job and change job_status to running
                queue_job(my_std, the_bei_url)
                # add my_std entry to std_list to track progress
                std_list.append(my_std)
                # write index file update
                print(f"write updated index {the_index_file}", flush=True)
                write_std_pipeline_index(the_index_file, std_list)
    ##########################################################################
    # run KWD corrections
    ##########################################################################
    # for my_kwd in new_correction_jobs_kwd.values():
    #     old_std: StandardizedData = find_historical_standardized_data(std_list, my_kwd.m_data_version, my_kwd.m_path_data, my_kwd.m_path_results)
    #     if old_std is not None:
    #         my_std: StandardizedData = StandardizedData(old_std.std_archive, my_kwd.m_data_version, '', '', '', '')
    #         # Handle Start New CORRECTION Job
    #         print(f'Start New CORRECTION Job {round(float(my_kwd.m_neg_log10_pvalue), 2)} {my_kwd.m_neg_log10_cutoff} {my_kwd.m_data_version} {my_kwd.m_path_data} {extract_batch(my_kwd.m_analysis_path)} {my_kwd.m_batches_called}', flush=True)
    #         # handle KWD vs DSC specific run versions
    #         create_job(my_std, the_bei_url, the_bei_dir, the_input_dir,
    #                    f'{the_run_version}-KWD', my_std.has_batch_info_p(the_input_dir), 'TCGA' in my_std.std_archive,
    #                    the_util_dir, the_run_source, the_sample_column_name,
    #                    extract_batch(my_kwd.m_analysis_path), "EB_withPara", "KWD")
    #         # write pipeline index std_list
    #         if my_std.job_status == 'created':
    #             # NOTE: std_list cannot directly track back to vib unfortunately
    #             # make HTTP call to run job and change job_status to running
    #             queue_job(my_std, the_bei_url)
    #             # add my_std entry to std_list to track progress
    #             std_list.append(my_std)
    #             # write index file update
    #             print(f"write updated index {the_index_file}", flush=True)
    #             write_std_pipeline_index(the_index_file, std_list)
    print('*************************************************', flush=True)
# pylint: enable=too-many-arguments,too-many-locals,too-many-statements,too-many-nested-blocks,too-many-branches


# noinspection DuplicatedCode
# pylint: disable=too-many-arguments,too-many-locals,too-many-statements,too-many-nested-blocks,too-many-branches
def check_correct(the_input_dir: str, the_output_dir: str, the_index_file: str, the_util_dir: str,
                  the_bei_url: str, the_bei_dir: str, the_run_version: str, the_run_source: str,
                  the_base_dir: str, the_sample_column_name: str) -> None:
    """
    Check the pipeline for correctable datasets, running MBatch on any unprocessed datasets
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
    # TODO: code is partially duplicated in execute_correct
    # read index files
    data_index_file: str = os.path.join(os.path.dirname(the_index_file), f"data_index_{the_run_source.lower()}.tsv")
    dsc_index_file: str = os.path.join(os.path.dirname(the_index_file), f"dsc_index_{the_run_source.lower()}.tsv")
    # kwd_index_file: str = os.path.join(os.path.dirname(the_index_file), f"kwd_index_{the_run_source.lower()}.tsv")
    vib: VisualIndexBase = VisualIndexBase(data_index_file, the_base_dir)
    vib.populate_index()
    vid: VisualIndexDsc = VisualIndexDsc(dsc_index_file, the_base_dir)
    vid.populate_index()
    # vik: VisualIndexKwd = VisualIndexKwd(kwd_index_file, the_base_dir)
    # vik.populate_index()
    # process pipeline
    # std_list: List[StandardizedData] = build_std_pipeline_index(the_input_dir, the_index_file, False)
    ##########################################################################
    # collect DSC correctable datasets
    ##########################################################################
    # new corrections to start
    new_correction_jobs_dsc: Dict[str, VisualIndexElementDsc] = {}
    print('*************************************************', flush=True)
    my_dsc: VisualIndexElementDsc
    for my_dsc in vid.m_ele_dict.values():
        if float(my_dsc.m_overall_dsc) > 0.30:
            if '<0.0005' == my_dsc.m_overall_dsc_pvalue:
                # check if already corrected
                corrected: bool = False
                my_data_index_entry: VisualIndexElementBase
                for my_data_index_entry in vib.m_ele_dict.values():
                    if my_data_index_entry.m_path_results == my_dsc.m_path_results:
                        if my_data_index_entry.m_data_version == my_dsc.m_data_version:
                            if "Adjusted-EB_withPara-DSC" == my_data_index_entry.job_type:
                                corrected = True
                if not corrected:
                    # add to list after filtering for largest DSC
                    if my_dsc.m_id not in new_correction_jobs_dsc:
                        new_correction_jobs_dsc[my_dsc.m_id] = my_dsc
                    else:
                        old_dsc: VisualIndexElementDsc = new_correction_jobs_dsc[my_dsc.m_id]
                        if float(my_dsc.m_overall_dsc) > float(old_dsc.m_overall_dsc):
                            new_correction_jobs_dsc[my_dsc.m_id] = my_dsc
    ##########################################################################
    # collect KWD correctable datasets
    ##########################################################################
    # new corrections to start
    # new_correction_jobs_kwd: Dict[str, VisualIndexElementKwd] = {}
    # print('*************************************************', flush=True)
    # my_kwd: VisualIndexElementKwd
    # for my_kwd in vik.m_ele_dict.values():
    #     if float(my_kwd.m_neg_log10_pvalue) >= float(my_kwd.m_neg_log10_cutoff):
    #         # second elements of analysis path is batch type
    #         # my_kwd.m_batches_called is (batch-id-1, batch-id-2)
    #         # check if already corrected
    #         corrected: bool = False
    #         my_data_index_entry: VisualIndexElementBase
    #         for my_data_index_entry in vib.m_ele_dict.values():
    #             if my_data_index_entry.m_path_results == my_kwd.m_path_results:
    #                 if my_data_index_entry.m_data_version == my_kwd.m_data_version:
    #                     if "Adjusted-EB_withPara-KWD" == my_data_index_entry.job_type:
    #                         corrected = True
    #         if not corrected:
    #             # add to list after filtering for largest DSC
    #             if my_kwd.m_id not in new_correction_jobs_kwd:
    #                 new_correction_jobs_kwd[my_kwd.m_id] = my_kwd
    #             else:
    #                 old_kwd: VisualIndexElementKwd = new_correction_jobs_kwd[my_kwd.m_id]
    #                 if float(my_kwd.m_neg_log10_pvalue) > float(old_kwd.m_neg_log10_pvalue):
    #                     new_correction_jobs_kwd[my_kwd.m_id] = my_kwd
    ##########################################################################
    # list correctable datasets
    ##########################################################################
    for my_dsc in new_correction_jobs_dsc.values():
        # Handle Start New CORRECTION Job
        print(f'Start New CORRECTION-DSC Job {round(float(my_dsc.m_overall_dsc),2)} {my_dsc.m_overall_dsc_pvalue} {my_dsc.m_data_version} {my_dsc.m_path_data} {extract_batch(my_dsc.m_analysis_path)}', flush=True)
    # for my_kwd in new_correction_jobs_kwd.values():
    #     # Handle Start New CORRECTION Job
    #     print(f'Start New CORRECTION-KWD Job {round(float(my_kwd.m_neg_log10_pvalue),2)} {my_kwd.m_neg_log10_cutoff} {my_kwd.m_data_version} {my_kwd.m_path_data} {extract_batch(my_kwd.m_analysis_path)} {my_kwd.m_batches_called}', flush=True)
    print('*************************************************', flush=True)
# pylint: enable=too-many-arguments,too-many-locals,too-many-statements,too-many-nested-blocks,too-many-branches
