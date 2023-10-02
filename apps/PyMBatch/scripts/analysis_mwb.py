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
from mbatch.pipeline.pipeline import execute_pipeline
# from mbatch.pipeline.pipeline import update_pipeline_data_index, update_pipeline_rebuild_results_archives
from mbatch.test.common import print_errors, print_warnings

dir_datafiles: str = '/BEA/DVLP/Pipeline-MWB/converted'

dir_pipeline_index: str = '/BEA/DVLP/Pipeline-MOB/index/pipline_index.tsv'
dir_pipeline_results: str = '/BEA/DVLP/Pipeline-MOB/results'
dir_pipeline_util: str = '/BEA/DVLP/Pipeline-MOB/util'

server_file: str = '/BEA/DVLP/Pipeline-MOB/pipeline_server.txt'
bei_url: str = ""
bei_dir: str = "/BEA/DVLP/MBA/OUTPUT"
index_base_dir: str = "/BEA/DAPI_MOB/DATA"

run_version: str = "2023_07_01_1200"
run_source: str = "MWB"

if __name__ == '__main__':
    # read URL file
    my_file: io.BufferedReader
    with open(server_file, 'r', encoding='utf-8') as my_file:
        bei_url = str(my_file.read().rstrip())
    # ########################################################
    # pipeline
    # ########################################################
    execute_pipeline(dir_datafiles, dir_pipeline_results, dir_pipeline_index, dir_pipeline_util,
                     bei_url, bei_dir, run_version, run_source, index_base_dir, "Sample")
    # below is used to update pipeline data index after problems.
    # update_pipeline_data_index(dir_datafiles, dir_pipeline_results, dir_pipeline_index, dir_pipeline_util,
    #                            bei_url, bei_dir, run_version, run_source, index_base_dir, "Sample")
    # below is used to rebuild the results ZIP files to fix the internal files
    # update_pipeline_rebuild_results_archives(dir_datafiles, dir_pipeline_results, dir_pipeline_index, dir_pipeline_util,
    #                                          bei_url, bei_dir, run_version, run_source, index_base_dir, "Sample")
    # ########################################################
    # print important warnings and errors
    # ########################################################
    print("-----------------------------------------", flush=True)
    print("-------------WARNING NOTES-------------", flush=True)
    print_warnings()
    print("-------------ERROR NOTES-------------", flush=True)
    print_errors()
    print("-------------pipeline_main:done-------------", flush=True)
