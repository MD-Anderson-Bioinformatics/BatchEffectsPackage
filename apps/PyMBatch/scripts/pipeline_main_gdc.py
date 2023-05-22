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

from mbatch.pipeline.pipeline import execute_pipeline, update_pipeline_data_index, update_pipeline_rebuild_results_archives
from mbatch.test.common import print_errors, print_warnings

dir_datafiles: str = '/your/data/path/Pipeline-GDC/converted/data'

dir_pipeline_index: str = '/your/data/path/Pipeline-MQA/index/pipline_index.tsv'
dir_pipeline_results: str = '/your/data/path/Pipeline-MQA/results'
dir_pipeline_util: str = '/your/data/path/Pipeline-MQA/util'

bei_url: str = "http://your.server.com/BEI/BEI/"
bei_dir: str = "/your/data/path/BEI/OUTPUT"
index_base_dir: str = "/DAPI_MQA/DATA"

run_version: str = "2022_12_28_1300"
run_source: str = "GDC"

if __name__ == '__main__':
    # ########################################################
    # pipeline
    # ########################################################
    execute_pipeline(dir_datafiles, dir_pipeline_results, dir_pipeline_index, dir_pipeline_util,
                     bei_url, bei_dir, run_version, run_source, index_base_dir, "aliquot_barcode")
    # below is used to update pipeline data index after problems.
    # update_pipeline_data_index(dir_datafiles, dir_pipeline_results, dir_pipeline_index, dir_pipeline_util,
    #                            bei_url, bei_dir, run_version, run_source, index_base_dir, "aliquot_barcode")
    # below is used to rebuild the results ZIP files to fix the internal files
    # update_pipeline_rebuild_results_archives(dir_datafiles, dir_pipeline_results, dir_pipeline_index, dir_pipeline_util,
    #                                         bei_url, bei_dir, run_version, run_source, index_base_dir, "aliquot_barcode")
    # ########################################################
    # print important warnings and errors
    # ########################################################
    print("-----------------------------------------", flush=True)
    print("-------------WARNING NOTES-------------", flush=True)
    print_warnings()
    print("-------------ERROR NOTES-------------", flush=True)
    print_errors()
    print("-------------pipeline_main:done-------------", flush=True)
