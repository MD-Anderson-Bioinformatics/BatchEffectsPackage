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

import io
from mbatch.pipeline.pipeline import execute_pipeline
from mbatch.test.common import print_errors, print_warnings

# update_pipeline_data_index, update_pipeline_rebuild_results_archives
# update this for internal paths when running for release pipeline

# NOTE: run version needs to be before pipeline index
RUN_VERSION: str = "2024_01_03_1400"
RUN_SOURCE: str = "GDC"

DIR_DATAFILES: str = '/BEA/DAPI_GDC/converted/data'

DIR_PIPELINE_INDEX: str = f'/BEA/DAPI_MQA/INDEXES/pipline_index_{RUN_VERSION}.tsv'
DIR_PIPELINE_RESULTS: str = '/BEA/DAPI_MQA/DATA'
DIR_PIPELINE_UTIL: str = '/BEA/DAPI_MQA/CONFIG'

SERVER_FILE: str = '/BEA/DAPI_MQA/pipeline_server.tsv'
BEI_URL: str = ''
BEI_DIR: str = "/BEA/BEI/OUTPUT"
INDEX_BASE_DIR: str = "/BEA/DAPI_MQA/DATA"


if __name__ == '__main__':
    # read URL file
    my_file: io.BufferedReader
    with open(SERVER_FILE, 'r', encoding='utf-8') as my_file:
        BEI_URL = str(my_file.read().rstrip())
    # ########################################################
    # pipeline
    # ########################################################
    execute_pipeline(DIR_DATAFILES, DIR_PIPELINE_RESULTS, DIR_PIPELINE_INDEX, DIR_PIPELINE_UTIL,
                     BEI_URL, BEI_DIR, RUN_VERSION, RUN_SOURCE,
                     INDEX_BASE_DIR, "aliquot_barcode", True)
    # ########################################################
    # print important warnings and errors
    # ########################################################
    print("-----------------------------------------", flush=True)
    print("-------------WARNING NOTES-------------", flush=True)
    print_warnings()
    print("-------------ERROR NOTES-------------", flush=True)
    print_errors()
    print("-------------pipeline_main:done-------------", flush=True)
