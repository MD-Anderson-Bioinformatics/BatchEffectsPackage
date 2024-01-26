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

from mbatch.pipeline.pipeline import reindex_pipeline
from mbatch.test.common import print_errors, print_warnings

# update this for internal paths when running for release pipeline
dir_result_data: str = '/BEA/DAPI_MQA/DATA'
dir_result_indexes: str = '/BEA/DAPI_MQA/INDEXES'
dir_temp: str = '/BEA/PIPELINE_TMP'

if __name__ == '__main__':
    # ########################################################
    # rebuild index.json files
    # TODO: collect data to rebuild other indexes too
    # ########################################################
    reindex_pipeline(dir_result_data, dir_temp, "_REINDEX_2023_08_07")
    # ########################################################
    # print important warnings and errors
    # ########################################################
    print("-----------------------------------------", flush=True)
    print("-------------WARNING NOTES-------------", flush=True)
    print_warnings()
    print("-------------ERROR NOTES-------------", flush=True)
    print_errors()
    print("-------------analysis_gdc_reindex:done-------------", flush=True)
