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


import unittest
import shutil
import os
from mbatch.pipeline.job import copy_json_original_files
from mbatch.index.index_original_data import read_json_original_data


# files for testing building json and ZIP from configout output
static_test_job_dir: str = "/BatchEffectsPackage_data/testing_static/PyMBatch/job"

# files for testing building json and ZIP from error output
dynamic_test_job_dir: str = "/BatchEffectsPackage_data/testing_dynamic/PyMBatch/job"


def copy_original_data_wrapper(the_dynamic_test_job_dir: str, the_static_test_job_dir: str) -> 'None':
    """
    Create new original_data.json from standardized data version.
    Then read file to makes sure format is correct.
    :param the_dynamic_test_job_dir: directory to write results
    :param the_static_test_job_dir: directory to read original file
    :return: nothing
    """
    std_json: str = os.path.join(the_static_test_job_dir, "index.zip")
    print(f"copy_original_data_wrapper read from std_json={std_json}", flush=True)
    copy_json_original_files(std_json, "test-source", "test-version", the_dynamic_test_job_dir)
    job_json: str = os.path.join(the_dynamic_test_job_dir, "original_data.json")
    print(f"copy_original_data_wrapper test from job_json={job_json}", flush=True)
    read_json_original_data(job_json)


# pylint: disable=too-many-instance-attributes
class TestJob(unittest.TestCase):
    """
    Class for setting up Job testing - clear/make directory for output
    Currently, only tests copying job data to new original_data.json
    """
    # do not set method variables, as they should be initialized in the init function
    # No local method variables

    def setUp(self: 'TestJob') -> None:
        """
        setup script to clear and re-populate test directory
        :return:
        """
        #############################
        # files for configout testing
        if self._testMethodName == 'test_process_configout_dir':
            print(f"TestJob::setUp dynamic_test_dir={dynamic_test_job_dir}", flush=True)
            if os.path.exists(dynamic_test_job_dir):
                shutil.rmtree(dynamic_test_job_dir)
            os.makedirs(dynamic_test_job_dir)
            print(f"TestJob::setUp static_test_dir={static_test_job_dir}", flush=True)
        #############################
        print("TestJob:setUp done", flush=True)

    def test_process_copy_original_data(self: 'TestJob') -> None:
        """
        test the create_index_archive function with the configout test directory
        :return: nothing
        """
        print("test_process_copy_original_data", flush=True)
        copy_original_data_wrapper(dynamic_test_job_dir, static_test_job_dir)
# pylint: enable=too-many-instance-attributes


if __name__ == '__main__':
    unittest.main()
