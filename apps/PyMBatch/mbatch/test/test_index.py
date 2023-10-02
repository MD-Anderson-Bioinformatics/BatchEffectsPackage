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
from typing import Tuple
from mbatch.test.common import copy_dirs_and_files
from mbatch.index.index import create_index_archive


# files for testing building json and ZIP from configout output
static_test_configout_dir: str = "/BEA/BatchEffectsPackage_data/testing_static/MBatchUtils/configout"
dynamic_test_configout_dir: str = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/index"
configout_result_dir: str = os.path.join(dynamic_test_configout_dir, "ZIP-RESULTS")
configout_info_dir: str = os.path.join(configout_result_dir, "ZIP-RESULTS")
configout_data_dir: str = os.path.join(dynamic_test_configout_dir, "ZIP-DATA")
configout_zip_dir: str = dynamic_test_configout_dir

# files for testing building json and ZIP from error output
static_test_error_dir: str = "/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/errors"
dynamic_test_error_dir: str = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/errors"
error_result_dir: str = os.path.join(dynamic_test_error_dir, "ZIP-RESULTS")
error_info_dir: str = os.path.join(error_result_dir, "ZIP-RESULTS")
error_data_dir: str = os.path.join(dynamic_test_error_dir, "ZIP-DATA")
error_zip_dir: str = dynamic_test_error_dir


def create_index_archive_wrapper(the_results_dir: str, the_data_dir: str, the_zip_dir: str, the_info_dir: str) -> Tuple[str, str]:
    """
    Build index and create zip archive
    :param the_results_dir: directory with MBatch results
    :param the_data_dir: directory with actual data
    :param the_zip_dir: directory in which to place ZIP file
    :param the_info_dir: directory with TEST_<version> entries
    :return: full pathname for ZIP file
    """
    return create_index_archive(the_results_dir, the_data_dir, the_zip_dir, the_info_dir, None, None)


# pylint: disable=too-many-instance-attributes
class TestIndex(unittest.TestCase):
    """
    Class for setting up Index testing - clear/make directory for output
    """
    # do not set method variables, as they should be initialized in the init function
    # No local method variables

    def setUp(self: 'TestIndex') -> None:
        """
        setup script to clear and re-populate test directory
        :return:
        """
        #############################
        # files for configout testing
        if self._testMethodName == 'test_process_configout_dir':
            print(f"TestIndex::setUp dynamic_test_dir={dynamic_test_configout_dir}", flush=True)
            if os.path.exists(dynamic_test_configout_dir):
                shutil.rmtree(dynamic_test_configout_dir)
            print(f"TestIndex::setUp static_test_dir={static_test_configout_dir}", flush=True)
            copy_dirs_and_files(static_test_configout_dir, dynamic_test_configout_dir)
        #############################
        # files for errors testing
        if self._testMethodName == 'test_process_error_dir':
            print(f"TestIndex::setUp dynamic_test_error_dir={dynamic_test_error_dir}", flush=True)
            if os.path.exists(dynamic_test_error_dir):
                shutil.rmtree(dynamic_test_error_dir)
            print(f"TestIndex::setUp static_test_error_dir={static_test_error_dir}", flush=True)
            copy_dirs_and_files(static_test_error_dir, dynamic_test_error_dir)
        #############################
        print("TestIndex:setUp done", flush=True)

    def test_process_configout_dir(self: 'TestIndex') -> None:
        """
        test the create_index_archive function with the configout test directory
        :return: nothing
        """
        print("test_process_configout_dir", flush=True)
        create_index_archive_wrapper(configout_result_dir, configout_data_dir, configout_zip_dir, configout_info_dir)

    def test_process_error_dir(self: 'TestIndex') -> None:
        """
        test the create_index_archive function with the error test directory
        :return: nothing
        """
        print("test_process_error_dir", flush=True)
        create_index_archive_wrapper(error_result_dir, error_data_dir, error_zip_dir, error_info_dir)
# pylint: enable=too-many-instance-attributes


if __name__ == '__main__':
    unittest.main()
