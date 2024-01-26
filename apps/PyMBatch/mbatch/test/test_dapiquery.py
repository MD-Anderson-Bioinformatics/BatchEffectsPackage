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


import unittest
from typing import List
import os
import shutil
import mbatch.dapi.query as dapi


url_base: str = os.environ.get('MBATCH_TEST_DAPI_URL', "")
# url_base = 'https://bioinformatics.mdanderson.org/MQA'
# url_base = 'http://localhost:8080/MQA'
download_dir: str = "/BEA/BatchEffectsPackage_data/testing_dynamic/dapiquery"


# pylint: disable=too-many-instance-attributes
class TestJob(unittest.TestCase):
    """
    Class for setting up Job testing - clear/make directory for output
    """

    # do not set method variables, as they should be initialized in the init function
    # No local method variables

    def setUp(self: 'TestJob') -> None:
        """
        setup script to clear and re-populate test directory
        :return:
        """
        #############################
        print(f"TestJob:setUp clear dir {download_dir}", flush=True)
        if os.path.exists(download_dir):
            shutil.rmtree(download_dir)
        os.makedirs(download_dir)
        print("TestJob:setUp done", flush=True)

    def test_dapi_query_initial(self: 'TestJob') -> None:
        """
        test the Data API
        :return: nothing
        """
        print("==== START ======================================", flush=True)
        if "" == url_base:
            print("no URL, skip test_dapi_query_initial", flush=True)
        else:
            print("test_dapi_query_initial before update", flush=True)
            dapi_query: dapi.DapiQuery = dapi.DapiQuery(url_base)
            dapi_query.update_from_selected()
            print("test_dapi_query_initial after update 1", flush=True)
            dapi_query.print_status(False)
            dapi_query.selected_projects.append('TCGA-LUSC')
            dapi_query.selected_jobtype.append('Original')
            dapi_query.selected_data.append('STAR - Counts')
            dapi_query.update_from_selected()
            print("test_dapi_query_initial after update 2", flush=True)
            dapi_query.print_status(True)
            assert 3 == len(dapi_query.available_datasets), "should have three datasets (may change over time)"
            dapi_entry: dapi.DapiEntry = dapi_query.available_datasets[0]
            dapi_entry.print_status(False)
            data_versions: List[str] = []
            ngchm_files: List[str] = []
            data_versions, ngchm_files = dapi_query.get_downloadable(dapi_entry)
            print(f"test_dapi_query_initial len(data_versions)={len(data_versions)}", flush=True)
            print(f"test_dapi_query_initial len(ngchm_files)={len(ngchm_files)}", flush=True)
            # assert original_str == decrypted_str, "original and decrypted strings are different"
        print("==== END ======================================", flush=True)

    def test_dapi_get_data_matrix(self: 'TestJob') -> None:
        """
        test the Data API
        :return: nothing
        """
        print("==== START ======================================", flush=True)
        if "" == url_base:
            print("no URL, skip test_dapi_get_data_matrix", flush=True)
        else:
            print("test_dapi_get_data_matrix start", flush=True)
            dapi_query: dapi.DapiQuery = dapi.DapiQuery(url_base)
            dapi_query.update_from_selected()
            dataset_id: str = '8ff96a845261f8455b5b2698b5db776f~2022-12-12~2022_12_28_1300'
            date: str = 'DATA_2022-12-12'
            download_file: str = f'{download_dir}matrix_original.tsv'
            # test_dapi_get_data_matrix
            # the_file_path: str, the_id: str, the_version: str, the_original
            dapi_query.download_data_matrix_to_file(download_file, dataset_id, date, True)
            download_file = f'{download_dir}matrix_pipeline.tsv'
            # test_dapi_get_data_matrix
            # the_file_path: str, the_id: str, the_version: str, the_original
            dapi_query.download_data_matrix_to_file(download_file, dataset_id, date, False)
        print("==== END ======================================", flush=True)

    def test_dapi_get_data_batches(self: 'TestJob') -> None:
        """
        test the Data API
        :return: nothing
        """
        print("==== START ======================================", flush=True)
        if "" == url_base:
            print("no URL, skip test_dapi_get_data_batches", flush=True)
        else:
            print("test_dapi_get_data_batches start", flush=True)
            dapi_query: dapi.DapiQuery = dapi.DapiQuery(url_base)
            dapi_query.update_from_selected()
            dataset_id: str = '8ff96a845261f8455b5b2698b5db776f~2022-12-12~2022_12_28_1300'
            date: str = 'DATA_2022-12-12'
            download_file: str = f'{download_dir}batches_original.tsv'
            # test_dapi_get_data_matrix
            # the_file_path: str, the_id: str, the_version: str, the_original
            dapi_query.download_data_batches_to_file(download_file, dataset_id, date, True)
            download_file = f'{download_dir}batches_pipeline.tsv'
            # test_dapi_get_data_matrix
            # the_file_path: str, the_id: str, the_version: str, the_original
            dapi_query.download_data_batches_to_file(download_file, dataset_id, date, False)
        print("==== END ======================================", flush=True)

    def test_dapi_get_ngchm_ngchm(self: 'TestJob') -> None:
        """
        test the Data API
        :return: nothing
        """
        print("==== START ======================================", flush=True)
        if "" == url_base:
            print("no URL, skip test_dapi_get_ngchm_ngchm", flush=True)
        else:
            print("test_dapi_get_ngchm_ngchm start", flush=True)
            dapi_query: dapi.DapiQuery = dapi.DapiQuery(url_base)
            dapi_query.update_from_selected()
            dataset_id: str = '8ff96a845261f8455b5b2698b5db776f~2022-12-12~2022_12_28_1300'
            zip_file: str = '/analysis/NGCHM/DATA_2022-12-12/TEST_2022_12_28_1300/All_ngchm.ngchm.html'
            download_file: str = f'{download_dir}batch_id_ngchm.ngchm'
            dapi_query.download_ngchm_ngchm_to_file(download_file, dataset_id, zip_file)
        print("==== END ======================================", flush=True)

    def test_dapi_get_ngchm_html(self: 'TestJob') -> None:
        """
        test the Data API
        :return: nothing
        """
        print("==== START ======================================", flush=True)
        if "" == url_base:
            print("no URL, skip test_dapi_get_ngchm_html", flush=True)
        else:
            print("test_dapi_get_ngchm_html start", flush=True)
            dapi_query: dapi.DapiQuery = dapi.DapiQuery(url_base)
            dapi_query.update_from_selected()
            dataset_id: str = '8ff96a845261f8455b5b2698b5db776f~2022-12-12~2022_12_28_1300'
            zip_file: str = '/analysis/SupervisedClustering/batch_id/DATA_2022-12-12/TEST_2022_12_28_1300/batch_id_ngchm.ngchm.html'
            download_file: str = f'{download_dir}batch_id_ngchm.html'
            dapi_query.download_ngchm_html_to_file(download_file, dataset_id, zip_file)
        print("==== END ======================================", flush=True)

# pylint: enable=too-many-instance-attributes


if __name__ == '__main__':
    unittest.main()
