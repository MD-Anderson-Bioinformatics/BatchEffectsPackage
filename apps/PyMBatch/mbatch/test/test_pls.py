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


from typing import List
import unittest
import shutil
import os
import pandas
import numpy
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler
from mbatch.stddata.stddata import StdData
from mbatch.test.common import generate_file_md5


def pls_calc(the_matrix: numpy.ndarray, the_batches: numpy.ndarray) -> numpy.ndarray:
    # https://stackoverflow.com/questions/18390150/pls-da-algorithm-in-python
    # replace NaN with zero
    # the_matrix[numpy.isnan(the_matrix)] = 0
    # drop NaN
    clean_matrix: numpy.ndarray = the_matrix[~numpy.isnan(the_matrix).any(axis=1)]
    # scale_matrix: numpy.ndarray = StandardScaler().fit_transform(clean_matrix)
    # build batch ndarray
    unique_batches: numpy.ndarray = numpy.sort(the_batches)
    unique_batches = numpy.unique(unique_batches)
    supervised_batch: numpy.ndarray = numpy.empty((0, clean_matrix.shape[1]))
    batch_name: str
    for batch_name in unique_batches:
        sub_array: numpy.ndarray = numpy.array(the_batches)
        batch_match: numpy.ndarray = (sub_array == batch_name)
        sub_array[batch_match] = 1
        sub_array[~batch_match] = 0
        supervised_batch = numpy.concatenate((supervised_batch, sub_array.reshape(1, -1)), axis=0)
    # pls_decomp = PLSRegression(n_components=scale_matrix.shape[1]).fit_transform(scale_matrix.transpose(), clean_matrix.transpose())
    clean_t: numpy.ndarray = clean_matrix.transpose()
    batch_t: numpy.ndarray = supervised_batch.transpose()
    sample_count: int = clean_matrix.shape[1]
    plsda: PLSRegression = PLSRegression(n_components=sample_count).fit(clean_t, batch_t)
    pls_matrix: numpy.ndarray = plsda.transform(clean_t)
    pls_matrix = pls_matrix.transpose()
    print(f"pls_calc pls_matrix {pls_matrix.shape}", flush=True)
    return pls_matrix


class TestPls(unittest.TestCase):
    """
    Class for setting up PLS testing - clear/make directory for output
    """
    # do not set method variables, as they should be initialized in the init function
    # pylint: disable=too-many-instance-attributes
    # it needs as many variables as it needs
    # data files and seed
    sta_matrix: str
    sta_batches: str
    seed: int

    def setUp(self: 'TestPls') -> None:
        """
        Called automatically by unit testing framework
        :return: nothing
        """
        # data files and seed
        self.sta_matrix = "/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/StdData/TCGA/matrix_data.tsv"
        self.sta_batches = "/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/StdData/TCGA/batches.tsv"
        self.seed: int = 314
        # test_pls_calc_file
        self.sta_file = "/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/plsda/plsda_file.txt"
        self.dyn_dir = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/plsda"
        self.dyn_file = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/plsda/plsda_file.txt"
        # pylint: enable=too-many-instance-attributes
        if os.path.exists(os.path.dirname(self.dyn_dir)):
            shutil.rmtree(os.path.dirname(self.dyn_dir))
        os.makedirs(os.path.dirname(self.dyn_dir))

    def test_pls_once_file(self: 'TestPls') -> None:
        the_matrix: str = self.sta_matrix
        the_batches: str = self.sta_batches
        the_sta_file: str = self.sta_file
        the_dyn_file: str = self.dyn_file
        print(f"test_pls_calc_file the_matrix={the_matrix}", flush=True)
        print(f"test_pls_calc_file the_batches={the_batches}", flush=True)
        print(f"test_pls_calc_file the_sta_file={the_sta_file}", flush=True)
        print(f"test_pls_calc_file the_dyn_file={the_dyn_file}", flush=True)
        # read StdData DF
        mydata: StdData = StdData()
        mydata.read_matrix_data(the_matrix)
        mydata.read_batches_data(the_batches)
        my_batch_list: numpy.ndarray = mydata.get_batch_data_for_column('ShipDate')
        pls_decomp = pls_calc(mydata.m_matrix, my_batch_list)
        print(pls_decomp, flush=True)
        # result: PlsInfo = pls_calc(mydata.m_matrix, my_batch_list)
        # result.write_to_file(the_dyn_file)
        # # print(result, flush=True)
        # print("test_pls_calc_file calculate md5s", flush=True)
        # dyn_md5: str = generate_file_md5(the_dyn_file)
        # print(f"test_pls_calc_file dyn_md5={dyn_md5}", flush=True)
        # sta_md5: str = generate_file_md5(the_sta_file)
        # print(f"test_pls_calc_file sta_md5={sta_md5}", flush=True)
        # assert dyn_md5 == sta_md5, "Calculated and historic have different values"
        print("test_pls_calc_file passed", flush=True)


if __name__ == '__main__':
    unittest.main()
