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
from mbatch.stddata.stddata import StdData
from mbatch.dsc.dsc_info import DscInfo
from mbatch.dsc.dsc_perm import DscPerm
from mbatch.dsc.dsc_calc import dsc_calc
from mbatch.test.common import generate_file_md5


# TOY DATASET dataframe
# 25 samples, 4 features
M_TOY_DATA: pandas.DataFrame = pandas.DataFrame({
    'Sample01': [1.41642285, 0.10694894, -1.23408485, -0.66788459],
    'Sample02': [-0.72350214, -0.39967497, -1.20617325, 2.02534118],
    'Sample03': [0.04824682, -0.23245758, 0.58790321, -1.62881050],
    'Sample04': [-1.04373111, -0.80310432, 1.20987098, -1.56504974],
    'Sample05': [-0.28991085, -0.18063447, 0.50197599, 1.38438714],
    'Sample06': [-0.28225272, -0.52935371, -0.79628805, -1.99388328],
    'Sample07': [0.61984947, -0.30598282, 0.46818999, -1.07130663],
    'Sample08': [-1.92119735, 0.09369464, -1.37164781, -1.54513036],
    'Sample09': [0.66259059, -0.09702015, -0.81402846, -0.19330046],
    'Sample10': [-1.46182162, 1.29318353, 1.18095836, 1.77052505],
    'Sample11': [1.25097229, 1.51530726, 0.44529141, -1.11088399],
    'Sample12': [-2.28993315, -1.17506637, -0.25687152, 0.53933316],
    'Sample13': [1.67527333, 1.74301675, -0.94686291, -1.46451769],
    'Sample14': [0.44926929, -0.73454580, -0.72579708, -1.33489991],
    'Sample15': [-0.46994309, 0.51241297, 1.12108049, -1.36887916],
    'Sample16': [1.66774028, 0.49402844, 0.36732663, 0.10993039],
    'Sample17': [-0.07856253, -0.20510748, 0.65800235, -0.03943472],
    'Sample18': [-0.57510140, -0.68024807, 2.12295046, -1.03431131],
    'Sample19': [0.24757499, 0.09078522, -0.47703802, 0.43342363],
    'Sample20': [0.29852667, -0.72552589, 0.67071872, 0.68971619],
    'Sample21': [-1.73853724, 0.96191800, 0.80947295, -1.06464023],
    'Sample22': [-0.87065928, -0.54730635, 1.32059294, 0.36257119],
    'Sample23': [2.14670453, -1.05164755, -0.60255080, -0.86457967],
    'Sample24': [0.86925408, -0.27599802, 1.68193731, 0.20491488],
    'Sample25': [-0.02002944, 0.64227494, -0.47946586, -0.15497779]},
    index=['Feature1', 'Feature2', 'Feature3', 'Feature4'])

# TOY DATASET batches
# used for dsc_calc
M_TOY_BATCHES: pandas.Series = pandas.Series([
    "a", "b", "c", "a", "b",
    "a", "b", "c", "a", "b",
    "a", "b", "c", "a", "b",
    "a", "b", "c", "a", "b",
    "a", "b", "c", "a", "b"],
    index=[
        'Sample01', 'Sample02', 'Sample03', 'Sample04', 'Sample05',
        'Sample06', 'Sample07', 'Sample08', 'Sample09', 'Sample10',
        'Sample11', 'Sample12', 'Sample13', 'Sample14', 'Sample15',
        'Sample16', 'Sample17', 'Sample18', 'Sample19', 'Sample20',
        'Sample21', 'Sample22', 'Sample23', 'Sample24', 'Sample25'])

# PERMUTATION ONLY TOY DATASET
# used for dsc_perm 8 samples, four features (Easier to test by eye)
M_TOY_PERM: pandas.DataFrame = pandas.DataFrame({
    'Sample01': [1.1, 1.2, 1.3, 1.4],
    'Sample02': [2.1, 2.2, 2.3, 2.4],
    'Sample03': [3.1, 3.2, 3.3, 3.4],
    'Sample04': [4.1, 4.2, 4.3, 4.4],
    'Sample05': [5.1, 5.2, 5.3, 5.4],
    'Sample06': [6.1, 6.2, 6.3, 6.4],
    'Sample07': [7.1, 7.2, 7.3, 7.4],
    'Sample08': [8.1, 8.2, 8.3, 8.4]},
    index=['Feature1', 'Feature2', 'Feature3', 'Feature4'])


class TestDsc(unittest.TestCase):
    """
    Class for setting up Dsc testing - clear/make directory for output
    """
    # do not set method variables, as they should be initialized in the init function
    # pylint: disable=too-many-instance-attributes
    # it needs as many variables as it needs
    # data files and seed
    sta_matrix: str
    sta_batches: str
    seed: int
    # test_dsc_calc_toy
    sta_toy: str
    dyn_toy: str
    # test_dsc_calc_file
    sta_file: str
    dyn_file: str
    # test_dsc_perm_toy
    sta_perm_toy: str
    dyn_perm_toy: str
    # test_dsc_perm_file
    sta_perm_file: str
    dyn_perm_file: str
    # test_dsc_count_toy
    sta_count_toy: str
    dyn_count_toy: str
    # test_dsc_count_file
    sta_count_file: str
    dyn_count_file: str

    def setUp(self: 'TestDsc') -> None:
        """
        Called automatically by unit testing framework
        :return: nothing
        """
        # data files and seed
        self.sta_matrix = "/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/StdData/TCGA/matrix_data.tsv"
        self.sta_batches = "/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/StdData/TCGA/batches.tsv"
        self.seed: int = 314
        # test_dsc_calc_toy
        self.sta_toy = "/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/dsc/dsc_toy.txt"
        self.dyn_toy = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/dsc/dsc_toy.txt"
        # test_dsc_calc_file
        self.sta_file = "/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/dsc/dsc_file.txt"
        self.dyn_file = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/dsc/dsc_file.txt"
        # test_dsc_perm_toy
        self.sta_perm_toy = "/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/dsc/perm_toy.txt"
        self.dyn_perm_toy = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/dsc/perm_toy.txt"
        # test_dsc_perm_file
        self.sta_perm_file = "/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/dsc/perm_file.txt"
        self.dyn_perm_file = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/dsc/perm_file.txt"
        # test_dsc_count_toy
        self.sta_count_toy = "/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/dsc/count_toy.txt"
        self.dyn_count_toy = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/dsc/count_toy.txt"
        # test_dsc_count_file
        self.sta_count_file = "/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/dsc/count_file.txt"
        self.dyn_count_file = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/dsc/count_file.txt"
        # pylint: enable=too-many-instance-attributes
        if os.path.exists(os.path.dirname(self.dyn_toy)):
            shutil.rmtree(os.path.dirname(self.dyn_toy))
        os.makedirs(os.path.dirname(self.dyn_toy))

    def test_perm_only_toy(self: 'TestDsc') -> None:
        the_sta_toy: str = self.sta_perm_toy
        the_dyn_toy: str = self.dyn_perm_toy
        the_seed: int = self.seed
        print("test_perm_only_toy", flush=True)
        print(f"test_perm_only_toy the_sta_toy={the_sta_toy}", flush=True)
        print(f"test_perm_only_toy the_dyn_toy={the_dyn_toy}", flush=True)
        print(f"test_perm_only_toy the_seed={the_seed}", flush=True)
        mydata: StdData = StdData(M_TOY_PERM)
        dpp: DscPerm = DscPerm(mydata.m_matrix, M_TOY_BATCHES.to_numpy(dtype=str), the_seed, 0, 1)
        dpp.perm_only()
        mydata.m_matrix = dpp.m_matrix
        mydata.write_matrix_data(the_dyn_toy)
        print("test_dsc_calc_count_toy calculate md5s", flush=True)
        dyn_md5: str = generate_file_md5(the_dyn_toy)
        print(f"test_dsc_calc_count_toy dyn_md5={dyn_md5}", flush=True)
        sta_md5: str = generate_file_md5(the_sta_toy)
        print(f"test_dsc_calc_count_toy sta_md5={sta_md5}", flush=True)
        assert dyn_md5 == sta_md5, "Calculated and historic have different values"
        print("test_dsc_perm_toy passed", flush=True)

    def test_perm_only_file(self: 'TestDsc') -> None:
        the_matrix: str = self.sta_matrix
        the_sta_file: str = self.sta_perm_file
        the_dyn_file: str = self.dyn_perm_file
        the_seed: int = self.seed
        print("test_dsc_perm_file", flush=True)
        print(f"test_dsc_perm_file the_sta_file={the_sta_file}", flush=True)
        print(f"test_dsc_perm_file the_dyn_toy={the_dyn_file}", flush=True)
        print(f"test_dsc_perm_file the_seed={the_seed}", flush=True)
        my_std_data: StdData = StdData()
        # read StdData DF
        my_std_data.read_matrix_data(the_matrix)
        dpp: DscPerm = DscPerm(my_std_data.m_matrix, numpy.empty((0, 0), str), the_seed, 0, 1)
        dpp.perm_only()
        my_std_data.m_matrix = dpp.m_matrix
        my_std_data.write_matrix_data(the_dyn_file)
        print("test_dsc_perm_file calculate md5s", flush=True)
        dyn_md5: str = generate_file_md5(the_dyn_file)
        print(f"test_dsc_perm_file dyn_md5={dyn_md5}", flush=True)
        sta_md5: str = generate_file_md5(the_sta_file)
        print(f"test_dsc_perm_file sta_md5={sta_md5}", flush=True)
        assert dyn_md5 == sta_md5, "Calculated and historic have different values"
        print("test_dsc_perm_file passed", flush=True)

    def test_dsc_once_toy(self: 'TestDsc') -> None:
        the_sta_toy: str = self.sta_toy
        the_dyn_toy: str = self.dyn_toy
        print("test_dsc_calc_toy", flush=True)
        print(f"test_dsc_calc_toy the_sta_toy={the_sta_toy}", flush=True)
        print(f"test_dsc_calc_toy the_dyn_toy={the_dyn_toy}", flush=True)
        mydata: StdData = StdData(M_TOY_DATA)
        # print(f"TOY DATASET 25 samples, 4 features", flush=True)
        # print(f"mydata.m_matrix.shape={mydata.m_matrix.shape}", flush=True)
        result: DscInfo = dsc_calc(mydata.m_matrix, M_TOY_BATCHES.to_numpy(dtype=str))
        result.write_to_file(the_dyn_toy)
        print("test_dsc_calc_toy calculate md5s", flush=True)
        dyn_md5: str = generate_file_md5(the_dyn_toy)
        print(f"test_dsc_calc_toy dyn_md5={dyn_md5}", flush=True)
        sta_md5: str = generate_file_md5(the_sta_toy)
        print(f"test_dsc_calc_toy sta_md5={sta_md5}", flush=True)
        assert dyn_md5 == sta_md5, "Calculated and historic have different values"
        print("test_dsc_calc_toy passed", flush=True)

    def test_dsc_once_file(self: 'TestDsc') -> None:
        the_matrix: str = self.sta_matrix
        the_batches: str = self.sta_batches
        the_sta_file: str = self.sta_file
        the_dyn_file: str = self.dyn_file
        print(f"test_dsc_calc_file the_matrix={the_matrix}", flush=True)
        print(f"test_dsc_calc_file the_batches={the_batches}", flush=True)
        print(f"test_dsc_calc_file the_sta_file={the_sta_file}", flush=True)
        print(f"test_dsc_calc_file the_dyn_file={the_dyn_file}", flush=True)
        # read StdData DF
        mydata: StdData = StdData()
        mydata.read_matrix_data(the_matrix)
        mydata.read_batches_data(the_batches)
        my_batch_list: numpy.ndarray = mydata.get_batch_data_for_column('ShipDate')
        result: DscInfo = dsc_calc(mydata.m_matrix, my_batch_list)
        result.write_to_file(the_dyn_file)
        # print(result, flush=True)
        print("test_dsc_calc_file calculate md5s", flush=True)
        dyn_md5: str = generate_file_md5(the_dyn_file)
        print(f"test_dsc_calc_file dyn_md5={dyn_md5}", flush=True)
        sta_md5: str = generate_file_md5(the_sta_file)
        print(f"test_dsc_calc_file sta_md5={sta_md5}", flush=True)
        assert dyn_md5 == sta_md5, "Calculated and historic have different values"
        print("test_dsc_calc_file passed", flush=True)

    # noinspection DuplicatedCode
    def test_dsc_multi_toy(self: 'TestDsc') -> None:
        the_sta_toy: str = self.sta_count_toy
        the_dyn_toy: str = self.dyn_count_toy
        the_seed: int = self.seed
        the_perms: int = 20
        print("test_dsc_multi_toy", flush=True)
        print(f"test_dsc_multi_toy the_sta_toy={the_sta_toy}", flush=True)
        print(f"test_dsc_multi_toy the_dyn_toy={the_dyn_toy}", flush=True)
        print(f"test_dsc_multi_toy the_seed={the_seed}", flush=True)
        print(f"test_dsc_multi_toy the_perms={the_perms}", flush=True)
        print("test_dsc_multi_toy call dsc_perm_calc_count", flush=True)
        mydata: StdData = StdData(M_TOY_DATA)
        dpp: DscPerm = DscPerm(mydata.m_matrix, M_TOY_BATCHES.to_numpy(dtype=str), the_seed, the_perms, 10)
        info_list: List[DscInfo] = dpp.perm_dsc_multi()
        info: DscInfo
        flag: str = 'w'
        print(f"test_dsc_multi_toy write to={the_dyn_toy}", flush=True)
        for info in info_list:
            # print(f"test_dsc_count_toy write flag={flag}", flush=True)
            # print(f"test_dsc_count_toy write info={info}", flush=True)
            info.write_to_file(the_dyn_toy, flag)
            if 'w' == flag:
                flag = 'a'
        print("test_dsc_multi_toy calculate md5s", flush=True)
        dyn_md5: str = generate_file_md5(the_dyn_toy)
        print(f"test_dsc_multi_toy dyn_md5={dyn_md5}", flush=True)
        sta_md5: str = generate_file_md5(the_sta_toy)
        print(f"test_dsc_multi_toy sta_md5={sta_md5}", flush=True)
        assert dyn_md5 == sta_md5, "Calculated and historic have different values"
        print("test_dsc_multi_toy passed", flush=True)

    # pylint: disable=too-many-arguments
    # noinspection DuplicatedCode
    def test_dsc_multi_file(self: 'TestDsc') -> None:
        the_matrix: str = self.sta_matrix
        the_batches: str = self.sta_batches
        the_sta_file: str = self.sta_count_file
        the_dyn_file: str = self.dyn_count_file
        the_seed: int = self.seed
        the_perms: int = 2
        print("test_dsc_multi_file", flush=True)
        print(f"test_dsc_multi_file the_matrix={the_matrix}", flush=True)
        print(f"test_dsc_multi_file the_batches={the_batches}", flush=True)
        print(f"test_dsc_multi_file the_sta_file={the_sta_file}", flush=True)
        print(f"test_dsc_multi_file the_dyn_file={the_dyn_file}", flush=True)
        print(f"test_dsc_multi_file the_seed={the_seed}", flush=True)
        print(f"test_dsc_multi_file the_perms={the_perms}", flush=True)
        print("test_dsc_multi_file call dsc_perm_calc_count", flush=True)
        # read StdData DF
        my_std_data: StdData = StdData()
        my_std_data.read_matrix_data(the_matrix)
        my_std_data.read_batches_data(the_batches)
        # build batch List
        my_batch_list: numpy.ndarray = my_std_data.get_batch_data_for_column('ShipDate')
        # call dsc permutations
        dpp: DscPerm = DscPerm(my_std_data.m_matrix, my_batch_list, the_seed, the_perms, 10)
        info_list: List[DscInfo] = dpp.perm_dsc_multi()
        info: DscInfo
        flag: str = 'w'
        print(f"test_dsc_multi_file write to={the_dyn_file}", flush=True)
        for info in info_list:
            # print(f"test_dsc_calc_count_file write flag={flag}", flush=True)
            # print(f"test_dsc_calc_count_file write info={info}", flush=True)
            info.write_to_file(the_dyn_file, flag)
            if 'w' == flag:
                flag = 'a'
        print("test_dsc_multi_file calculate md5s", flush=True)
        dyn_md5: str = generate_file_md5(the_dyn_file)
        print(f"test_dsc_multi_file dyn_md5={dyn_md5}", flush=True)
        sta_md5: str = generate_file_md5(the_sta_file)
        print(f"test_dsc_multi_file sta_md5={sta_md5}", flush=True)
        assert dyn_md5 == sta_md5, "Calculated and historic have different values"
        print("test_dsc_multi_file passed", flush=True)
    # pylint: enable=too-many-arguments


if __name__ == '__main__':
    unittest.main()
