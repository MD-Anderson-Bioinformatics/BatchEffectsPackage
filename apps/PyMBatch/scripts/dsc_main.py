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

from mbatch.test.test_dsc import test_dsc_multi_file, test_dsc_multi_toy
from mbatch.test.test_dsc import test_dsc_once_toy, test_dsc_once_file
from mbatch.test.test_dsc import test_perm_only_file, test_perm_only_toy

# data files and seed
sta_matrix: str = "/BatchEffectsPackage_data/testing_static/PyMBatch/StdData/TCGA/matrix_data.tsv"
sta_batches: str = "/BatchEffectsPackage_data/testing_static/PyMBatch/StdData/TCGA/batches.tsv"
seed: int = 314
# test_dsc_calc_toy
sta_toy: str = "/BatchEffectsPackage_data/testing_static/PyMBatch/dsc/dsc_toy.txt"
dyn_toy: str = "/BatchEffectsPackage_data/testing_dynamic/PyMBatch/dsc/dsc_toy.txt"
# test_dsc_calc_file
sta_file: str = "/BatchEffectsPackage_data/testing_static/PyMBatch/dsc/dsc_file.txt"
dyn_file: str = "/BatchEffectsPackage_data/testing_dynamic/PyMBatch/dsc/dsc_file.txt"
# test_dsc_perm_toy
sta_perm_toy: str = "/BatchEffectsPackage_data/testing_static/PyMBatch/dsc/perm_toy.txt"
dyn_perm_toy: str = "/BatchEffectsPackage_data/testing_dynamic/PyMBatch/dsc/perm_toy.txt"
# test_dsc_perm_file
sta_perm_file: str = "/BatchEffectsPackage_data/testing_static/PyMBatch/dsc/perm_file.txt"
dyn_perm_file: str = "/BatchEffectsPackage_data/testing_dynamic/PyMBatch/dsc/perm_file.txt"
# test_dsc_count_toy
sta_count_toy: str = "/BatchEffectsPackage_data/testing_static/PyMBatch/dsc/count_toy.txt"
dyn_count_toy: str = "/BatchEffectsPackage_data/testing_dynamic/PyMBatch/dsc/count_toy.txt"
# test_dsc_count_file
sta_count_file: str = "/BatchEffectsPackage_data/testing_static/PyMBatch/dsc/count_file.txt"
dyn_count_file: str = "/BatchEffectsPackage_data/testing_dynamic/PyMBatch/dsc/count_file.txt"


if __name__ == '__main__':
    # single dsc calculation on toy dataset
    test_dsc_once_toy(sta_toy, dyn_toy)
    # single dsc calculation on dataset from file
    test_dsc_once_file(sta_matrix, sta_batches, sta_file, dyn_file)
    # permute dataframe of toy dataset
    test_perm_only_toy(sta_perm_toy, dyn_perm_toy, seed)
    # permute dataframe of dataset from file
    test_perm_only_file(sta_matrix, sta_perm_file, dyn_perm_file, seed)
    # do 20 permutations of toy dataset, calculate DSC, and write to file
    test_dsc_multi_toy(sta_count_toy, dyn_count_toy, seed, 20)
    # do 20 permutations of file-based dataset, calculate DSC, and write to file
    test_dsc_multi_file(sta_matrix, sta_batches, sta_count_file, dyn_count_file, seed, 2)
