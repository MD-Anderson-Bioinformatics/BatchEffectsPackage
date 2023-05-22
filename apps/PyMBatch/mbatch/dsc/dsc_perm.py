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

from typing import List
import multiprocessing
import numpy
from mbatch.dsc.dsc_info import DscInfo
from mbatch.dsc.dsc_calc import dsc_calc
from mbatch.test.common import handle_error


class DscPerm:
    """
    Class to encapsulate doing permutations of DSC values.
    Permute each row of the matrix and calculate the DSC values.
    Do this the_perms times. Use the_seed as random number generator seed.
    MEMBER VARIABLES
    self.m_seed: int = the_seed - random number generator seed
    self.m_random_num: numpy.random.Generator - random number generator
    self.m_matrix: numpy.ndarray = the_matrix - StdData values
    self.m_batches: numpy.ndarray = the_batches - StdData batches
    self.m_perms: int = the_perms - number of permutations of data to calculate
    self.m_cores: int = the_threads - number of cores or threads to use
    """
    # do not set method variables, as they should be initialized in the init function
    m_seed: int
    m_random_num: numpy.random.Generator
    m_matrix: numpy.ndarray
    m_batches: numpy.ndarray
    m_perms: int
    m_cores: int
    m_counter: int

    # pylint: disable=too-many-arguments
    def __init__(self: 'DscPerm', the_matrix: numpy.ndarray, the_batches: numpy.ndarray, the_seed: int, the_perms: int, the_threads: int) -> None:
        """
        initialize values -- member variables described at class
        :param the_matrix: two-dimensional StdData of values
        :param the_batches: one dimensional array of batches for samples
        :param the_seed: random number generator seed
        :param the_perms: number of permutations of data to calculate
        :param the_threads: number of threads/cores to use
        """
        super().__init__()
        self.m_seed = the_seed
        self.m_random_num = numpy.random.default_rng(seed=the_seed)
        self.m_matrix = the_matrix
        self.m_batches = the_batches
        self.m_perms = the_perms
        self.m_cores = the_threads
        self.m_counter = 0
    # pylint: enable=too-many-arguments

    def perm_dsc_once(self: 'DscPerm') -> DscInfo:
        """
        calculate DSC after permuting matrix, return DscInfo with results
        not thread safe.
        :return: DscInfo with results
        """
        perm_df: numpy.ndarray = self.m_random_num.permuted(self.m_matrix, axis=1)
        info: DscInfo = dsc_calc(perm_df, self.m_batches)
        return info

    def perm_dsc_multi(self: 'DscPerm') -> List[DscInfo]:
        """
        Using member variables, do m_perm number of permutations and DSC calculations.
        Return results as list of DscInfo
        Thread safe.
        :return: list of DscInfo (of permuted matrix)
        """
        info_list: List[DscInfo] = []
        # setup for multiprocessing (real parallel for Python)
        # more than 10 needs a very beefy server
        pool: multiprocessing.Pool
        # pool_manager: multiprocessing.managers.SyncManager = multiprocessing.Manager()
        with multiprocessing.Pool(processes=self.m_cores) as pool:
            # launching multiple evaluations asynchronously *may* use more processes
            _: int
            # with apply_async first argument, function, is copied into async object at base state
            # arguments within square brackets are also copied into async object at base state (unused)
            # see information about how pickle is used. Use ValueProxy to pass updatable, lockable object
            updated_pools: List[multiprocessing.Pool] = \
                [pool.apply_async(dsc_calc, [self.m_random_num.permuted(self.m_matrix, axis=1),
                                             self.m_batches, False],
                                  error_callback=handle_error) for _ in range(self.m_perms)]
            started_proc: multiprocessing.Pool
            for started_proc in updated_pools:
                # print(f"perm_calc_count started_proc={started_proc}", flush=True)
                # no timeout
                info: DscInfo = started_proc.get(None)
                info_list.append(info)
        return info_list

    def perm_only(self: 'DscPerm') -> None:
        """
        Just permute the matrix -- no calculations.
        Not thread safe.
        :return: nothing
        """
        # axis=1 means to permute the values within each row
        # use permuted (permutation and shuffle do different modifications)
        self.m_matrix = self.m_random_num.permuted(self.m_matrix, axis=1)


def dsc_perm_calc_count(the_df: numpy.ndarray, the_batches: numpy.ndarray, the_seed: int, the_perms: int, the_threads: int) -> List[DscInfo]:
    """
    This is used to do the_perms number of permutations of the given dataframe,
    do a DSC calculation on each permutation,
    and return a list of results
    :param the_df: matrix to permute
    :param the_batches: batches for samples
    :param the_seed: seed for random number generator
    :param the_perms: number of permutations to do
    :param the_threads: number of threads/cores to use
    :return: list of DscInfo of permuted matrix
    """
    # print(f"dsc_perm_calc_count the_df={the_df}", flush=True)
    # print(f"dsc_perm_calc_count the_batches={the_batches}", flush=True)
    print(f"dsc_perm_calc_count the_seed={the_seed}", flush=True)
    print(f"dsc_perm_calc_count the_perms={the_perms}", flush=True)
    print(f"dsc_perm_calc_count the_threads={the_threads}", flush=True)
    # NOT THREAD/MULTI-CORE SAFE
    dpp: DscPerm = DscPerm(the_df, the_batches, the_seed, the_perms, the_threads)
    info_list: List[DscInfo] = dpp.perm_dsc_multi()
    return info_list
