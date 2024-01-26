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


import csv
import typing
import pandas
import numpy


class StdData:
    """
    Class to hold StdData with row and column names.
    But using numpy for faster math routines.
    MEMBER VALUES
    self.m_matrix: numpy.ndarray - floating point data
    self.m_samples: numpy.ndarray - sample ids
    self.m_features: numpy.ndarray - feature ids
    self.m_batches: numpy.ndarray - batch data with each row is a sample
    self.m_columns: numpy.ndarray = column labels for batches. First column is Sample
    """
    # do not set method variables, as they should be initialized in the init function
    m_matrix: numpy.ndarray
    m_samples: numpy.ndarray
    m_features: numpy.ndarray
    # Batch information - first column is Sample
    # rows are sorted by values in Sample column
    m_batches: numpy.ndarray
    m_columns: numpy.ndarray

    def __init__(self: 'StdData', the_matrix: typing.Optional[pandas.DataFrame] = None) -> None:
        """
        init empty/nan values.
        type hinting for class or none is from https://stackoverflow.com/questions/19202633/python-3-type-hinting-for-none
        Members described at class level
        """
        # Matrix data fields
        self.m_matrix = numpy.empty((0, 0), float)
        self.m_samples = numpy.empty((0, 0), str)
        self.m_features = numpy.empty((0, 0), str)
        # Batch information - first column is Sample
        # rows are sorted by values in Sample column
        self.m_batches = numpy.empty((0, 0), str)
        self.m_columns = numpy.empty((0, 0), str)
        # check if the_matrix is passed in
        if the_matrix is not None:
            self.process_panda_matrix(the_matrix)

    def process_panda_matrix(self: 'StdData', the_matrix: pandas.DataFrame) -> None:
        """
        Convert Pandas dataframe to numpy
        :param the_matrix: pandas DataFrame
        :return: nothing
        """
        the_matrix.apply(pandas.to_numeric)
        # print(f"the_matrix.shape={the_matrix.shape}", flush=True)
        # index_sort() sorts by row names
        # index_sort(axis=1) sorts by column names
        the_matrix = the_matrix.sort_index().sort_index(axis=1)
        self.m_matrix = the_matrix.to_numpy()
        # print(f"self.m_matrix.shape={self.m_matrix.shape}", flush=True)
        self.m_samples = the_matrix.columns.to_numpy(dtype=str)
        self.m_features = the_matrix.index.to_numpy(dtype=str)

    def read_matrix_data(self: 'StdData', the_file: str) -> None:
        """
        Load standardized data from file into StdData object.
        Uses Pandas for easy loading. But stores in numpy
        :param the_file: matrix file to load. Samples as columns with leading tab, features as rows
        :return: nothing
        """
        my_matrix: pandas.DataFrame = pandas.read_csv(the_file, sep='\t', quoting=csv.QUOTE_NONE,
                                                      encoding='utf-8', index_col=0)
        self.process_panda_matrix(my_matrix)

    def read_batches_data(self: 'StdData', the_file: str, the_sample_id_col: str = "Sample") -> None:
        """
        Load standardized data from file into StdData object.
        Uses Pandas for easy loading. But stores in numpy
        :param the_file: batches data file to load, with Sample as first column, and rows sorted by sample
        :param the_sample_id_col: string name of sample id column
        :return: nothing
        """
        my_batches: pandas.DataFrame = pandas.read_csv(the_file, sep='\t', quoting=csv.QUOTE_NONE,
                                                       encoding='utf-8', index_col=0, dtype=str)
        my_batches = my_batches.sort_values(by=the_sample_id_col)
        self.m_batches = my_batches.to_numpy(dtype=str)
        self.m_columns = my_batches.columns.to_numpy(dtype=str)
        # sample is not in .columns so needs to be added
        self.m_columns = numpy.insert(self.m_columns, 0, the_sample_id_col)

    def write_matrix_data(self: 'StdData', the_file: str) -> None:
        """
        Write file from standardized data from **numpy** objects in StdData object.
        Uses Pandas for easy writing. But stores in numpy
        :param the_file: full path to matrix file to write. Samples as columns with leading tab, features as rows
        :return: nothing
        """
        my_matrix: pandas.DataFrame = pandas.DataFrame(self.m_matrix, self.m_features, self.m_samples)
        my_matrix.to_csv(the_file, sep='\t', encoding='utf-8')

    def get_batch_data_for_column(self: 'StdData', the_value: str) -> numpy.array:
        """
        Batch type to get values for
        :param the_value: batch type name
        :return: numpy.ndarray of batch values
        """
        index: int = numpy.where(self.m_columns == the_value)[0][0]
        # subtract one, since numpy doesn't use Sample as a "real" column
        index = index - 1
        values: numpy.ndarray = self.m_batches[:, index]
        return values
