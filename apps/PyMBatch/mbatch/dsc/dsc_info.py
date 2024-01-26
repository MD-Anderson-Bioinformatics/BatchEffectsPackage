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
from textwrap import dedent
import math
import io


def convert_to_list(the_string: str) -> List[float]:
    """
    Convert a string version of a float list [ x1, x2, ... x3 ] into a real list.
    Used for reading DscInfo objects written to files
    :param the_string: string version of a float list [ x1, x2, ... x3 ]
    :return: list of float built from string
    """
    result: List[float] = []
    # trim [ and ] off front and back
    the_string = the_string[1:]
    the_string = the_string[:-2]
    # parse and iterate over line converting to float
    float_str: str
    for float_str in the_string.split(", "):
        result.append(float(float_str))
    return result


def epsilon_zero_check_value(the_value: float) -> float:
    """
    if the value is less than or equal 1x10^-7, use 0.0, otherwise use the value
    :param the_value: value to test
    :return: if the value is less than or equal 1x10^-7, use 0.0, otherwise use the value
    """
    if (not math.isnan(the_value)) & (abs(the_value) <= 0.0000001):
        the_value = 0.0
    return the_value


class DscInfo:
    """
    Class to hold values computed for DSC
    For more details see: https://bioinformatics.mdanderson.org/public-software/tcga-batch-effects/#the-dsc-metric
    Dispersion Separability Criterion
    Higher means more dispersion (more chance of batch effects)
    MEMBER VALUES
    m_list_of_feature_dsc: List[float] - dispersion by feature
    m_list_of_feature_db: List[float] - dispersion between batches by feature
    m_list_of_feature_dw: List[float] - dispersion within batches by feature
    m_dsc: float - dispersion for dataset
    m_db: float - dispersion between batches
    m_dw: float - dispersion within batches
    """
    # do not set method variables, as they should be initialized in the init function
    m_list_of_feature_dsc: List[float]
    m_list_of_feature_db: List[float]
    m_list_of_feature_dw: List[float]
    m_dsc: float
    m_db: float
    m_dw: float

    def __init__(self: 'DscInfo') -> None:
        """
        init empty/nan values.
        Members described at class level
        """
        super().__init__()
        self.m_list_of_feature_dsc = []
        self.m_list_of_feature_db = []
        self.m_list_of_feature_dw = []
        self.m_dsc = float('nan')
        self.m_db = float('nan')
        self.m_dw = float('nan')

    def __str__(self: 'DscInfo') -> str:
        """
        tostring makes a string with data
        :return: String version of object and its data
        """
        return dedent(f"""
            {super().__str__()}
            m_dsc = {self.m_dsc}
            m_db = {self.m_db}
            m_dw = {self.m_dw}
            m_list_of_feature_dsc = {self.m_list_of_feature_dsc[:3]}
            m_list_of_feature_db = {self.m_list_of_feature_db[:3]}
            m_list_of_feature_dw = {self.m_list_of_feature_dw[:3]}""")

    def epsilon_zero_check(self: 'DscInfo') -> None:
        """
        Check member variables for "close enough to zero"
        """
        self.m_dsc = epsilon_zero_check_value(self.m_dsc)
        self.m_db = epsilon_zero_check_value(self.m_db)
        self.m_dw = epsilon_zero_check_value(self.m_dw)
        # use this method (rather than list/map), this keeps type checking properly
        value: float
        self.m_list_of_feature_dsc = [epsilon_zero_check_value(value) for value in self.m_list_of_feature_dsc]
        self.m_list_of_feature_db = [epsilon_zero_check_value(value) for value in self.m_list_of_feature_db]
        self.m_list_of_feature_dw = [epsilon_zero_check_value(value) for value in self.m_list_of_feature_dw]

    def __eq__(self: 'DscInfo', the_other: 'DscInfo') -> bool:
        """
        Check that all values are equal between object.
        :param the_other: another DscInfo object
        :return: True is all values match
        """
        is_equal: bool = True
        if self.m_dsc != the_other.m_dsc:
            is_equal = False
        elif self.m_db != the_other.m_db:
            is_equal = False
        elif self.m_dw != the_other.m_dw:
            is_equal = False
        elif self.m_list_of_feature_dsc != the_other.m_list_of_feature_dsc:
            is_equal = False
        elif self.m_list_of_feature_db != the_other.m_list_of_feature_db:
            is_equal = False
        elif self.m_list_of_feature_dw != the_other.m_list_of_feature_dw:
            is_equal = False
        return is_equal

    def write_to_file(self: 'DscInfo', the_file: str, the_open_flag: str = 'w') -> str:
        """
        Write object to disk file, in manner comparable between runs
        :param the_file: full path to File to which to write this DscInfo object
        :param the_open_flag: flag to 'w' truncate or 'a' append, argument to open
        :return: the_file (file written)
        """
        elem: float
        out_file: io.TextIOWrapper
        with open(the_file, the_open_flag, encoding='utf-8') as out_file:
            out_file.write("<mbatch.dsc.dsc_info.DscInfo>\n")
            out_file.write(f"{self.m_dsc:.8f}\n")
            out_file.write(f"{self.m_db:.8f}\n")
            out_file.write(f"{self.m_dw:.8f}\n")
            # from https://stackoverflow.com/a/45693874
            # pylint: disable=consider-using-f-string
            out_file.write(f"{'[{:s}]'.format(', '.join(['{:.8f}'.format(elem) for elem in self.m_list_of_feature_dsc]))}\n")
            out_file.write(f"{'[{:s}]'.format(', '.join(['{:.8f}'.format(elem) for elem in self.m_list_of_feature_db]))}\n")
            out_file.write(f"{'[{:s}]'.format(', '.join(['{:.8f}'.format(elem) for elem in self.m_list_of_feature_dw]))}\n")
            # pylint: enable=consider-using-f-string
        return the_file

    def read_from_file(self: 'DscInfo', the_file: str) -> str:
        """
        Read file written by write_to_file, populate self from that file
        :param the_file: full path to file to read
        :return: the_file (file read)
        """
        print(f"read_from_file the_file={the_file}", flush=True)
        in_file: io.TextIOWrapper
        line: str
        count: int = 0
        with open(the_file, 'r', encoding='utf-8') as in_file:
            for line in in_file:
                if 1 == count:
                    self.m_dsc = float(line)
                elif 2 == count:
                    self.m_db = float(line)
                elif 3 == count:
                    self.m_dw = float(line)
                elif 4 == count:
                    self.m_list_of_feature_dsc = convert_to_list(line)
                elif 5 == count:
                    self.m_list_of_feature_db = convert_to_list(line)
                elif 6 == count:
                    self.m_list_of_feature_dw = convert_to_list(line)
                count += 1
        return the_file
