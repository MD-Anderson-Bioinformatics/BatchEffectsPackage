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
import os
import csv
import typing
import pandas


def write_converted_data(the_matrix: pandas.DataFrame, the_sample_list: List[str],
                         the_program_batches: str, the_program_clinical: str, the_output_dir: str) -> None:
    """
    Write the_matrix to the_filename within the_output_dir.
    Sort the index before writing.
    Matrix is passed in, but dataset is used to read program/project level
    biospecimen information and filter into just this batch data.
    :param the_matrix: pandas.DataFrame to write
    :param the_sample_list: list of samples used in matrix.
    :param the_program_batches: Filename with path for batch information to read.
    :param the_program_clinical: Filename with path for clinical information to read.
    :param the_output_dir: Full path to directory into which to write file. Create if needed.
    :return: Nothing
    """
    if not os.path.exists(the_output_dir):
        os.makedirs(the_output_dir)
    # sort rows and columns of DataFrame
    the_matrix.sort_index(axis=0, inplace=True)
    the_matrix.sort_index(axis=1, inplace=True)
    the_sample_list.sort()
    # write matrix file
    matrix_tsv_file: str = os.path.join(the_output_dir, "matrix.tsv")
    the_matrix.to_csv(matrix_tsv_file, sep='\t', encoding='utf-8', na_rep='NA')
    #
    # write batches file
    if the_program_batches is not None:
        if os.path.exists(the_program_batches):
            clinical_tsv_file: str = os.path.join(the_output_dir, "clinical.tsv")
            batch_tsv_file: str = os.path.join(the_output_dir, "batches.tsv")
            my_batches: pandas.DataFrame = pandas.read_csv(the_program_batches, sep='\t', quoting=csv.QUOTE_NONE, encoding='utf-8', index_col=0, dtype=str)
            my_batches = my_batches.loc[my_batches.index.isin(the_sample_list)]
            my_batches.to_csv(batch_tsv_file, sep='\t', encoding='utf-8', na_rep='Unknown')
            my_patients_series: pandas.Series = my_batches['patient_barcode']
            my_patients: List[str] = my_patients_series.to_list()
            # write clinical file
            if the_program_clinical is not None:
                if os.path.exists(the_program_clinical):
                    my_clinical: pandas.DataFrame = pandas.read_csv(the_program_clinical, sep='\t', quoting=csv.QUOTE_NONE, encoding='utf-8', index_col=0, dtype=str)
                    my_clinical = my_clinical.loc[my_clinical.index.isin(my_patients)]
                    my_clinical.to_csv(clinical_tsv_file, sep='\t', encoding='utf-8', na_rep='Unknown')
    print(f"write_converted_data wrote {the_output_dir}", flush=True)


def write_converted_dataframe(the_matrix: pandas.DataFrame, the_output_dir: str, the_filename: str,
                              the_sample_column: typing.Optional[str]) -> None:
    """
    Write the_matrix to the_filename within the_output_dir.
    Sort the index before writing.
    :param the_matrix: pandas.DataFrame to write
    :param the_output_dir: Full path to directory into which to write file. Create if needed.
    :param the_filename: Filename (no path) to write.
    :param the_sample_column: string column name to move to front of dataframe
    :return: Nothing
    """
    if not os.path.exists(the_output_dir):
        os.makedirs(the_output_dir)
    # sort rows and columns of DataFrame
    the_matrix.sort_index(axis=0, inplace=True)
    the_matrix.sort_index(axis=1, inplace=True)
    # move sample column to be first, if provided
    if the_sample_column is not None:
        column_to_move: pandas.Series = the_matrix.pop(the_sample_column)
        the_matrix.insert(0, the_sample_column, column_to_move)
    # write matrix file
    matrix_tsv_file: str = os.path.join(the_output_dir, the_filename)
    # index=False because index has been added as a normal column
    the_matrix.to_csv(matrix_tsv_file, sep='\t', encoding='utf-8', na_rep='NA', index=False)
