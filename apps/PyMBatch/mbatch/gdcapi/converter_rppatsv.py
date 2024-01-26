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


import io
import zipfile
import typing
from typing import Dict, List, Tuple
import pandas
from mbatch.gdcapi.convert_dataset import Dataset
from mbatch.gdcapi.download_datafile import GdcApiDatafile


# pylint: disable=too-many-arguments
def convert_rppatsv(the_dataset: Dataset, the_sample_dir: str, the_download_dir: str) -> Tuple[pandas.DataFrame, List[str]]:
    """
    Build the data matrix for this dataset and return it and the list of sample ids (barcodes)
    :param the_dataset: Dataset object describing matrix to build.
    :param the_sample_dir: Full path to sample directory inside biospecimen directory.
    :param the_download_dir: Full path to download/data directory.
    :return: A tuple of the DataFrame (matrix) and List of barcodes (sample ids) in the matrix.
    """
    print(f"convert_rppatsv size={len(the_dataset.files)} for {the_dataset.get_dataset_path('/')}", flush=True)
    # my matrix to be populated
    my_matrix: pandas.DataFrame = pandas.DataFrame({}, dtype='str')
    sample_list: List[str] = []
    # my_batches: pandas.DataFrame = pandas.DataFrame({}, dtype='str')
    # my_clinical: pandas.DataFrame = pandas.DataFrame({}, dtype='str')
    # read headers from RPPA file (which is inside a ZIP archive)
    my_datafile: GdcApiDatafile
    for my_datafile in the_dataset.files.values():
        zip_file: str = my_datafile.get_file_archive(the_download_dir)
        # print(f"convert_rppatsv zip_file={zip_file} file_name={my_datafile.file_name}", flush=True)
        my_datafile.read_sample_file(the_sample_dir)
        assert 1 == len(my_datafile.sample_dict), "convert_rppatsv - only expected one sample"
        barcode: str = list(my_datafile.sample_dict.values())[0].sample_barcode
        my_matrix = read_and_process_file(my_matrix, barcode, zip_file, my_datafile.file_name)
        sample_list.append(barcode)
    # now we have a complete matrix my_matrix
    # and a list of samples in that matrix sample_list
    return my_matrix, sample_list
# pylint: enable=too-many-arguments


# pylint: disable=too-many-arguments,too-many-locals,consider-using-min-builtin,too-many-nested-blocks
def read_and_process_file(the_matrix: pandas.DataFrame, the_barcode: str,
                          the_file_zip: str, the_filename: str) -> pandas.DataFrame:
    """
    Populate the_matrix with next column (sample) of data from file inside ZIP archive.
    :param the_matrix: pandas.DataFrame to which to add a new column
    :param the_barcode: sample id of new data.
    :param the_file_zip: Full path and file name of ZIP file with data file.
    :param the_filename: File name of data inside ZIP.
    :return: Return the_matrix with new column (sample) of data added.
    """
    ret_value: pandas.DataFrame = the_matrix
    zip_file: zipfile.ZipFile
    in_file: io.TextIOWrapper
    headers: typing.Optional[List[str]] = None
    with zipfile.ZipFile(the_file_zip, 'r') as zip_file:
        with zip_file.open(the_filename, mode="r") as in_file:
            value_dict: Dict[str, str] = {}
            print(f"convert_rppatsv {the_barcode}", flush=True, end='')
            counter: int = 0
            for line in io.TextIOWrapper(in_file, encoding="utf-8"):
                if 0 == counter % 1000:
                    print(".", flush=True, end='')
                counter += 1
                line = line.rstrip('\n')
                if not line.startswith('#'):
                    if headers is None:
                        # populate headers
                        headers = line.split("\t")
                    else:
                        # process star count data line
                        tsv_dict: Dict[str, str] = dict(zip(headers, line.split("\t")))
                        # uuid: str = tsv_dict['GDC_Aliquot']
                        id_a: str = tsv_dict['peptide_target']
                        id_b: str = tsv_dict['AGID']
                        feature: str = id_a + "|" + id_b
                        val_str: str = tsv_dict['protein_expression']
                        assert feature not in value_dict, f"convert_rppatsv - (1) found pre-existing feature {feature}"
                        value_dict[feature] = val_str
            print("-", flush=True, end='')
            if len(value_dict) > 0:
                my_sample_matrix: pandas.DataFrame = pandas.DataFrame.from_dict(value_dict, orient='index', dtype='str', columns=[the_barcode])
                ret_value = pandas.concat([the_matrix, my_sample_matrix], axis=1)
            print("-", flush=True)
    return ret_value
# pylint: enable=too-many-arguments,too-many-locals,consider-using-min-builtin,too-many-nested-blocks
