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


import io
import zipfile
from typing import Dict, List, Tuple
import pandas
from mbatch.gdcapi.convert_dataset import Dataset
from mbatch.gdcapi.download_datafile import GdcApiDatafile
from mbatch.gdcapi.convert_hg38_meth import Hg38Meth


# pylint: disable=too-many-arguments
def convert_sesamemeth(the_dataset: Dataset, the_meth450_map: Dict[str, Hg38Meth],
                       the_sample_dir: str, the_download_dir: str,
                       the_no_xy_flag: bool) -> Tuple[pandas.DataFrame, List[str]]:
    """
    Build the data matrix for this dataset and return it and the list of sample ids (barcodes)
    :param the_dataset: Dataset object describing matrix to build.
    :param the_meth450_map: Dictionary of Methylation Probes used in conversion.
    :param the_sample_dir: Full path to sample directory inside biospecimen directory.
    :param the_download_dir: Full path to download/data directory.
    :param the_no_xy_flag: if True, exclude sex-chromosomes.
    :return: A tuple of the DataFrame (matrix) and List of barcodes (sample ids) in the matrix.
    """
    if len(the_dataset.files) > 1000:
        print(f"SKIP!! convert_sesamemeth size={len(the_dataset.files)} for {the_dataset.get_dataset_path('/')}", flush=True)
        return pandas.DataFrame(), []
    print(f"convert_sesamemeth size={len(the_dataset.files)} for {the_dataset.get_dataset_path('/')}", flush=True)
    # my matrix to be populated
    my_matrix: pandas.DataFrame = pandas.DataFrame({}, dtype='str')
    sample_list: List[str] = []
    # my_batches: pandas.DataFrame = pandas.DataFrame({}, dtype='str')
    # my_clinical: pandas.DataFrame = pandas.DataFrame({}, dtype='str')
    # read headers from meth file (which is inside a ZIP archive)
    my_datafile: GdcApiDatafile
    for my_datafile in the_dataset.files.values():
        zip_file: str = my_datafile.get_file_archive(the_download_dir)
        # print(f"convert_sesamemeth zip_file={zip_file} file_name={my_datafile.file_name}", flush=True)
        my_datafile.read_sample_file(the_sample_dir)
        assert 1 == len(my_datafile.sample_dict), "convert_sesamemeth - only expected one sample"
        barcode: str = list(my_datafile.sample_dict.values())[0].sample_barcode
        my_matrix = read_and_process_file(my_matrix, barcode,
                                          zip_file, my_datafile.file_name,
                                          the_meth450_map, the_no_xy_flag)
        sample_list.append(barcode)
    # now we have a complete matrix my_matrix
    # and a list of samples in that matrix sample_list
    return my_matrix, sample_list
# pylint: enable=too-many-arguments


# pylint: disable=too-many-arguments,too-many-locals,consider-using-min-builtin,too-many-nested-blocks
def read_and_process_file(the_matrix: pandas.DataFrame, the_barcode: str,
                          the_file_zip: str, the_filename: str,
                          the_meth450_map: Dict[str, Hg38Meth], the_no_xy_flag: bool) -> pandas.DataFrame:
    """
    Build the data matrix for this dataset and return it and the list of sample ids (barcodes)
    :param the_matrix: pandas.DataFrame to which to add a new column
    :param the_barcode: sample id of new data.
    :param the_file_zip: Full path and file name of ZIP file with data file.
    :param the_filename: File name of data inside ZIP.
    :param the_meth450_map: Dictionary of Methylation Probes used in conversion.
    :param the_no_xy_flag: if True, exclude sex-chromosomes.
    :return: Return the_matrix with new column (sample) of data added.
    """
    ret_value: pandas.DataFrame = the_matrix
    zip_file: zipfile.ZipFile
    in_file: io.TextIOWrapper
    with zipfile.ZipFile(the_file_zip, 'r') as zip_file:
        with zip_file.open(the_filename, mode="r") as in_file:
            print(f"convert_sesamemeth {the_barcode}", flush=True, end='')
            value_dict: Dict[str, str] = {}
            counter: int = 0
            for line in io.TextIOWrapper(in_file, encoding="utf-8"):
                # THIS FILE HAS NO HEADERS PROBE\tVALUE
                if 0 == counter % 10000:
                    print(".", flush=True, end='')
                counter += 1
                line = line.rstrip('\n')
                # process sesame meth data line
                line_splitted: List[str] = line.split("\t")
                probe: str = line_splitted[0]
                val_str: str = line_splitted[1]
                use_probe: bool = False
                hg38_probe: Hg38Meth = the_meth450_map.get(probe)
                # the check (hg38_probe is not None) filters mitochondria probes,
                # which are present in the data file but not in the HG38 methylation probe file.
                # Will also filter unrecognized probes, regardless of XY flag.
                if hg38_probe is not None:
                    if the_no_xy_flag:
                        if not "y" == hg38_probe.chromosome.lower():
                            if not "x" == hg38_probe.chromosome.lower():
                                use_probe = True
                    else:
                        use_probe = True
                if use_probe:
                    assert probe not in value_dict, f"convert_sesamemeth - (1) found pre-existing probe {probe}"
                    value_dict[probe] = val_str
            print("-", flush=True, end='')
            if len(value_dict) > 0:
                my_sample_matrix: pandas.DataFrame = pandas.DataFrame.from_dict(value_dict, orient='index', dtype='str', columns=[the_barcode])
                ret_value = pandas.concat([the_matrix, my_sample_matrix], axis=1)
            print("-", flush=True)
    return ret_value
# pylint: enable=too-many-arguments,too-many-locals,consider-using-min-builtin,too-many-nested-blocks
