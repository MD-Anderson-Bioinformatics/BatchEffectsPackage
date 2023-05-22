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
import re
import zipfile
import typing
from typing import Dict, List, Tuple
import pandas
from mbatch.gdcapi.convert_dataset import Dataset
from mbatch.gdcapi.download_datafile import GdcApiDatafile


# pylint: disable=too-many-arguments
def convert_mutationmaf(the_dataset: Dataset, the_sample_dir: str, the_download_dir: str) -> \
        Tuple[pandas.DataFrame, List[str], pandas.DataFrame, List[str]]:
    """
    Build the TWO data matrixes for this dataset and return them and their lists of sample ids (barcodes)
    One matrix is non-silent mutation count. One is rationalized, simplified MAF-like file of all mutation calls.
    :param the_dataset: Dataset object describing matrix to build.
    :param the_sample_dir: Full path to sample directory inside biospecimen directory.
    :param the_download_dir: Full path to download/data directory.
    :return: A tuple of the DataFrame (matrix) and List of barcodes (sample ids) in the matrix and the MAF rewrite DataFrame and list of barcodes.
    """
    print(f"convert_mutationmaf size={len(the_dataset.files)} for {the_dataset.get_dataset_path('/')}", flush=True)
    # my matrix to be populated (use int as it counts non-silent mutations)
    my_matrix: pandas.DataFrame = pandas.DataFrame({}, dtype='int')
    # my dataframe to be populated (nicer version of MAF with mutation information)
    my_dataframe: pandas.DataFrame = pandas.DataFrame({}, dtype='str')
    # list of samples used
    sample_list_matrix: List[str] = []
    sample_list_dataframe: List[str] = []
    my_datafile: GdcApiDatafile
    for my_datafile in the_dataset.files.values():
        zip_file: str = my_datafile.get_file_archive(the_download_dir)
        # read sample information (may not be needed)
        my_datafile.read_sample_file(the_sample_dir)
        # do not need to pass barcode, is included in MAF file
        # remove .gz from end of filename
        filename: str = my_datafile.file_name.removesuffix(".gz")
        barcode: str
        my_matrix, barcode = read_and_process_file_matrix(my_matrix, zip_file, filename)
        if barcode is not None:
            sample_list_matrix.append(barcode)
        my_dataframe, barcode = read_and_process_file_dataframe(my_dataframe, zip_file, filename)
        if barcode is not None:
            sample_list_dataframe.append(barcode)
    # now we have a complete matrix my_matrix
    # and a list of samples in that matrix sample_list
    return my_matrix, sample_list_matrix, my_dataframe, sample_list_dataframe
# pylint: enable=too-many-arguments


# pylint: disable=too-many-arguments,too-many-locals,consider-using-min-builtin,too-many-nested-blocks
def read_and_process_file_matrix(the_matrix: pandas.DataFrame, the_file_zip: str, the_filename: str) -> Tuple[pandas.DataFrame, str]:
    """
    Read a MAF file inside a ZIP archive and add column of data to the matrix.
    In this case the_matrix is a matrix af non-silent mutation counts.
    We return the list, since unlike other converters, MAF files may have multiple samples included,
    but are for a single patient.
    :param the_matrix: pandas.DataFrame to which to add a new column
    :param the_file_zip: Full path and file name of ZIP file with data file.
    :param the_filename: File name of data inside ZIP.
    :return: A tuple with the dataframe a list of samples in the dataframe
    """
    ret_value: pandas.DataFrame = the_matrix
    barcode_check: typing.Optional[str] = None
    zip_file: zipfile.ZipFile
    in_file: io.TextIOWrapper
    headers: typing.Optional[List[str]] = None
    with zipfile.ZipFile(the_file_zip, 'r') as zip_file:
        with zip_file.open(the_filename, mode="r") as in_file:
            value_dict: Dict[str, int] = {}
            counter: int = 0
            for line in io.TextIOWrapper(in_file, encoding="utf-8"):
                line = line.rstrip('\n')
                # no known comments but just in case...
                if not line.startswith('#'):
                    if headers is None:
                        # populate headers
                        headers = line.split("\t")
                    else:
                        # process mutation data line
                        tsv_dict: Dict[str, str] = dict(zip(headers, line.split("\t")))
                        barcode: str = tsv_dict['Tumor_Sample_Barcode']
                        if barcode_check is None:
                            barcode_check = barcode
                            print(f"convert_mutationmaf matrix {barcode_check}", flush=True, end='')
                        assert barcode_check == barcode, f"Barcodes do not match. {barcode_check}!={barcode} File should only have one barcode."
                        if 0 == counter % 1000:
                            print(".", flush=True, end='')
                        counter += 1
                        build: str = tsv_dict['NCBI_Build']
                        if "38" in build:
                            gene: str = tsv_dict['Hugo_Symbol']
                            variant: str = tsv_dict['Variant_Classification']
                            if not "silent" == variant.lower():
                                if gene in value_dict:
                                    value_dict[gene] = value_dict[gene] + 1
                                else:
                                    value_dict[gene] = 1
            print("-", flush=True, end='')
            if len(value_dict) > 0:
                assert barcode_check is not None, "Barcodes should not be empty."
                my_sample_matrix: pandas.DataFrame = pandas.DataFrame.from_dict(value_dict, orient='index', dtype='int', columns=[barcode_check])
                ret_value = pandas.concat([ret_value, my_sample_matrix], axis=1)
            print("-", flush=True)
    return ret_value, barcode_check
# pylint: enable=too-many-arguments,too-many-locals,consider-using-min-builtin,too-many-nested-blocks


# pylint: disable=too-many-arguments,too-many-locals,consider-using-min-builtin,too-many-nested-blocks,too-many-branches,too-many-statements
def read_and_process_file_dataframe(the_matrix: pandas.DataFrame, the_file_zip: str, the_filename: str) -> Tuple[pandas.DataFrame, str]:
    """
    Read a MAF file inside a ZIP archive and add column of data to the matrix.
    In this case the_matrix is an actual dataframe of various mutations.
    The dataframe has column names rationalized/simplified.
    (Some of those changes are legacy from older versions of data with different column ids.)
    We return the list, since unlike other converters, MAF files may have multiple samples included,
    but are for a single patient.
    :param the_matrix: pandas.DataFrame to which to add a new column
    :param the_file_zip: Full path and file name of ZIP file with data file.
    :param the_filename: File name of data inside ZIP.
    :return: A tuple with the dataframe a list of samples in the dataframe
    """
    ret_value: pandas.DataFrame = the_matrix
    barcode_check: typing.Optional[str] = None
    zip_file: zipfile.ZipFile
    in_file: io.TextIOWrapper
    headers: typing.Optional[List[str]] = None
    with zipfile.ZipFile(the_file_zip, 'r') as zip_file:
        with zip_file.open(the_filename, mode="r") as in_file:
            value_dict: Dict[str, str] = {}
            counter: int = 0
            for line in io.TextIOWrapper(in_file, encoding="utf-8"):
                line = line.rstrip('\n')
                # no known comments but just in case...
                if not line.startswith('#'):
                    if headers is None:
                        # populate headers
                        headers = line.split("\t")
                    else:
                        # process mutation data line
                        tsv_dict: Dict[str, str] = dict(zip(headers, line.split("\t")))
                        barcode: str = tsv_dict['Tumor_Sample_Barcode']
                        if barcode_check is None:
                            barcode_check = barcode
                            print(f"convert_mutationmaf dataframe {barcode_check}", flush=True, end='')
                        assert barcode_check == barcode, f"Barcodes do not match. {barcode_check}!={barcode} File should only have one barcode."
                        if 0 == counter % 1000:
                            print(".", flush=True, end='')
                        counter += 1
                        build: str = tsv_dict['NCBI_Build']
                        if "38" in build:
                            # "Tumor_Sample_Barcode"
                            value_dict["Tumor_Sample_Barcode"] = barcode
                            # "Gene"
                            value_dict["Tumor_Sample_Barcode"] = tsv_dict['Hugo_Symbol']
                            # "EntrezId"
                            value_dict["EntrezId"] = tsv_dict['Entrez_Gene_Id']
                            # "Variant_Classification"
                            value_dict["Variant_Classification"] = tsv_dict['Variant_Classification']
                            # "Position_Chromosome"
                            value_dict["Position_Chromosome"] = tsv_dict['Chromosome']
                            # "Position_Start"
                            value_dict["Position_Start"] = tsv_dict['Start_Position']
                            # "Position_End"
                            value_dict["Position_End"] = tsv_dict['End_Position']
                            # "Position_Strand"
                            value_dict["Position_Strand"] = tsv_dict['Strand']
                            # "TranscriptId"
                            value_dict["TranscriptId"] = tsv_dict['Transcript_ID']
                            # "Tumor_Depth"
                            value_dict["Tumor_Depth"] = tsv_dict['t_depth']
                            # "Tumor_Reference_Count"
                            value_dict["Tumor_Reference_Count"] = tsv_dict['t_ref_count']
                            # "Tumor_Variant_Count"
                            value_dict["Tumor_Variant_Count"] = tsv_dict['t_alt_count']
                            # "Normal_Depth"
                            value_dict["Normal_Depth"] = tsv_dict['n_depth']
                            # "Normal_Reference_Count"
                            value_dict["Normal_Variant_Count"] = tsv_dict['n_ref_count']
                            # "Normal_Variant_Count"
                            value_dict["Normal_Variant_Count"] = tsv_dict['n_alt_count']
                            # "HGVSp_Short"
                            value_dict["HGVSp_Short"] = tsv_dict['HGVSp_Short']
                            # "amino_acid_position"
                            amino_acid_position: str = tsv_dict['HGVSp_Short']
                            if '' != amino_acid_position:
                                amino_acid_position = amino_acid_position.removeprefix("p.")
                                matches: List[str] = re.findall("(\\d+)", amino_acid_position)
                                amino_acid_position = matches[0]
                            value_dict["amino_acid_position"] = amino_acid_position
                            # "amino_acid_normal"
                            amino_acid_normal: str = tsv_dict['HGVSp_Short']
                            if '' != amino_acid_normal:
                                amino_acid_normal = amino_acid_normal.removeprefix("p.")
                                if ">" in amino_acid_normal:
                                    matches: List[str] = re.findall("(\\D+)>", amino_acid_normal)
                                    amino_acid_normal = matches[0]
                                else:
                                    matches: List[str] = re.findall("(\\D+)(\\d+)", amino_acid_normal)
                                    amino_acid_normal = matches[0][0]
                            value_dict["amino_acid_normal"] = amino_acid_normal
                            # "amino_acid_tumor"
                            amino_acid_tumor: str = tsv_dict['HGVSp_Short']
                            if '' != amino_acid_tumor:
                                amino_acid_tumor = amino_acid_tumor.removeprefix("p.")
                                if ">" in amino_acid_tumor:
                                    matches: List[str] = re.findall(">(\\D+)$", amino_acid_tumor)
                                    amino_acid_tumor = matches[0]
                                else:
                                    # can also have an asterisk
                                    matches: List[str] = re.findall("(?:[A-Za-z*]+)(?:\\d+)(?:_?[\\d]+)?(.+)?$", amino_acid_tumor)
                                    amino_acid_tumor = matches[0]
                            value_dict["amino_acid_tumor"] = amino_acid_tumor
                            # "NCBI_Build"
                            value_dict["NCBI_Build"] = build
            print("-", flush=True, end='')
            if len(value_dict) > 0:
                assert barcode_check is not None, "Barcodes should not be empty."
                my_sample_matrix: pandas.DataFrame = pandas.DataFrame.from_dict(value_dict, orient='index', dtype='str', columns=[barcode_check])
                ret_value = pandas.concat([ret_value, my_sample_matrix], axis=1)
            print("-", flush=True)
    return ret_value, barcode_check
# pylint: enable=too-many-arguments,too-many-locals,consider-using-min-builtin,too-many-nested-blocks,too-many-branches,too-many-statements
