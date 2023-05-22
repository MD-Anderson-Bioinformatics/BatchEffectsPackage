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


from typing import Dict, List
import os
import typing
import gc
import pandas
from mbatch.gdcapi.convert_dataset import Dataset
from mbatch.gdcapi.standardized_data import write_converted_data, write_converted_dataframe
from mbatch.gdcapi.convert_hg38_gene import Hg38Gene, read_hg38_gene_file
from mbatch.gdcapi.convert_hg38_meth import Hg38Meth, read_hg38_meth_file
from mbatch.gdcapi.converter_snp6txt import convert_snp6txt
from mbatch.gdcapi.converter_starcounts import convert_starcounts
from mbatch.gdcapi.converter_sesamemeth import convert_sesamemeth
from mbatch.gdcapi.converter_rppatsv import convert_rppatsv
from mbatch.gdcapi.converter_mirnatxt import convert_mirnatxt
from mbatch.gdcapi.converter_copynumtsv import convert_copynumtsv
from mbatch.gdcapi.converter_mutationmaf import convert_mutationmaf
from mbatch.gdcapi.converter_scdiff import convert_scdiff
from mbatch.test.common import add_error, delete_directory_contents, extract_zip_to_dir, archive_dir_contents_to_zip, add_warnings


#
# CONVERT DATASET FILES INTO MATRIX FILES
#

# defined here to prevent circular references
# pylint: disable=too-many-arguments
def convert_datasets(the_datasets: Dict[str, Dataset],
                     the_util_dir: str, the_sample_dir: str,
                     the_download_dir: str, the_biospecimen_dir: str, the_clinical_dir: str,
                     the_convert_dir: str, the_temp_dir: str) -> None:
    """
    Take list of datasets that need conversion, and convert them
    :param the_datasets:  Dictionary of Dataset objects
    :param the_util_dir: full directory path to util dir, that has gene and probe maps
    :param the_sample_dir: full directory path to sample dir inside biospecimen
    :param the_download_dir: full directory path to download/data dir
    :param the_biospecimen_dir: full directory path to converted/biospecimen dir
    :param the_clinical_dir: full directory path to converted/clinical dir
    :param the_convert_dir: full directory path to converted/data dir
    :param the_temp_dir: full path to temp directory used to build uncompressed files before ZIPing them
    :return: Nothing
    """
    # the_biospecimens: Dict[str, GdcApiBiospecimen], the_clinicals: Dict[str, GdcApiClinical],
    print("convert_datasets start", flush=True)
    # read HG38_Genes file into a Dict with ensembl-id as the key.
    # note that ensembl-id is spelled ensembl_id for compatibility with Python syntax
    hg38_genes: Dict[str, Hg38Gene] = read_hg38_gene_file(os.path.join(the_util_dir, "HG38_Genes.tsv"), 'ensembl_id')
    hg38_meth: Dict[str, Hg38Meth] = read_hg38_meth_file(os.path.join(the_util_dir, "HG38_SSM450.tsv"))
    my_dataset: Dataset
    for my_dataset in the_datasets.values():
        # my_biospecimen: GdcApiBiospecimen = my_dataset.get_biospecimen(the_biospecimens)
        # my_clinical: GdcApiClinical = my_dataset.get_clinical(the_clinicals)
        convert_case(my_dataset, hg38_genes, hg38_meth, the_sample_dir, the_download_dir,
                     the_biospecimen_dir, the_clinical_dir, the_convert_dir, the_temp_dir)
    print("convert_datasets done", flush=True)
# pylint: enable=too-many-arguments


# pylint: disable=too-many-arguments,too-many-branches,too-many-statements
def convert_case(the_dataset: Dataset, the_hg38_map: Dict[str, Hg38Gene], the_hg38_meth: Dict[str, Hg38Meth],
                 the_sample_dir: str, the_download_dir: str,
                 the_biospecimen_dir: str, the_clinical_dir: str, the_convert_dir: str, the_temp_dir: str) -> None:
    """
    Determine convert type and perform conversion.
    Strings must match those from convert_dataset::build_datasets.
    TODO: find better way to encapsulate convert matching strings
    :param the_dataset:  the Dataset object to convert
    :param the_hg38_map: HG38 Gene Maps for simple gene scoring algorithm
    :param the_hg38_meth: HG38 Meth probe data for SeSAME mapping
    :param the_sample_dir: full directory path to sample dir inside biospecimen
    :param the_download_dir: full directory path to download/data dir
    :param the_biospecimen_dir: full directory path to converted/biospecimen dir
    :param the_clinical_dir: full directory path to converted/clinical dir
    :param the_convert_dir: full directory path to converted/data dir
    :param the_temp_dir: full path to temp directory used to build uncompressed files before ZIPing them
    :return: nothing
    """
    # determine what kind of convert to use
    my_matrix: pandas.DataFrame
    sample_list: List[str]
    if the_dataset.conversion_match("Genotyping Array", "Copy Number Segment", "DNAcopy", "Copy-Number-with-CNV"):
        # Snp6TXT processing (historical name, based on first dataset conversion worked on)
        my_matrix, sample_list = convert_snp6txt(the_dataset, the_hg38_map, the_sample_dir, the_download_dir)
        write_converted_data_case(sample_list, my_matrix, the_dataset, the_biospecimen_dir, the_clinical_dir, the_convert_dir, the_temp_dir, None)
    elif the_dataset.conversion_match("Genotyping Array", "Masked Copy Number Segment", "DNAcopy", "Copy-Number-no-CNV"):
        # Snp6TXT processing (historical name, based on first dataset conversion worked on)
        my_matrix, sample_list = convert_snp6txt(the_dataset, the_hg38_map, the_sample_dir, the_download_dir)
        write_converted_data_case(sample_list, my_matrix, the_dataset, the_biospecimen_dir, the_clinical_dir, the_convert_dir, the_temp_dir, None)
    elif the_dataset.conversion_match("RNA-Seq", "Gene Expression Quantification", "STAR - Counts", "RNASeq-TPM"):
        # StarCounts processing
        my_matrix, sample_list = convert_starcounts(the_dataset, the_sample_dir, the_download_dir, "tpm_unstranded")
        write_converted_data_case(sample_list, my_matrix, the_dataset, the_biospecimen_dir, the_clinical_dir, the_convert_dir, the_temp_dir, None)
    elif the_dataset.conversion_match("RNA-Seq", "Gene Expression Quantification", "STAR - Counts", "RNASeq-FPKM"):
        # StarCounts processing
        my_matrix, sample_list = convert_starcounts(the_dataset, the_sample_dir, the_download_dir, "fpkm_unstranded")
        write_converted_data_case(sample_list, my_matrix, the_dataset, the_biospecimen_dir, the_clinical_dir, the_convert_dir, the_temp_dir, None)
    elif the_dataset.conversion_match("RNA-Seq", "Gene Expression Quantification", "STAR - Counts", "RNASeq-FPKM-UQ"):
        # StarCounts processing
        my_matrix, sample_list = convert_starcounts(the_dataset, the_sample_dir, the_download_dir, "fpkm_uq_unstranded")
        write_converted_data_case(sample_list, my_matrix, the_dataset, the_biospecimen_dir, the_clinical_dir, the_convert_dir, the_temp_dir, None)
    elif the_dataset.conversion_match("Methylation Array", "Methylation Beta Value", "SeSAMe Methylation Beta Estimation", "Methylation-With-Sex-Chromosomes"):
        # Sesame Methylation processing
        my_matrix, sample_list = convert_sesamemeth(the_dataset, the_hg38_meth, the_sample_dir, the_download_dir, False)
        write_converted_data_case(sample_list, my_matrix, the_dataset, the_biospecimen_dir, the_clinical_dir, the_convert_dir, the_temp_dir, None)
    elif the_dataset.conversion_match("Methylation Array", "Methylation Beta Value", "SeSAMe Methylation Beta Estimation", "Methylation-No-Sex-Chromosomes"):
        # Sesame Methylation processing
        my_matrix, sample_list = convert_sesamemeth(the_dataset, the_hg38_meth, the_sample_dir, the_download_dir, True)
        if len(sample_list) > 0:
            write_converted_data_case(sample_list, my_matrix, the_dataset, the_biospecimen_dir, the_clinical_dir, the_convert_dir, the_temp_dir, None)
        else:
            add_warnings(f"Methylation data empty for {the_dataset.get_dataset_path('/')}")
    elif the_dataset.conversion_match("Reverse Phase Protein Array", "Protein Expression Quantification", "Protein Analysis", "RPPA"):
        # RPPA processing
        my_matrix, sample_list = convert_rppatsv(the_dataset, the_sample_dir, the_download_dir)
        write_converted_data_case(sample_list, my_matrix, the_dataset, the_biospecimen_dir, the_clinical_dir, the_convert_dir, the_temp_dir, None)
    elif the_dataset.conversion_match("miRNA-Seq", "Isoform Expression Quantification", "BCGSC miRNA Profiling", "miRNA-Isoform"):
        # miRNA processing
        my_matrix, sample_list = convert_mirnatxt(the_dataset, the_sample_dir, the_download_dir)
        write_converted_data_case(sample_list, my_matrix, the_dataset, the_biospecimen_dir, the_clinical_dir, the_convert_dir, the_temp_dir, None)
    elif the_dataset.conversion_match("miRNA-Seq", "miRNA Expression Quantification", "BCGSC miRNA Profiling", "miRNA-Genes"):
        # miRNA processing
        my_matrix, sample_list = convert_mirnatxt(the_dataset, the_sample_dir, the_download_dir)
        write_converted_data_case(sample_list, my_matrix, the_dataset, the_biospecimen_dir, the_clinical_dir, the_convert_dir, the_temp_dir, None)
    elif the_dataset.conversion_match("Genotyping Array", "Gene Level Copy Number", "ASCAT2", "Copy-Number-With-Sex-Chromosomes"):
        # Copy Number processing
        my_matrix, sample_list = convert_copynumtsv(the_dataset, the_sample_dir, the_download_dir, False, the_biospecimen_dir)
        write_converted_data_case(sample_list, my_matrix, the_dataset, the_biospecimen_dir, the_clinical_dir, the_convert_dir, the_temp_dir, None)
    elif the_dataset.conversion_match("Genotyping Array", "Gene Level Copy Number", "ASCAT2", "Copy-Number-No-Sex-Chromosomes"):
        # Copy Number processing
        my_matrix, sample_list = convert_copynumtsv(the_dataset, the_sample_dir, the_download_dir, True, the_biospecimen_dir)
        write_converted_data_case(sample_list, my_matrix, the_dataset, the_biospecimen_dir, the_clinical_dir, the_convert_dir, the_temp_dir, None)
    elif the_dataset.conversion_match("WGS", "Gene Level Copy Number", "AscatNGS", "Copy-Number-With-Sex-Chromosomes"):
        # Copy Number processing
        my_matrix, sample_list = convert_copynumtsv(the_dataset, the_sample_dir, the_download_dir, False, the_biospecimen_dir)
        write_converted_data_case(sample_list, my_matrix, the_dataset, the_biospecimen_dir, the_clinical_dir, the_convert_dir, the_temp_dir, None)
    elif the_dataset.conversion_match("WGS", "Gene Level Copy Number", "AscatNGS", "Copy-Number-No-Sex-Chromosomes"):
        # Copy Number processing
        my_matrix, sample_list = convert_copynumtsv(the_dataset, the_sample_dir, the_download_dir, True, the_biospecimen_dir)
        write_converted_data_case(sample_list, my_matrix, the_dataset, the_biospecimen_dir, the_clinical_dir, the_convert_dir, the_temp_dir, None)
    elif "Mutation-Calling" == the_dataset.details:
        # mutations processing
        my_dataframe: pandas.DataFrame
        sample_list_dataframe: List[str]
        my_matrix, sample_list, my_dataframe, sample_list_dataframe = convert_mutationmaf(the_dataset, the_sample_dir, the_download_dir)
        # not really used
        sample_list_dataframe.sort()
        write_converted_data_case(sample_list, my_matrix, the_dataset, the_biospecimen_dir, the_clinical_dir, the_convert_dir, the_temp_dir, my_dataframe)
    elif the_dataset.conversion_match("scRNA-Seq", "Differential Gene Expression", "Seurat - 10x Chromium", "scRNA-Differential"):
        # sc RNA processing
        my_matrix, sample_list = convert_scdiff(the_dataset, the_sample_dir, the_download_dir)
        write_converted_data_case(sample_list, my_matrix, the_dataset, the_biospecimen_dir, the_clinical_dir, the_convert_dir, the_temp_dir, None)
    else:
        add_error(f"Unrecognized conversion dataset={the_dataset.get_dataset_path('')}")
    gc.collect()
# pylint: enable=too-many-arguments,too-many-branches,too-many-statements


# pylint: disable=too-many-arguments
def write_converted_data_case(the_sample_list: List[str], the_matrix: pandas.DataFrame, the_dataset: Dataset,
                              the_biospecimen_dir: str, the_clinical_dir: str, the_convert_dir: str, the_temp_dir: str,
                              the_mutations_dataframe: typing.Optional[pandas.DataFrame]) -> None:
    """
    After conversion, write the converted matrix and batch information.
    Matrix is passed in, but dataset is used to read program/project level
    biospecimen information and filter into just this batch data.
    :param the_sample_list: list of samples in this file
    :param the_matrix: matrix of data to write. (Columns are samples. Features are rows.)
    :param the_dataset: Dataset object, used to read batch data
    :param the_biospecimen_dir: full directory path to converted/biospecimen dir
    :param the_clinical_dir: full directory path to converted/clinical dir
    :param the_convert_dir: full directory path to converted/data dir
    :param the_temp_dir: full path to temp directory used to build uncompressed files before ZIPing them
    :param the_mutations_dataframe: dataframe of mutations data. write if not none
    :return: nothing
    """
    # call write for standardized data files
    # updated to support working in temp unZIPped directory the_temp_dir
    if len(the_sample_list) > 0:
        # empty existing contents of temp dir, if any
        print(f"Delete contents of {the_temp_dir}", flush=True)
        delete_directory_contents(the_temp_dir)
        # paths to places in temp directory structure
        versions_dir: str = os.path.join(the_temp_dir, 'versions')
        versions_data_dir: str = os.path.join(versions_dir, f"DATA_{the_dataset.version}")
        # extract existing ZIP archive (if any)
        archive_path: str = the_dataset.get_archive_path(the_convert_dir)
        archive_file: str = the_dataset.get_archive_filename()
        full_archive: str = os.path.join(archive_path, archive_file)
        if os.path.exists(full_archive):
            # extract archive to temp dir, if archive exists
            print(f"Extract {full_archive}", flush=True)
            extract_zip_to_dir(full_archive, the_temp_dir)
        else:
            # if archive does not exist, create the infrastructure in the temp dir
            print(f"Initialize {versions_dir}", flush=True)
            os.makedirs(versions_dir)
            # write dataset info
            the_dataset.dump_to_dict_file(os.path.join(the_temp_dir, 'index.json'))
        # new directories for new data
        os.makedirs(versions_data_dir)
        # write versions/DATA_<GDC-History-Date>/original/batches.tsv and matrix.tsv
        program_batch_file: typing.Optional[str] = the_dataset.get_batches_tsv_path(the_biospecimen_dir)
        program_clinical_file: typing.Optional[str] = the_dataset.get_clinical_tsv_path(the_clinical_dir)
        # This means versions/DATA_<GDC-History-Date>/original/batches.tsv and matrix.tsv
        print(f"Write Data To {versions_data_dir}", flush=True)
        write_converted_data(the_matrix, the_sample_list, program_batch_file, program_clinical_file, versions_data_dir)
        if the_mutations_dataframe is not None:
            print(f"Write Mutations To {versions_data_dir}", flush=True)
            write_converted_dataframe(the_mutations_dataframe, versions_data_dir, "mutations.tsv", None)
        # historical, copy contents of directory A to B
        # shutil.copytree(dir_source_a, dir_dest_b, dirs_exist_ok=True)
        # update to compress into and replace/copy-to ZIP archive
        print(f"Archive {full_archive}", flush=True)
        if os.path.exists(full_archive):
            os.unlink(full_archive)
        archive_dir_contents_to_zip(full_archive, the_temp_dir)
    else:
        add_warnings(f"No samples found for conversion of dataset {the_dataset.get_dataset_path('/')}")
# pylint: enable=too-many-arguments
