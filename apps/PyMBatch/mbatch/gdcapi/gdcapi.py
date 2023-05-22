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


import os.path
from typing import List, Dict, Tuple
from mbatch.gdcapi.download_history_mixin import update_history_from_gdc
from mbatch.gdcapi.download_biospecimen import GdcApiBiospecimen, read_biospecimen_index, get_biospecimen_list_from_gdc, write_biospecimen_index
from mbatch.gdcapi.download_clinical import GdcApiClinical, read_clinical_index, get_clinical_list_from_gdc, write_clinical_index
from mbatch.gdcapi.download_datafile import GdcApiDatafile, read_datafile_index, get_datafile_list_from_gdc
from mbatch.gdcapi.download_datafile import write_datafile_index, read_map_datafiles_by_program_projects
from mbatch.gdcapi.convert_dataset import Dataset, build_datasets, get_unconverted_datasets
from mbatch.gdcapi.convert_case import convert_datasets
from mbatch.gdcapi.convert_aliquot import convert_biospecimen_batches
from mbatch.gdcapi.convert_clinical import convert_clinical_batches
from mbatch.test.common import add_error


def update_clinical_index(the_clinical_index: str) -> None:
    """
    Update list of known clinical files, and their history (release) status.
    :param the_clinical_index: full path to clinical index file, with file name.
    :return: None
    """
    print("*************START:update_clinical_index*************", flush=True)
    # read existing index files (if any)
    my_entries: Dict[str, GdcApiClinical] = {}
    if os.path.exists(the_clinical_index):
        my_entries = read_clinical_index(the_clinical_index)
    # update list of files, adding any new ones (new ones will not have history data)
    fresh_entries: List[GdcApiClinical] = get_clinical_list_from_gdc(my_entries)
    # save the current index files (if needed for restart)
    write_clinical_index(the_clinical_index, fresh_entries)
    # update history
    if update_history_from_gdc(fresh_entries, 5):
        # save index files if anything changed
        write_clinical_index(the_clinical_index, fresh_entries)
    print("*************DONE:update_clinical_index*************", flush=True)


def download_clinical_index(the_clinical_index: str, the_base_dir: str) -> None:
    """
    Download un-downloaded clinical files into appropriate base dir and file-specific path.
    Resulting files are un-gz'd if needed, and are then ZIPped.
    :param the_clinical_index: full path to clinical index file.
    :param the_base_dir: full path to base directory.
    :return: None
    """
    print("*************START:download_clinical_index*************", flush=True)
    if os.path.exists(the_clinical_index):
        my_entries: Dict[str, GdcApiClinical] = read_clinical_index(the_clinical_index)
        entry: GdcApiClinical
        for entry in my_entries.values():
            entry.download_gdcapi_file(the_base_dir)
    print("*************DONE:download_clinical_index*************", flush=True)


def update_biospecimen_index(the_biospecimen_index: str) -> None:
    """
    Update list of known biospecimen files, and their history (release) status.
    :param the_biospecimen_index: full path to biospecimen index file, with file name.
    :return: None
    """
    print("*************START:update_biospecimen_index*************", flush=True)
    # read existing index files (if any)
    my_entries: Dict[str, GdcApiBiospecimen] = {}
    if os.path.exists(the_biospecimen_index):
        my_entries = read_biospecimen_index(the_biospecimen_index)
    # update list of files, adding any new ones (new ones will not have history data)
    fresh_entries: List[GdcApiBiospecimen] = get_biospecimen_list_from_gdc(my_entries)
    # save the current index files (if needed for restart)
    write_biospecimen_index(the_biospecimen_index, fresh_entries)
    # update history
    if update_history_from_gdc(fresh_entries, 5):
        # save index files if anything changed
        write_biospecimen_index(the_biospecimen_index, fresh_entries)
    print("*************DONE:update_biospecimen_index*************", flush=True)


def download_biospecimen_index(the_biospecimen_index: str, the_base_dir: str) -> None:
    """
    Download un-downloaded biospecimen files into appropriate base dir and file-specific path.
    Resulting files are un-gz'd if needed, and are then ZIPped.
    :param the_biospecimen_index: full path to clinical index file.
    :param the_base_dir: full path to base directory.
    :return: None
    """
    print("*************START:download_biospecimen_index*************", flush=True)
    if os.path.exists(the_biospecimen_index):
        my_entries: Dict[str, GdcApiBiospecimen] = read_biospecimen_index(the_biospecimen_index)
        entry: GdcApiBiospecimen
        for entry in my_entries.values():
            entry.download_gdcapi_file(the_base_dir)
    print("*************DONE:download_biospecimen_index*************", flush=True)


def update_datafile_index(the_file_datafiles: str, the_sample_dir: str) -> None:
    """
    Update list of known datafile files, and their history (release) status.
    :param the_file_datafiles: full path to datafile index file, with file name.
    :param the_sample_dir: full path to directory for index files
    :return: None
    """
    print("*************START:update_datafile_index*************", flush=True)
    # read existing index files (if any)
    my_entries: Dict[str, GdcApiDatafile] = {}
    if os.path.exists(the_file_datafiles):
        my_entries = read_datafile_index(the_file_datafiles)
    # update list of files, adding any new one (new ones will not have history data)
    all_entries: List[GdcApiDatafile] = get_datafile_list_from_gdc(my_entries)
    # save the current index files (if needed for restart)
    write_datafile_index(the_file_datafiles, the_sample_dir, all_entries)
    # update history
    if update_history_from_gdc(all_entries, 5):
        # save index files if anything changed
        write_datafile_index(the_file_datafiles, the_sample_dir, all_entries)
    print("*************DONE:update_datafile_index*************", flush=True)


def download_datafile_index(the_datafile_index: str, the_download_dir: str, the_sample_dir: str) -> None:
    """
    Download un-downloaded datafile files into appropriate base dir and file-specific path.
    Resulting files are un-gz'd if needed, and are then ZIPped.
    :param the_datafile_index: full path to datafile index file.
    :param the_download_dir: full path to base directory.
    :param the_sample_dir: full path to directory with sample index files.
    :return: None
    """
    print("download_datafile_index start", flush=True)
    print(f"download_datafile_index the_datafile_index={the_datafile_index}", flush=True)
    print(f"download_datafile_index the_sample_dir={the_sample_dir}", flush=True)
    print(f"download_datafile_index the_download_dir={the_download_dir}", flush=True)
    # read known files/samples
    my_entries: Dict[str, GdcApiDatafile] = read_datafile_index(the_datafile_index)
    # failed downloads
    failed: List[GdcApiDatafile] = []
    # iterate over known files
    entry: GdcApiDatafile
    index: int = 1
    for entry in my_entries.values():
        if 0 == index % 10:
            print(f"download_datafile_index {index} / {len(my_entries.values())}", flush=True)
        index += 1
        # get list of available files to download
        # download: str = entry.get_download_file(the_download_dir)
        if not entry.download_gdcapi_file(the_download_dir):
            failed.append(entry)
    for entry in failed:
        add_error(f"download_datafile_index failed={entry.uuid}")
    print("download_datafile_index done", flush=True)


# pylint: disable=too-many-arguments
def convert_update_datasets(the_datafile_index: str, the_biospecimen_index: str, the_clinical_index: str,
                            the_download_dir: str, the_biospecimen_dir: str, the_clinical_dir: str,
                            the_convert_dir: str, the_sample_dir: str, the_util_dir: str, the_temp_dir: str) -> None:
    """
    Perform conversion on the dataset files to create matrix output
    :param the_datafile_index: full path and filename to index file of dataset files
    :param the_biospecimen_index: biospecimen index (check if needed for conversion)
    :param the_clinical_index: full path and filename to index file of clinical files (check if needed for conversion)
    :param the_download_dir: full path to download directory for dataset files
    :param the_biospecimen_dir: full path to biospecimen directory (where batch files are stored) (check if needed for conversion)
    :param the_clinical_dir: full path to clinical directory (where batch files are stored) (check if needed for conversion)
    :param the_convert_dir: full path to convert directory, where matrix files go
    :param the_sample_dir: full path to sample dir, which are associated with the dataset index (REQUIRED for some conversions)
    :param the_util_dir: full path to util dir, containing info like gene maps (REQUIRED for some conversions)
    :param the_temp_dir: full path to temp directory used to build uncompressed files before ZIPing them
    :return: Nothing
    """
    print("*************START:convert_update_datasets*************", flush=True)
    print(f"convert_update_datasets the_datafile_index={the_datafile_index}", flush=True)
    print(f"convert_update_datasets the_biospecimen_index={the_biospecimen_index}", flush=True)
    print(f"convert_update_datasets the_clinical_index={the_clinical_index}", flush=True)
    print(f"convert_update_datasets the_datafile_index={the_datafile_index}", flush=True)
    print(f"convert_update_datasets the_download_dir={the_download_dir}", flush=True)
    print(f"convert_update_datasets the_biospecimen_dir={the_biospecimen_dir}", flush=True)
    print(f"convert_update_datasets the_clinical_dir={the_clinical_dir}", flush=True)
    print(f"convert_update_datasets the_convert_dir={the_convert_dir}", flush=True)
    print(f"convert_update_datasets the_sample_dir={the_sample_dir}", flush=True)
    print(f"convert_update_datasets the_util_dir={the_util_dir}", flush=True)
    print(f"convert_update_datasets the_temp_dir={the_temp_dir}", flush=True)
    # read datafile index
    my_datafiles: Dict[str, GdcApiDatafile] = read_datafile_index(the_datafile_index)
    # read biospecimen index
    # my_biospecimen: Dict[str, GdcApiBiospecimen] = read_biospecimen_index(the_biospecimen_index)
    # read clinical index
    # my_clinical: Dict[str, GdcApiClinical] = read_clinical_index(the_clinical_index)
    # get list of currently "released" datasets
    my_datasets: Dict[str, Dataset] = build_datasets(my_datafiles, the_biospecimen_dir)
    # reduce list to unconverted datasets or ones with new versions
    my_datasets = get_unconverted_datasets(my_datasets, the_convert_dir)
    # build new/updated datasets
    convert_datasets(my_datasets, the_util_dir, the_sample_dir, the_download_dir,
                     the_biospecimen_dir, the_clinical_dir, the_convert_dir, the_temp_dir)
    print("*************DONE:convert_update_datasets*************", flush=True)
# pylint: enable=too-many-arguments


# pylint: disable=too-many-locals
# noinspection DuplicatedCode
def convert_biospecimen_files(the_biospecimen_index: str,
                              the_file_datafiles: str, the_sample_dir: str,
                              the_download_dir: str, the_convert_dir: str) -> None:
    """
    Convert biospecimen files to batch files by Project and Program.
    Uses dataset/sample index information to get link between barcode and UUIDs, since some
    biospecimen information only has barcodes.
    :param the_biospecimen_index: full path and filename to biospecimen index
    :param the_file_datafiles: full path and filename to index file of dataset files
    :param the_sample_dir: full path to sample dir, which are associated with the dataset index (REQUIRED for some conversions)
    :param the_download_dir: full path to directory with downloaded biospecimen files
    :param the_convert_dir: full path to directory to place Program/Project batch files.
    :return: Nothing
    """
    print("*************START:convert_biospecimen_files*************", flush=True)
    print(f"convert_biospecimen_files the_biospecimen_index={the_biospecimen_index}", flush=True)
    print(f"convert_biospecimen_files the_file_datafiles={the_file_datafiles}", flush=True)
    print(f"convert_biospecimen_files the_download_dir={the_download_dir}", flush=True)
    print(f"convert_biospecimen_files the_convert_dir={the_convert_dir}", flush=True)
    # mapping uses (program, project) as key
    mapping_datafiles: Dict[Tuple[str, str], List[GdcApiDatafile]] = read_map_datafiles_by_program_projects(the_file_datafiles)
    # use tuple key (program, project) for organizing files to be converted
    files_to_convert: Dict[(str, str), List[GdcApiBiospecimen]] = {}
    # read index of biospecimen files to convert
    my_entries: Dict[str, GdcApiBiospecimen] = read_biospecimen_index(the_biospecimen_index)
    print(f"convert_biospecimen_files len(my_entries)={len(my_entries)}", flush=True)
    # divide entries up into program-project groups of only released data
    biospecimen: GdcApiBiospecimen
    index: int
    for index, biospecimen in enumerate(my_entries.values()):
        if 0 == index % 100:
            print(f"convert_biospecimen_files biospecimen index={index}", flush=True)
        if 'released' == biospecimen.history_status:
            converter: str = biospecimen.get_convertable_type()
            if converter.startswith("Unknown"):
                add_error(f"convert_biospecimen_files unknown file type {biospecimen.program} - {biospecimen.project} - {biospecimen.file_name}")
            elif converter.startswith("convert"):
                key: Tuple[str, str] = (biospecimen.program, biospecimen.project)
                biospecimen_list: List[GdcApiBiospecimen] = files_to_convert.get(key, [])
                biospecimen_list.append(biospecimen)
                files_to_convert[key] = biospecimen_list
    print(f"convert_biospecimen_files len(files_to_convert)={len(files_to_convert)}", flush=True)
    # iterate through each project-program list
    key_var: Tuple[str, str]
    val_list: List[GdcApiBiospecimen]
    for key_var, val_list in files_to_convert.items():
        print(f"convert_biospecimen_files key_var={key_var}", flush=True)
        # find largest release date (use as version)
        version: str = max([biospecimen.history_release_date for biospecimen in val_list])
        batches_file: str = os.path.join(the_convert_dir, key_var[0], key_var[1], version, 'batches.tsv')
        if not os.path.exists(batches_file):
            convert_biospecimen_batches(the_download_dir, batches_file, val_list,
                                        mapping_datafiles.get(key_var, []), the_sample_dir)
    print("*************DONE:convert_biospecimen_files*************", flush=True)
# pylint: enable=too-many-locals


# pylint: disable=too-many-locals
# noinspection DuplicatedCode
def convert_clinical_files(the_clinical_index: str, the_download_dir: str, the_convert_dir: str) -> None:
    """
    Convert public clinical files into annotation type batch files.
    :param the_clinical_index: full path and filename of clinical index file.
    :param the_download_dir: full path to download directory for clinical files.
    :param the_convert_dir: full path to convert directory for resulting annotation type batch files.
    :return: Nothing
    """
    print("*************START:convert_clinical_files*************", flush=True)
    print(f"convert_clinical_files the_clinical_index={the_clinical_index}", flush=True)
    print(f"convert_clinical_files the_download_dir={the_download_dir}", flush=True)
    print(f"convert_clinical_files the_convert_dir={the_convert_dir}", flush=True)
    # first key is program, second key is project
    files_to_convert: Dict[str, Dict[str, List[GdcApiClinical]]] = {}
    # divide entries up into program-project groups of only released data
    my_entries: Dict[str, GdcApiClinical] = read_clinical_index(the_clinical_index)
    print(f"convert_clinical_files len(my_entries)={len(my_entries)}", flush=True)
    # populate files_to_convert
    clinical: GdcApiClinical
    index: int
    for index, clinical in enumerate(my_entries.values()):
        if 0 == index % 100:
            print(f"convert_clinical_files clinical index={index}", flush=True)
        if 'released' == clinical.history_status:
            converter: str = clinical.get_convertable_type()
            if converter.startswith("Unknown"):
                add_error(f"convert_clinical_files unknown file type {clinical.program} - {clinical.project} - {clinical.file_name}")
            elif converter.startswith("convert"):
                program: str = clinical.program
                project: str = clinical.project
                project_dict: Dict[str, List[GdcApiClinical]] = {}
                project_list: List[GdcApiClinical] = []
                if program in files_to_convert:
                    project_dict = files_to_convert[program]
                    if project in project_dict:
                        project_list = project_dict[project]
                project_list.append(clinical)
                project_dict[project] = project_list
                files_to_convert[program] = project_dict
    print(f"convert_clinical_files len(files_to_convert)={len(files_to_convert)}", flush=True)
    # iterate through each project-program
    my_program: str
    project_dict: Dict[str, List[GdcApiClinical]]
    for my_program, project_dict in files_to_convert.items():
        print(f"convert_clinical_files my_program={my_program}", flush=True)
        my_project: str
        project_list: List[GdcApiClinical]
        for my_project, project_list in project_dict.items():
            print(f"convert_clinical_files my_project={my_project}", flush=True)
            # find largest release date (use as version)
            version: str = max([clinical.history_release_date for clinical in project_list])
            clinical_file: str = os.path.join(the_convert_dir, my_program, my_project, version, 'clinical.tsv')
            if not os.path.exists(clinical_file):
                convert_clinical_batches(the_download_dir, clinical_file, project_list)
    print("*************DONE:convert_clinical_files*************", flush=True)
# pylint: enable=too-many-locals
