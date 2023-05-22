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


import os
import io
import json
import hashlib
import typing
from typing import Dict, List, Set
import zipfile
from mbatch.gdcapi.download_datafile import GdcApiDatafile
from mbatch.gdcapi.download_biospecimen import GdcApiBiospecimen
from mbatch.gdcapi.download_clinical import GdcApiClinical
from mbatch.test.common import get_newest_dir
from mbatch.test.common import add_error, add_warnings


# pylint: disable=too-many-arguments
def calc_md5_dataset(the_source: str, the_project: str, the_sub_project: str, the_category: str,
                     the_platform: str, the_data: str, the_details: str) -> str:
    """
    Calcualte the MD5 sum for a dataset
    :param the_source: source for dataset (currently GDC or MW)
    :param the_project: project (program) name for dataset
    :param the_sub_project: project name for dataset
    :param the_category: dataset category
    :param the_platform: platform category
    :param the_data: data category
    :param the_details: details category
    :return: MD5 sum from above strings
    """
    my_hash: hashlib.md5 = hashlib.md5()
    my_hash.update(the_source.encode('utf-8'))
    my_hash.update(the_project.encode('utf-8'))
    my_hash.update(the_sub_project.encode('utf-8'))
    my_hash.update(the_category.encode('utf-8'))
    my_hash.update(the_platform.encode('utf-8'))
    my_hash.update(the_data.encode('utf-8'))
    my_hash.update(the_details.encode('utf-8'))
    return my_hash.hexdigest()
# pylint: enable=too-many-arguments


# pylint: disable=too-many-instance-attributes,too-many-arguments,too-few-public-methods
class Dataset:
    """
    Class to represent information about a single dataset.
    A dataset has multiple files going in, but generates a single matrix TSV

    Datasets are defined in the following way:
    datafile attribute ----> menu entry
    ==============================
    program ----> project
    project ----> sub-project
    experimental_strategy ----> category
    data_type_display ----> platform
    workflow ----> data
    history_version ----> Use as new "data-version" entry?
    """
    # declare but do not set member attributes
    unique: str
    source: str
    project: str
    sub_project: str
    category: str
    platform: str
    data: str
    details: str
    files: Dict[str, str]
    version: str

    def __init__(self: 'Dataset', the_source: str, the_project: str, the_sub_project: str, the_category: str,
                 the_platform: str, the_data: str, the_details: str) -> None:
        """
        Constructor for Dataset object.
        Files are set to empty.
        MD5 sum is calculated.
        Version is set later.
        :param the_source: source for dataset (currently GDC or MW)
        :param the_project: program name for dataset
        :param the_sub_project: project name for dataset
        :param the_category: dataset category
        :param the_platform: platform category
        :param the_data: data category
        :param the_details: datails category
        """
        super().__init__()
        self.unique: str = calc_md5_dataset(the_source, the_project, the_sub_project, the_category, the_platform, the_data, the_details)
        self.source: str = the_source
        self.project: str = the_project
        self.sub_project: str = the_sub_project
        self.category: str = the_category
        self.platform: str = the_platform
        self.data: str = the_data
        self.details: str = the_details
        self.files: Dict[str, str] = {}
        self.version: str = ""

    def add_with_check(self: 'Dataset', the_datafile: GdcApiDatafile) -> bool:
        """
        Add GdcApiDatafile if it belongs to this Dataset.
        For GdcApiDatafile object, get MD5 sum and compare it against the Dataset MD5.
        If it matched, check that GdcApiDatafile has a "released" status.
        Then add to dictionary of files.
        :param the_datafile: GdcApiDatafile to potentially add to file Dictionary
        :return: True is GdcApiDatafile was added. Otherwise False.
        """
        added: bool = False
        md5dataset: str = calc_md5_dataset(the_datafile.source,
                                           the_datafile.program, the_datafile.project,
                                           the_datafile.experimental_strategy,
                                           the_datafile.data_type_display,
                                           the_datafile.workflow, self.details)
        if self.unique == md5dataset:
            if 'released' == the_datafile.history_status:
                added = True
                self.files[the_datafile.uuid] = the_datafile
        return added

    def determine_version(self: 'Dataset') -> None:
        """
        Determine version of Dataset by finding newest release data
        from component GdcApiDatafile objects.
        Set version to that.
        :return: Nothing
        """
        max_release: str = ""
        my_file: GdcApiDatafile
        for my_file in self.files.values():
            if my_file.history_release_date > max_release:
                max_release = my_file.history_release_date
        self.version = max_release

    def get_dataset_path(self: 'Dataset', the_convert_dir: dir) -> str:
        """
        Build output directory for writing converted Dataset file.
        :param the_convert_dir: Full path to base convert directory
        :return: New path
        """
        return os.path.join(the_convert_dir, self.project, self.sub_project, self.category,
                            self.platform, self.data, self.details, self.version)

    def get_biospecimen(self: 'Dataset', the_biospecimens: Dict[str, GdcApiBiospecimen]) -> GdcApiBiospecimen:
        """
        Find newest GdcApiBiospecimen object (using history_release_date).
        :param the_biospecimens: Dictionary of GdcApiBiospecimen objects
        :return: newest GdcApiBiospecimen object (using history_release_date)
        """
        my_biospecimen: typing.Optional[GdcApiBiospecimen] = None
        tmp: GdcApiBiospecimen
        for tmp in the_biospecimens.values():
            if self.project == tmp.program:
                if self.sub_project == tmp.project:
                    if my_biospecimen is None:
                        my_biospecimen = tmp
                    elif tmp.history_release_date > my_biospecimen.history_release_date:
                        my_biospecimen = tmp
        return my_biospecimen

    def get_clinical(self: 'Dataset', the_clinicals: Dict[str, GdcApiClinical]) -> GdcApiClinical:
        """
        Find newest GdcApiClinical object (using history_release_date).
        :param the_clinicals: Dictionary of GdcApiClinical objects
        :return: newest GdcApiClinical object (using history_release_date)
        """
        my_clinical: typing.Optional[GdcApiClinical] = None
        tmp: GdcApiClinical
        for tmp in the_clinicals.values():
            if self.project == tmp.program:
                if self.sub_project == tmp.project:
                    if my_clinical is None:
                        my_clinical = tmp
                    elif tmp.history_release_date > my_clinical.history_release_date:
                        my_clinical = tmp
        return my_clinical

    def conversion_match(self: 'Dataset', the_experimental_strategy: str, the_data_type_display: str,
                         the_workflow: str, the_details: str) -> bool:
        """
        Check if Dataset's category, platform, and data match the ones passed in
        :param the_experimental_strategy: compared against self.category
        :param the_data_type_display: compared against self.platform
        :param the_workflow: compared against self.data
        :param the_details: compared against self.details
        :return: True is the three values match. False otherwise.
        """
        matched: bool = False
        if the_experimental_strategy == self.category:
            if the_data_type_display == self.platform:
                if the_workflow == self.data:
                    if the_details == self.details:
                        matched = True
        return matched

    def get_batches_tsv_path(self: 'Dataset', the_convert_dir: dir) -> typing.Optional[str]:
        """
        Build full path to batches.tsv for the dataset
        :param the_convert_dir: Full path to base/biospecimen convert directory
        :return: New path and batches.tsv filename
        """
        version_path: str = get_newest_dir(os.path.join(the_convert_dir, self.project, self.sub_project))
        if version_path is not None:
            version_path = os.path.join(version_path, "batches.tsv")
        return version_path

    def get_clinical_tsv_path(self: 'Dataset', the_convert_dir: dir) -> typing.Optional[str]:
        """
        Build full path to clinical.tsv for the dataset
        :param the_convert_dir: Full path to base/clinical convert directory
        :return: New path and clinical.tsv filename
        """
        version_path: str = get_newest_dir(os.path.join(the_convert_dir, self.project, self.sub_project))
        if version_path is not None:
            version_path = os.path.join(version_path, "clinical.tsv")
        return version_path

    def get_archive_path(self: 'Dataset', the_convert_dir: dir) -> str:
        """
        Build output directory for storing ZIP archive.
        :param the_convert_dir: Full path to base convert directory
        :return: New path
        """
        return os.path.join(the_convert_dir, self.project, self.sub_project, self.category,
                            self.platform, self.data, self.details)

    def get_archive_filename(self: 'Dataset') -> str:
        """
        Build ZIP archive filename.
        :return: File name of sub_project plus "_" plus details
        """
        return self.sub_project + "_" + self.details + ".zip"

    def get_archive_internal_data_path(self: 'Dataset') -> str:
        """
        Build ZIP archive filename. Include slash on end, which ZIP file uses.
        :return: versions/DATA_xxx/
        """
        return "versions/DATA_" + self.version + "/"

    def dump_to_dict_file(self: 'Dataset', the_file: str) -> None:
        """
        Write to a JSON file
        :param the_file: full path and filename to which to dump class as dictionary
        :return: Nothing
        """
        my_dict: Dict[str, str] = \
            {
                'program': self.project,
                'project': self.sub_project,
                'category': self.category,
                'platform': self.platform,
                'data': self.data,
                'details': self.details
            }
        out_file: io.TextIOWrapper
        with open(the_file, 'w', encoding='utf-8') as out_file:
            json.dump(my_dict, out_file)
# pylint: enable=too-many-instance-attributes,too-many-arguments,too-few-public-methods


# pylint: disable=too-many-branches,too-many-statements
def build_datasets(the_datafiles: Dict[str, GdcApiDatafile], the_biospecimen_dir: str) -> Dict[str, Dataset]:
    """
    Build a Dataset dictionary from a dictionary of GdcApiDatafile objects.
    GdcApiDatafile specifies a single GDC file and associated metadata.
    Dataset is a single matrix.tsv in StandardizedData,
    with usually multiple GdcApiDatafiles associated with it.
    Dataset dictionary breaks the GdcApiDatafile objects into datasets
    based on program, project, workflow, platform, etc.
    :param the_datafiles: Dictionary of GdcApiDatafile objects
    :param the_biospecimen_dir: full path to directory with biospecimen index files
    :return: Binned dictionary of Dataset objects, each object representing a matrix.tsv file.
    """
    unknown_conversions: Set[str] = set()
    dataset_dict: Dict[str, Dataset] = {}
    my_datafile: GdcApiDatafile
    for my_datafile in the_datafiles.values():
        if 'released' == my_datafile.history_status:
            dataset_list: List[Dataset] = []
            compare: List[str] = [my_datafile.experimental_strategy, my_datafile.data_type_display, my_datafile.workflow]
            workflow: str = my_datafile.workflow
            if ["Genotyping Array", "Copy Number Segment", "DNAcopy"] == compare:
                dataset_list.append(Dataset(my_datafile.source, my_datafile.program, my_datafile.project,
                                            my_datafile.experimental_strategy,
                                            my_datafile.data_type_display,
                                            my_datafile.workflow, "Copy-Number-with-CNV"))
            elif ["Genotyping Array", "Masked Copy Number Segment", "DNAcopy"] == compare:
                dataset_list.append(Dataset(my_datafile.source, my_datafile.program, my_datafile.project,
                                            my_datafile.experimental_strategy,
                                            my_datafile.data_type_display,
                                            my_datafile.workflow, "Copy-Number-no-CNV"))
            elif ["RNA-Seq", "Gene Expression Quantification", "STAR - Counts"] == compare:
                dataset_list.append(Dataset(my_datafile.source, my_datafile.program, my_datafile.project,
                                            my_datafile.experimental_strategy,
                                            my_datafile.data_type_display,
                                            my_datafile.workflow, "RNASeq-TPM"))
                dataset_list.append(Dataset(my_datafile.source, my_datafile.program, my_datafile.project,
                                            my_datafile.experimental_strategy,
                                            my_datafile.data_type_display,
                                            my_datafile.workflow, "RNASeq-FPKM"))
                dataset_list.append(Dataset(my_datafile.source, my_datafile.program, my_datafile.project,
                                            my_datafile.experimental_strategy,
                                            my_datafile.data_type_display,
                                            my_datafile.workflow, "RNASeq-FPKM-UQ"))
            elif ["Methylation Array", "Methylation Beta Value", "SeSAMe Methylation Beta Estimation"] == compare:
                dataset_list.append(Dataset(my_datafile.source, my_datafile.program, my_datafile.project,
                                            my_datafile.experimental_strategy,
                                            my_datafile.data_type_display,
                                            my_datafile.workflow, "Methylation-With-Sex-Chromosomes"))
                dataset_list.append(Dataset(my_datafile.source, my_datafile.program, my_datafile.project,
                                            my_datafile.experimental_strategy,
                                            my_datafile.data_type_display,
                                            my_datafile.workflow, "Methylation-No-Sex-Chromosomes"))
            elif ["Reverse Phase Protein Array", "Protein Expression Quantification", "Protein Analysis"] == compare:
                dataset_list.append(Dataset(my_datafile.source, my_datafile.program, my_datafile.project,
                                            my_datafile.experimental_strategy,
                                            my_datafile.data_type_display,
                                            my_datafile.workflow, "RPPA"))
            elif ["miRNA-Seq", "Isoform Expression Quantification", "BCGSC miRNA Profiling"] == compare:
                dataset_list.append(Dataset(my_datafile.source, my_datafile.program, my_datafile.project,
                                            my_datafile.experimental_strategy,
                                            my_datafile.data_type_display,
                                            my_datafile.workflow, "miRNA-Isoform"))
            elif ["miRNA-Seq", "miRNA Expression Quantification", "BCGSC miRNA Profiling"] == compare:
                dataset_list.append(Dataset(my_datafile.source, my_datafile.program, my_datafile.project,
                                            my_datafile.experimental_strategy,
                                            my_datafile.data_type_display,
                                            my_datafile.workflow, "miRNA-Genes"))
            elif ["Genotyping Array", "Gene Level Copy Number", "ASCAT2"] == compare:
                if not Dataset(my_datafile.source, my_datafile.program, my_datafile.project, my_datafile.experimental_strategy, my_datafile.data_type_display, my_datafile.workflow, "null").get_batches_tsv_path(the_biospecimen_dir) is None:
                    dataset_list.append(Dataset(my_datafile.source, my_datafile.program, my_datafile.project,
                                                my_datafile.experimental_strategy,
                                                my_datafile.data_type_display,
                                                my_datafile.workflow, "Copy-Number-With-Sex-Chromosomes"))
                    dataset_list.append(Dataset(my_datafile.source, my_datafile.program, my_datafile.project,
                                                my_datafile.experimental_strategy,
                                                my_datafile.data_type_display,
                                                my_datafile.workflow, "Copy-Number-No-Sex-Chromosomes"))
                else:
                    add_warnings(f"Do not process Genotyping Array - Gene Level Copy Number - ASCAT2 - no biospecimen index for {my_datafile.program} - {my_datafile.project} - {my_datafile.history_release_date}")
            elif ["Genotyping Array", "Allele-specific Copy Number Segment", "ASCAT2"] == compare:
                add_warnings("Do not process Allele-specific Copy Number Segment dataset - not useful")
            elif ["WGS", "Gene Level Copy Number", "AscatNGS"] == compare:
                # check if biospecimen file exists
                if not Dataset(my_datafile.source, my_datafile.program, my_datafile.project, my_datafile.experimental_strategy, my_datafile.data_type_display, my_datafile.workflow, "null").get_batches_tsv_path(the_biospecimen_dir) is None:
                    dataset_list.append(Dataset(my_datafile.source, my_datafile.program, my_datafile.project,
                                                my_datafile.experimental_strategy,
                                                my_datafile.data_type_display,
                                                my_datafile.workflow, "Copy-Number-With-Sex-Chromosomes"))
                    dataset_list.append(Dataset(my_datafile.source, my_datafile.program, my_datafile.project,
                                                my_datafile.experimental_strategy,
                                                my_datafile.data_type_display,
                                                my_datafile.workflow, "Copy-Number-No-Sex-Chromosomes"))
                else:
                    add_warnings(f"Do not process WGS - Gene Level Copy Number - AscatNGS dataset - no biospecimen index for {my_datafile.program} - {my_datafile.project} - {my_datafile.history_release_date}")
            elif ["WGS", "Copy Number Segment", "AscatNGS"] == compare:
                add_warnings("Do not process WGS Copy Number Segment dataset - not useful")
            elif "Aliquot Ensemble Somatic Variant Merging and Masking" == workflow:
                dataset_list.append(Dataset(my_datafile.source, my_datafile.program, my_datafile.project,
                                            my_datafile.experimental_strategy,
                                            my_datafile.data_type_display,
                                            my_datafile.workflow, "Mutation-Calling"))
            elif ["scRNA-Seq", "Single Cell Analysis", "Seurat - 10x Chromium"] == compare:
                add_warnings("Do not process Single Cell Analysis dataset - this is PCA results")
            elif ["scRNA-Seq", "Differential Gene Expression", "Seurat - 10x Chromium"] == compare:
                dataset_list.append(Dataset(my_datafile.source, my_datafile.program, my_datafile.project,
                                            my_datafile.experimental_strategy,
                                            my_datafile.data_type_display,
                                            my_datafile.workflow, "scRNA-Differential"))
            else:
                # added to add_error later
                unknown_conversions.add(my_datafile.experimental_strategy + "|" + my_datafile.data_type_display + "|" + my_datafile.workflow)
            my_dataset: Dataset
            for my_dataset in dataset_list:
                if my_dataset.unique in dataset_dict:
                    # dataset exists, get it
                    my_dataset = dataset_dict[my_dataset.unique]
                # add file to dataset
                my_dataset.add_with_check(my_datafile)
                dataset_dict[my_dataset.unique] = my_dataset
    my_dataset: Dataset
    for my_dataset in dataset_dict.values():
        my_dataset.determine_version()
    my_str: str
    for my_str in unknown_conversions:
        add_error(f"Unknown Conversion {my_str}")
    return dataset_dict
# pylint: enable=too-many-branches,too-many-statements


def get_unconverted_datasets(the_datasets: Dict[str, Dataset], the_convert_dir: str) -> Dict[str, Dataset]:
    """
    Check if the dataset converted path for a Dataset exists, return a dictionary
    containing only Dataset objects that have not been converted.
    :param the_datasets: Dictionary of Dataset objects that may or many not have been converted
    :param the_convert_dir: Full path to base directory for converted files
    :return: Dictionary of Dataset objects are not converted
    """
    # redo get_unconverted_datasets to support updating existing and creating new
    # new version should work with ZIP archives
    unconverted: Dict[str, Dataset] = {}
    my_dataset: Dataset
    for my_dataset in the_datasets.values():
        archive_path: str = my_dataset.get_archive_path(the_convert_dir)
        archive_name: str = my_dataset.get_archive_filename()
        archive_filepath: str = os.path.join(archive_path, archive_name)
        # check if archive exists
        if not os.path.exists(archive_filepath):
            # if archive does not exist at all, add it to list to be converted
            unconverted[my_dataset.unique] = my_dataset
        else:
            # if archive exists, check if this version exists inside it
            #  inside archive is versions/DATA_<GDC history date version>
            found: bool = False
            int_path: str = my_dataset.get_archive_internal_data_path()
            zip_file: zipfile.ZipFile
            with zipfile.ZipFile(archive_filepath, 'r') as zip_file:
                int_file: str
                for int_file in zip_file.namelist():
                    # print(f"{int_file}", flush=True)
                    if int_file == int_path:
                        found = True
            if not found:
                # if archive DATA VERSION does not exist at all, add it to list to be converted
                unconverted[my_dataset.unique] = my_dataset
    return unconverted
