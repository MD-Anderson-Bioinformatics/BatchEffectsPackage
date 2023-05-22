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

from mbatch.gdcapi.gdcapi import update_datafile_index, download_datafile_index
from mbatch.gdcapi.gdcapi import update_biospecimen_index, download_biospecimen_index
from mbatch.gdcapi.gdcapi import update_clinical_index, download_clinical_index
from mbatch.gdcapi.gdcapi import convert_update_datasets, convert_biospecimen_files, convert_clinical_files
from mbatch.test.common import print_errors, print_warnings

index_file_clinical: str = '/your/data/path/Pipeline-GDC/indexes/gdc_clinical.tsv'
index_file_biospecimen: str = '/your/data/path/Pipeline-GDC/indexes/gdc_biospecimen.tsv'
index_file_datafiles: str = '/your/data/path/Pipeline-GDC/indexes/gdc_datafiles.tsv'
index_file_samples_dir: str = '/your/data/path/Pipeline-GDC/indexes/samples'

dir_download_datafiles: str = '/your/data/path/Pipeline-GDC/downloaded/data'
dir_download_biospecimen: str = '/your/data/path/Pipeline-GDC/downloaded/biospecimen'
dir_download_clinical: str = '/your/data/path/Pipeline-GDC/downloaded/clinical'

dir_convert_datafiles: str = '/your/data/path/Pipeline-GDC/converted/data'
dir_convert_biospecimen: str = '/your/data/path/Pipeline-GDC/converted/biospecimen'
dir_convert_clinical: str = '/your/data/path/Pipeline-GDC/converted/clinical'

dir_convert_utils: str = '/your/data/path/Pipeline-GDC/util'
dir_convert_temp: str = '/your/data/path/Pipeline-GDC/CONVERT_TEMP'

if __name__ == '__main__':
    # ########################################################
    # download information about biospecimen files
    # ########################################################
    update_biospecimen_index(index_file_biospecimen)
    # ########################################################
    # download information about (public) clinical files
    # ########################################################
    update_clinical_index(index_file_clinical)
    # ########################################################
    # download information about datafiles and their samples
    # ########################################################
    update_datafile_index(index_file_datafiles, index_file_samples_dir)
    # ########################################################
    # download new clinical files
    # ########################################################
    download_clinical_index(index_file_clinical, dir_download_clinical)
    # ########################################################
    # download new biospecimen files
    # ########################################################
    download_biospecimen_index(index_file_biospecimen, dir_download_biospecimen)
    # ########################################################
    # download new datafiles
    # ########################################################
    download_datafile_index(index_file_datafiles, dir_download_datafiles, index_file_samples_dir)
    # ########################################################
    # convert biospecimen files
    # ########################################################
    convert_biospecimen_files(index_file_biospecimen, index_file_datafiles,
                              index_file_samples_dir,
                              dir_download_biospecimen, dir_convert_biospecimen)
    # ########################################################
    # convert clinical files
    # ########################################################
    convert_clinical_files(index_file_clinical, dir_download_clinical, dir_convert_clinical)
    # ########################################################
    # convert datasets
    # ########################################################
    convert_update_datasets(index_file_datafiles, index_file_biospecimen, index_file_clinical,
                            dir_download_datafiles, dir_convert_biospecimen, dir_convert_clinical, dir_convert_datafiles,
                            index_file_samples_dir, dir_convert_utils, dir_convert_temp)
    # ########################################################
    # print important warnings and errors
    # ########################################################
    print("-----------------------------------------", flush=True)
    print("-------------WARNING NOTES-------------", flush=True)
    print_warnings()
    print("-------------ERROR NOTES-------------", flush=True)
    print_errors()
    print("-------------gdcapi_main:done-------------", flush=True)
