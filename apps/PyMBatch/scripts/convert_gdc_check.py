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
from mbatch.gdcapi.download_biospecimen import GdcApiBiospecimen
from mbatch.gdcapi.download_clinical import GdcApiClinical
from mbatch.gdcapi.download_datafile import GdcApiDatafile
from mbatch.gdcapi.gdcapi import update_biospecimen_index
from mbatch.gdcapi.gdcapi import update_clinical_index
from mbatch.gdcapi.gdcapi import update_datafile_index
from mbatch.test.common import print_errors, print_warnings

# update this for internal paths when running for release pipeline
index_file_clinical: str = '/BEA/DAPI_GDC/indexes/gdc_clinical.tsv'
index_file_biospecimen: str = '/BEA/DAPI_GDC/indexes/gdc_biospecimen.tsv'
index_file_datafiles: str = '/BEA/DAPI_GDC/indexes/gdc_datafiles.tsv'
index_file_samples_dir: str = '/BEA/DAPI_GDC/indexes/samples'


if __name__ == '__main__':
    # ########################################################
    # get updates for biospecimen files
    # ########################################################
    new_biospecimen: List[GdcApiBiospecimen] = update_biospecimen_index(index_file_biospecimen, False)
    # ########################################################
    # download information about (public) clinical files
    # ########################################################
    new_clinical: List[GdcApiClinical] = update_clinical_index(index_file_clinical, False)
    # ########################################################
    # download information about datafiles and their samples
    # ########################################################
    new_files: List[GdcApiDatafile] = update_datafile_index(index_file_datafiles, index_file_samples_dir, False)
    # ########################################################
    # print changes
    # ########################################################
    print("#### New Biospecimen Entries", flush=True)
    bio_entry: GdcApiBiospecimen
    for bio_entry in new_biospecimen:
        print(f"{bio_entry.project} {bio_entry.program} {bio_entry.file_name} {bio_entry.uuid} ", sep='\n', flush=True)
    print("#### New Clinical Entries", flush=True)
    cli_entry: GdcApiClinical
    for cli_entry in new_clinical:
        print(f"{cli_entry.project} {cli_entry.program} {cli_entry.file_name} {cli_entry.uuid} ", sep='\n', flush=True)
    print("#### New DataFile Entries", flush=True)
    data_entry: GdcApiDatafile
    for data_entry in new_files:
        print(f"{data_entry.project} {data_entry.program} {data_entry.file_name} {data_entry.uuid} ", sep='\n', flush=True)
    print(f"#### New Biospecimen Entries Count={len(new_biospecimen)}", flush=True)
    print(f"#### New Clinical Entries Count={len(new_clinical)}", flush=True)
    print(f"#### New DataFile Entries Count={len(new_files)}", flush=True)
    # ########################################################
    # print important warnings and errors
    # ########################################################
    print("-----------------------------------------", flush=True)
    print("-------------WARNING NOTES-------------", flush=True)
    print_warnings()
    print("-------------ERROR NOTES-------------", flush=True)
    print_errors()
    print("-------------check4_main:done-------------", flush=True)
