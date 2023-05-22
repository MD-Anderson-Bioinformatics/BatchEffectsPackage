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
import os
import json
from typing import List, Dict
import requests
from mbatch.test.common import write_tsv_list, read_headers
from mbatch.gdcapi.download_history_mixin import GdcApiHistoryMixin, write_headers_history
from mbatch.gdcapi.gdcapi_globals import GLOBAL_MAX_COUNT, GLOBAL_DLD_COUNT


# ########################################################
# Files Endpoint, getting downloadable data files
# ########################################################

# GDC API endpoint for files API
BIOSPECIMEN_ENDPOINT: str = 'https://api.gdc.cancer.gov/files'

# GDC API filter for release open-access data. JSON format.
BIOSPECIMEN_FILTER: json = {
    'op': 'and',
    'content':
        [
            {
                'op': '=',
                'content':
                    {
                        'field': 'state',
                        'value': 'released'
                    }
            },
            {
                'op': '!=',
                'content':
                    {
                        'field': 'type',
                        'value': 'clinical_supplement'
                    }
            },
            {
                'op': '=',
                'content':
                    {
                        'field': 'data_category',
                        'value': 'Biospecimen'
                    }
            },
            {
                'op': '=',
                'content':
                    {
                        'field': 'access',
                        'value': 'open'
                    }
            },
            {
                'op': 'in',
                'content':
                    {
                        'field': 'data_format',
                        'value': ['XLSX', 'TSV', 'BCR XML']
                    }
            }
        ]
    }

# GDC API list of fields to include in results
BIOSPECIMEN_FIELDS: str = '' + \
    'file_id,' + \
    'md5sum,' + \
    'file_name,' + \
    'data_type,' + \
    'data_format,' + \
    'data_category,' + \
    'cases.project.project_id,' + \
    'cases.project.program.name'

# GDC API parameters for HTTPS request.
# added expand so entries get populated even without end data
BIOSPECIMEN_PARAMS: Dict[str, str] = {
    'pretty': 'true',
    'size': '50',
    'filters': json.dumps(BIOSPECIMEN_FILTER),
    'fields': BIOSPECIMEN_FIELDS,
    'expand': '',
    'from': '0'
}

BIOSPECIMEN_HEADERS: List[str] = [
    'data_format',
    'project',
    'program'
]


# pylint: disable=too-many-instance-attributes,too-many-arguments,too-few-public-methods
class GdcApiBiospecimen(GdcApiHistoryMixin):
    """
    Biospecimen files from GDC, have data_format, project, and program data at this level.
    Uses history details inherited from GdcApiHistoryMixin.
    data_format - BCR XML
    program - top level descriptor, like TCGA, TARGET, etc.
    project - similar to disease type. Usually, program plus disease or cohort
    """
    # declare but do not set member attributes
    data_format: str
    project: str
    program: str

    def __init__(self: 'GdcApiBiospecimen', the_dict: Dict, the_json_flag: bool) -> None:
        """
        Initialize GdcApiBiospecimen from JSON or TSV built Dictionary.
        JSON data will be for non-History labeled attributes.
        :param the_dict: dictionary with attribute values
        :param the_json_flag: if True Dictionary is from JSON, otherwise from TSV
        """
        super().__init__(the_dict, the_json_flag)
        self.data_format: str = ""
        self.project: str = ""
        self.program: str = ""
        if the_json_flag:
            # process JSON Dictionary
            self.data_format = the_dict['data_format']
            self.project = ''
            self.program = ''
            case_entry: List[Dict] = the_dict['cases']
            case_dict: Dict
            for case_dict in case_entry:
                if '' == self.project:
                    self.project = case_dict['project']['project_id']
                    self.program = case_dict['project']['program']['name']
        else:
            # process TSV Dict
            entry: str
            for entry in BIOSPECIMEN_HEADERS:
                setattr(self, entry, the_dict.get(entry))

    def write_biospecimen(self: 'GdcApiBiospecimen', the_out_file: io.TextIOWrapper) -> None:
        """
        Write out just the biospecimen portion of a file element to a file,
        after first calling parent history to write out its portion.
        Writes a newline.
        :param the_out_file: stream to which to write strings
        :return: None
        """
        data_list: List[str] = []
        for attr in BIOSPECIMEN_HEADERS:
            data_list.append(getattr(self, attr))
        self.write_history(the_out_file)
        write_tsv_list(the_out_file, data_list, True, False)

    def get_download_file(self: 'GdcApiBiospecimen', the_base_dir: str) -> str:
        """
        Build the "biospecimen" version of a file path.
        Builds /base-dir/program/project/version/file.name path.
        :param the_base_dir: full path to base directory to add specialized path and filename
        :return: full path with filename to which to download the file
        """
        # print("GdcApiBiospecimen::get_download_file", flush=True)
        return os.path.join(the_base_dir,
                            self.program, self.project, self.history_version,
                            self.file_name)

    def get_file_archive(self: 'GdcApiBiospecimen', the_base_dir: str) -> str:
        """
        Build the "datafile" version of a file path for the ZIP archive.
        Builds /base-dir/program/project/experimental_strategy/data_type_display/workflow/version/uuid.zip
        :param the_base_dir: full path to base directory to add specialized path and filename
        :return: full path with filename which is the ZIP archive
        """
        # print("GdcApiDatafile::get_download_file", flush=True)
        return os.path.join(the_base_dir,
                            self.program, self.project, self.history_version,
                            self.uuid + ".zip")

    def get_convertable_type(self: 'GdcApiBiospecimen') -> str:
        """
        Determine conversion type (as a string) for this GdcApiBiospecimen.
        Conversion based on data_format and file_names.
        (Some file_name checking looks for program names.)
        :return: String giving conversion method. 'Unknown' if no matches.
        """
        ret_val: str = 'Unknown'
        # order matters for this!
        if ('BCR XML' == self.data_format) & (self.file_name.endswith('.xml')):
            ret_val = 'convert-bcr_xml'
        elif ('XLSX' == self.data_format) & (self.file_name.endswith('.xlsx')) & \
                (self.file_name.startswith('CGCI')) & ('biospecimenSupplementSpreadsheet' in self.file_name):
            ret_val = 'convert-cgci_xlsx'
        elif ('XLSX' == self.data_format) & (self.file_name.endswith('.xlsx')) & \
                (self.file_name.startswith('CGCI')) & ('biospecimenSupplementData' in self.file_name):
            ret_val = 'convert-cgci_xlsx'
        elif ('TSV' == self.data_format) & (self.file_name.endswith('.tsv')) & \
                (self.file_name.startswith('FM')) & ('_Biospecimen' in self.file_name):
            ret_val = 'convert-fm_tsv'
        elif ('XLSX' == self.data_format) & (self.file_name.endswith('.xlsx')) & \
                (self.file_name.startswith('TARGET')) & ('SampleMatrix' in self.file_name):
            ret_val = 'convert-target_xlsx'
        return ret_val

# ########################################################
# retrieve data from GDC HTTP API
# ########################################################


def get_biospecimen_list_from_gdc(the_entries: Dict[str, GdcApiBiospecimen]) -> List[GdcApiBiospecimen]:
    """
    Call the file endpoint, looping until we get no results on the new page, populating
    list of GdcApiBiospecimen. Use existing object from dictionary if present.
    :param the_entries: Dictionary of GdcApiBiospecimen objects.
    :return: Updated list of GdcApiBiospecimen objects.
    """
    results: List[GdcApiBiospecimen] = []
    # development max size limit (0 does not limit download)
    max_size: int = GLOBAL_MAX_COUNT
    size: int = GLOBAL_DLD_COUNT
    from_int: int = 0
    continue_while: bool = True
    print("get_biospecimen_list_from_gdc start", flush=True)
    while continue_while:
        print(f"loop size={size} from_int={from_int}", flush=True)
        BIOSPECIMEN_PARAMS['size'] = f'{size}'
        BIOSPECIMEN_PARAMS['from'] = f'{from_int}'
        response: requests = requests.get(BIOSPECIMEN_ENDPOINT, params=BIOSPECIMEN_PARAMS, timeout=60, allow_redirects=True)
        print(f"get_biospecimen_list_from_gdc response={response}", flush=True)
        my_str: str = json.dumps(response.json(), indent=2)
        my_dict: Dict = json.loads(my_str)
        datafile_list: List[Dict] = my_dict['data']['hits']
        my_entry: Dict
        for my_entry in datafile_list:
            new_item: GdcApiBiospecimen = GdcApiBiospecimen(my_entry, True)
            if new_item.uuid in the_entries:
                # if existing entry found, use existing, saves some calls to history
                new_item = the_entries[new_item.uuid]
            results.append(new_item)
        len_datafile_list: int = len(datafile_list)
        print(f"get_biospecimen_list_from_gdc len(len_datafile_list)={len_datafile_list}", flush=True)
        if len_datafile_list < 1:
            continue_while = False
        # development short version
        if max_size > 0:
            if len(results) >= max_size:
                continue_while = False
        from_int += size
    # end of loop
    print("get_biospecimen_list_from_gdc done", flush=True)
    return results

# ########################################################
# writing files
# ########################################################


def write_headers_biospecimen(the_out_file: io.TextIOWrapper) -> None:
    """
    Write the headers for biospecimen portion of object,
    after calling to write History portion.
    Writes newline at the end.
    :param the_out_file: stream to which to write headers
    :return: None
    """
    write_headers_history(the_out_file)
    write_tsv_list(the_out_file, BIOSPECIMEN_HEADERS, True, False)


def write_biospecimen_index(the_index_file: str, the_entries: List[GdcApiBiospecimen]) -> None:
    """
    Write the headers and line items for index files.
    Replaces existing index file.
    :param the_index_file: full path to index file to write.
    :param the_entries: list of GdcApiBiospecimen to write to index.
    :return: None
    """
    print(f"update_biospecimens the_index_file={the_index_file}", flush=True)
    the_entries.sort(key=lambda my_biospecimen: (my_biospecimen.program, my_biospecimen.project, my_biospecimen.history_version), reverse=False)
    out_file: io.TextIOWrapper
    with open(the_index_file, 'w', encoding='utf-8') as out_file:
        write_headers_biospecimen(out_file)
        my_entry: GdcApiBiospecimen
        for my_entry in the_entries:
            my_entry.write_biospecimen(out_file)

# ########################################################
# reading files
# ########################################################


def read_biospecimen_index(the_index_file: str) -> Dict[str, GdcApiBiospecimen]:
    """
    Read index file and create Dictionary of GdcApiBiospecimen object.
    :param the_index_file: full path and filename of index file to read.
    :return: Dictionary with uuid as keys and GdcApiBiospecimen objects as values.
    """
    print(f"read_biospecimen_index the_index_file={the_index_file}", flush=True)
    in_file: io.TextIOWrapper
    my_entries: Dict[str, GdcApiBiospecimen] = {}
    with open(the_index_file, 'r', encoding='utf-8') as in_file:
        keys: List[str] = read_headers(in_file)
        line: str = in_file.readline().rstrip('\n')
        index: int = 1
        while '' != line:
            if 0 == index % 1000:
                print(f"read_biospecimen_index biospecimen index={index}", flush=True)
            index += 1
            values: List[str] = line.split("\t")
            tsv_dict: Dict[str, str] = dict(zip(keys, values))
            entry: GdcApiBiospecimen = GdcApiBiospecimen(tsv_dict, False)
            my_entries[entry.uuid] = entry
            line = in_file.readline().rstrip('\n')
    return my_entries

# ########################################################
#
# ########################################################
