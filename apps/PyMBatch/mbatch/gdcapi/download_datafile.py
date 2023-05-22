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
import typing
from typing import Dict, List, Optional, Tuple
import pandas
import requests
from mbatch.gdcapi.download_history_mixin import GdcApiHistoryMixin, write_headers_history
from mbatch.gdcapi.download_sample import GdcApiSample, write_sample_index_file, read_sample_index_file
from mbatch.gdcapi.gdcapi_globals import GLOBAL_MAX_COUNT, GLOBAL_DLD_COUNT
from mbatch.test.common import write_tsv_list, read_headers, add_error


# ########################################################
# Files Endpoint, getting downloadable data files
# ########################################################

# GDC API endpoint for files API
DATAFILE_ENDPOINT: str = 'https://api.gdc.cancer.gov/files'

# GDC API filter for release open-access data. JSON format.
DATAFILE_FILTER: json = {
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
                'op': '!=',
                'content':
                    {
                        'field': 'type',
                        'value': 'biospecimen_supplement'
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
                'op': 'or',
                'content':
                [
                    {
                        'op': '=',
                        'content':
                        {
                            'field': 'data_format',
                            'value': 'TSV'
                        }
                    },
                    {
                        'op': '=',
                        'content':
                        {
                            'field': 'data_format',
                            'value': 'MAF'
                        }
                    },
                    {
                        'op': '=',
                        'content':
                        {
                            'field': 'data_format',
                            'value': 'TXT'
                        }
                    }
                ]
            }
        ]
    }

# GDC API list of fields to include in results
DATAFILE_FIELDS: str = '' + \
    'file_id,' + \
    'md5sum,' + \
    'file_name,' + \
    'data_type,' + \
    'data_format,' + \
    'data_category,' + \
    'experimental_strategy,' + \
    'submitter_id,' + \
    'analysis.workflow_type,' + \
    'cases.samples.portions.analytes.aliquots.aliquot_id,' + \
    'cases.samples.portions.analytes.aliquots.submitter_id,' + \
    'cases.project.project_id,' + \
    'cases.project.program.name'

# GDC API parameters for HTTPS request.
# added expand so entries get populated even without end data
DATAFILE_PARAMS: Dict[str, str] = {
    'pretty': 'true',
    'size': '50',
    'filters': json.dumps(DATAFILE_FILTER),
    'fields': DATAFILE_FIELDS,
    'expand': 'cases.samples.portions,cases.samples.portions.analytes,cases.samples.portions.analytes.aliquots',
    'from': '0'
}

DATAFILE_HEADERS: List[str] = [
    'data_type_display',
    'data_format',
    'experimental_strategy',
    'workflow',
    'project',
    'program',
]


# pylint: disable=too-many-instance-attributes,too-many-arguments,too-few-public-methods
class GdcApiDatafile(GdcApiHistoryMixin):
    """
    Datafile files from GDC, have data_format, project, and program data at this level.
    Uses history details inherited from GdcApiHistoryMixin.
    data_format - BCR XML
    program - top level descriptor, like TCGA, TARGET, etc.
    project - similar to disease type. Usually, program plus disease or cohort
    experimental_strategy - strings such as Methylation Array, WGS, RNA-Seq, WXS
    data_type_display - strings such as Methylation Beta Value, Gene Level Copy Number, Gene Expression Quantification, Masked Somatic Mutation
    workflow - strings such as SeSAMe Methylation Beta Estimation, AscatNGS, STAR - Counts, Aliquot Ensemble Somatic Variant Merging and Masking
    data_format - currently, object are TXT, TSV, MAF
    .
    Data hierarchy order for menus:
    attribute ----> menu entry
    ==============================
    program ----> project
    project ----> sub-project
    experimental_strategy ----> category
    data_type_display ----> platform
    workflow ----> data
    history_version ----> Use as new "data-version" entry?
    """
    # declare but do not set member attributes
    source: str
    data_type_display: str
    data_format: str
    experimental_strategy: str
    workflow: str
    project: str
    program: str
    sample_dict: Dict[str, GdcApiSample]
    _dirty_samples: Optional[bool]

    # pylint: disable=too-many-nested-blocks,too-many-branches,too-many-locals
    def __init__(self: 'GdcApiDatafile', the_dict: Dict, the_json_flag: bool) -> None:
        """
        Initialize GdcApiDatafile from JSON or TSV built Dictionary.
        JSON data will be for non-History labeled attributes.
        :param the_dict: dictionary with attribute values
        :param the_json_flag: if True Dictionary is from JSON, otherwise from TSV
        """
        # the_uuid: str, the_dict: Dict, the_json_flag: bool
        super().__init__(the_dict, the_json_flag)
        # declare instance variables
        # DO NOT DECLARE ones in parent class
        # this hardcoding is safe as long as it is only used for GDC download/convert
        self.source: str = 'GDC'
        self.data_type_display: str = ""
        self.data_format: str = ""
        self.experimental_strategy: str = ""
        self.workflow: str = ""
        self.project: str = ""
        self.program: str = ""
        # when reading from disk, sample_dict is not populated until it is needed
        # since the thousands of sample files otherwise take too long to load
        # (too long being 12 hours and it was still in the middle of TCGA files)
        self.sample_dict: Dict[str, GdcApiSample] = {}
        # internal flag for sample list being dirty and needing writing
        # True - dirty and needs to be written
        # False - clean and loaded from TSV file
        # None - clean and not loaded from TSV file
        self._dirty_samples: Optional[bool] = None
        # set instance variable values
        if the_json_flag:
            # set samples to dirty, since JSON means information is from GDC and therefore new
            self._dirty_samples = True
            self.data_type_display = the_dict['data_type']
            self.data_format = the_dict['data_format']
            self.experimental_strategy = the_dict['experimental_strategy']
            case_entry: List[Dict] = the_dict['cases']
            case_dict: Dict
            # use cases to find project and program
            for case_dict in case_entry:
                if '' == self.project:
                    self.project = case_dict['project']['project_id']
                    self.program = case_dict['project']['program']['name']
            # back to top level for RPPA vs other samples
            if "Protein Expression Quantification" == self.data_type_display:
                self.workflow = "Protein Analysis"
                rppa_barcode = the_dict['submitter_id']
                # remove _RPPA off end
                rppa_barcode = rppa_barcode.removesuffix("_RPPA")
                sample: GdcApiSample = GdcApiSample(None, rppa_barcode, rppa_barcode)
                self.sample_dict[sample.sample_uuid] = sample
            else:
                self.workflow = the_dict['analysis']['workflow_type']
                # use cases again to go through cases down to aliquot (sample)
                for case_dict in case_entry:
                    sample_entry_list: List[Dict] = case_dict['samples']
                    samples_dict: Dict
                    for samples_dict in sample_entry_list:
                        portion_list: List[Dict] = samples_dict['portions']
                        portion: Dict
                        for portion in portion_list:
                            if 'analytes' in portion:
                                analyte_list: List[Dict] = portion['analytes']
                                analyte: Dict
                                for analyte in analyte_list:
                                    aliquot_list: List[Dict] = analyte['aliquots']
                                    aliquot: Dict
                                    for aliquot in aliquot_list:
                                        sample_uuid: str = aliquot['aliquot_id']
                                        sample_barcode: str = aliquot['submitter_id']
                                        sample: GdcApiSample = GdcApiSample(None, sample_barcode, sample_uuid)
                                        self.sample_dict[sample.sample_uuid] = sample
        else:
            # set samples to None, since TSV means information is from TSV file and therefore already downloaded
            # but we are not loading the data yet
            self._dirty_samples = None
            entry: str
            for entry in DATAFILE_HEADERS:
                setattr(self, entry, the_dict.get(entry))
    # pylint: enable=too-many-nested-blocks,too-many-branches,too-many-locals

    def write_datafile(self: 'GdcApiDatafile', the_out_file: io.TextIOWrapper) -> None:
        """
        Write out just the GdcApiDatafile portion of a file element to a file,
        after first calling parent history to write out its portion.
        Writes a newline.
        :param the_out_file: stream to which to write strings
        :return: None
        """
        data_list: List[str] = []
        for attr in DATAFILE_HEADERS:
            data_list.append(getattr(self, attr))
        self.write_history(the_out_file)
        write_tsv_list(the_out_file, data_list, True, False)

    def write_sample_index(self: 'GdcApiDatafile', the_sample_dir: str) -> None:
        """
        Write the sample index file, only if this object is marked as dirty.
        :param the_sample_dir: Full path to base sample index output directory
        :return: Nothing
        """
        # only write if samples is dirty
        # False and None will both not write the index
        if self._dirty_samples:
            index_sub_dir: str = os.path.join(the_sample_dir, self.program, self.project)
            os.makedirs(index_sub_dir, exist_ok=True)
            index_file: str = os.path.join(index_sub_dir, self.uuid + ".tsv")
            write_sample_index_file(index_file, list(self.sample_dict.values()))
            self._dirty_samples = False

    def get_download_file(self: 'GdcApiDatafile', the_base_dir: str) -> str:
        """
        Build the "datafile" version of a file path.
        Builds /base-dir/program/project/experimental_strategy/data_type_display/workflow/version/file.name path.
        :param the_base_dir: full path to base directory to add specialized path and filename
        :return: full path with filename to which to download the file
        """
        # print("GdcApiDatafile::get_download_file", flush=True)
        return os.path.join(the_base_dir,
                            self.program, self.project, self.experimental_strategy,
                            self.data_type_display, self.workflow, self.history_version,
                            self.file_name)

    def get_file_archive(self: 'GdcApiDatafile', the_base_dir: str) -> str:
        """
        Build the "datafile" version of a file path for the ZIP archive.
        Builds /base-dir/program/project/experimental_strategy/data_type_display/workflow/version/uuid.zip
        :param the_base_dir: full path to base directory to add specialized path and filename
        :return: full path with filename which is the ZIP archive
        """
        # print("GdcApiDatafile::get_download_file", flush=True)
        return os.path.join(the_base_dir,
                            self.program, self.project, self.experimental_strategy,
                            self.data_type_display, self.workflow, self.history_version,
                            self.uuid + ".zip")

    def read_sample_file(self: 'GdcApiDatafile', the_sample_dir: str) -> None:
        """
        Read the sample file and populate sample_dict for this object
        :param the_sample_dir: full path to base directory to add specialized path and filename
        :return: Nothing
        """
        if self._dirty_samples is None:
            # set dirty flag to clean and loaded (False)
            self._dirty_samples = False
            index_sub_dir: str = os.path.join(the_sample_dir, self.program, self.project)
            index_file: str = os.path.join(index_sub_dir, self.uuid + ".tsv")
            read_sample_index_file(index_file, self.sample_dict)

    def find_tumor_barcode(self: 'GdcApiDatafile', the_batches: pandas.DataFrame) -> str:
        """
        Find barcode in sample_dict with sample type code 01
        Return corresponding barcode or empty string.
        :param the_batches: DataFrame of batches.tsv file, for determining type based on sample id(s).
        :return: Return corresponding barcode or empty string.
        """
        barcode: str = ""
        sample: GdcApiSample
        for sample in self.sample_dict.values():
            is_tumor: bool = self.is_a_tumor_barcode(sample.sample_uuid, sample.sample_barcode, the_batches)
            if is_tumor:
                barcode = sample.sample_barcode
        if '' == barcode:
            # list_control: List[str] = ['20', '21', '22', '23', '24', '25', '26', '27', '28', '29' ]
            # list_normal: List[str] = ['10', '11', '12', '13', '14', '15', '16', '17', '18', '19' ]
            list_tumor: List[str] = ['01', '02', '03', '04', '05', '06', '07', '08', '09',
                                     '1', '2', '3', '4', '5', '6', '7', '8', '9',
                                     '85', '86']
            for sample in self.sample_dict.values():
                if '' == barcode:
                    print(f"Attempt manual parsing of barcode not in biospecimen files {sample.sample_barcode} with UUID {sample.sample_uuid}", flush=True)
                    # see if it can be manually parsed
                    barcode_split: List[str] = sample.sample_barcode.split("-")
                    if (len(barcode_split) == 7) & ('TCGA' == barcode_split[0]):
                        # TCGA Barcode
                        type_str: str = barcode_split[3]
                        type_str = type_str[:-1]
                        if type_str in list_tumor:
                            barcode = sample.sample_barcode
                        add_error(f"TCGA Barcode not in biospecimen files {sample.sample_barcode} with UUID {sample.sample_uuid}")
                    elif (len(barcode_split) == 8) & ('HCM' == barcode_split[0]):
                        # HCMI barcode - most likely to need this processing
                        # found barcodes in GDC but not in biospecimen files
                        type_str: str = barcode_split[4]
                        type_str = type_str[:-1]
                        if type_str in list_tumor:
                            barcode = sample.sample_barcode
                        add_error(f"HCMI Barcode not in biospecimen files {sample.sample_barcode} with UUID {sample.sample_uuid}")
        assert '' != barcode, f"Tumor Barcode not found for {self.get_file_archive('/')}"
        return barcode

    # pylint: disable=too-many-branches
    def is_a_tumor_barcode(self: 'GdcApiDatafile', the_uuid: str, the_barcode: str, the_batches: pandas.DataFrame) -> bool:
        """
        Determine if barcode represents a tumor by looking up values in the batch information
        :param the_uuid: UUID of sample/aliquot
        :param the_barcode: Barcode of the sample/aliquot
        :param the_batches: Batch information
        :return: True is barcode is found as a tumor
        """
        is_tumor: bool = False
        # check if sample_type_name and sample_type_id are found
        check_type_name: bool = False
        check_type_id: bool = False
        if 'sample_type_name' in the_batches.columns:
            check_type_name = True
        if 'sample_type_id' in the_batches.columns:
            check_type_id = True

        assert 'Unknown' != the_barcode, f"Unknown Barcode sent for {self.get_file_archive('/')}"
        matched_row: pandas.DataFrame = the_batches.loc[the_batches['aliquot_barcode'].isin([the_barcode])]
        my_type_a: typing.Optional[str] = None
        my_type_b: typing.Optional[str] = None
        if check_type_name:
            if matched_row['sample_type_name'].size > 0:
                my_type_a = matched_row['sample_type_name'].values[0].lower()
        if check_type_id:
            if matched_row['sample_type_id'].size > 0:
                my_type_b = matched_row['sample_type_id'].values[0].lower()
        if check_for_tumor(my_type_a, my_type_b):
            is_tumor = True
        if not is_tumor:
            assert 'Unknown' != the_uuid, f"Unknown UUID sent for {self.get_file_archive('/')}"
            matched_row: pandas.DataFrame = the_batches.loc[the_batches['aliquot_uuid'].isin([the_uuid.upper()])]
            if matched_row.size < 1:
                # most of the time biospecimen files are upper case.
                # sometimes they are lower.
                matched_row = the_batches.loc[the_batches['aliquot_uuid'].isin([the_uuid.lower()])]
            my_type_a = None
            my_type_b = None
            if check_type_name:
                if matched_row['sample_type_name'].size > 0:
                    my_type_a = matched_row['sample_type_name'].values[0].lower()
            if check_type_id:
                if matched_row['sample_type_id'].size > 0:
                    my_type_b = matched_row['sample_type_id'].values[0].lower()
            if check_for_tumor(my_type_a, my_type_b):
                is_tumor = True
        return is_tumor
    # pylint: disable=too-many-branches

    def check_filenames_for_uuid(self: 'GdcApiDatafile') -> str:
        """
        Check file_name for uuid from sample_dict.
        Return corresponding barcode or empty string.
        :return: Return corresponding barcode or empty string.
        """
        barcode: str = ""
        sample: GdcApiSample
        for sample in self.sample_dict.values():
            if sample.sample_uuid in self.file_name:
                barcode = sample.sample_barcode
        return barcode
# pylint: enable=too-many-instance-attributes,too-many-arguments,too-few-public-methods


# pylint: disable=too-many-branches
def check_for_tumor(the_val_a: typing.Optional[str], the_val_b: typing.Optional[str]) -> bool:
    """
    Value a and b have had .lower() called before being sent here. May be None.
    :param the_val_a: usually sample_type_name - used to determine if it represents a tumor.
    :param the_val_b: usually sample_type_id - used to determine if it represents a tumor.
    :return: True if this looks like a tumor.
    """
    normal: bool = False
    tumor: bool = False
    control: bool = False
    if the_val_a is not None:
        if 'normal' in the_val_a:
            normal = True
        if 'granulocytes' in the_val_a:
            normal = True
        if 'tumor' in the_val_a:
            tumor = True
        if 'control' in the_val_a:
            control = True
        if 'metastatic' == the_val_a:
            tumor = True
        if '01' == the_val_a:
            tumor = True
        if '1' == the_val_a:
            tumor = True
    if (not normal) & (not tumor) & (not control):
        if the_val_b is not None:
            # if 'normal' in the_val_b:
            #     normal = True
            # if 'granulocytes' in the_val_b:
            #     normal = True
            # if 'control' in the_val_b:
            #     control = True
            if 'tumor' in the_val_b:
                tumor = True
            if 'metastatic' == the_val_b:
                tumor = True
            if '85' == the_val_b:
                tumor = True
            if '86' == the_val_b:
                tumor = True
            if '01' == the_val_b:
                tumor = True
            if '06' == the_val_b:
                tumor = True
            if '1' == the_val_b:
                tumor = True
            if '31' == the_val_b:
                tumor = True
    # do not check for problems here, since we parse barcodes outside this
    return tumor
# pylint: enable=too-many-branches

# ########################################################
# retrieve data from GDC HTTP API
# ########################################################


def get_datafile_list_from_gdc(the_entries: Dict[str, GdcApiDatafile]) -> List[GdcApiDatafile]:
    """
    Call the file endpoint, looping until we get no results on the new page, populating
    list of GdcApiDatafile. Use existing object from dictionary if present.
    :param the_entries: Dictionary of GdcApiDatafile objects.
    :return: Updated list of GdcApiDatafile objects.
    """
    results: List[GdcApiDatafile] = []
    # development max size limit (0 does not limit download)
    max_size: int = GLOBAL_MAX_COUNT
    size: int = GLOBAL_DLD_COUNT
    from_int: int = 0
    continue_while: bool = True
    print("get_file_list_from_gdc start", flush=True)
    while continue_while:
        print(f"loop size={size} from_int={from_int}", flush=True)
        DATAFILE_PARAMS['size'] = f'{size}'
        DATAFILE_PARAMS['from'] = f'{from_int}'
        response: requests = requests.get(DATAFILE_ENDPOINT, params=DATAFILE_PARAMS, timeout=60, allow_redirects=True)
        print(f"get_file_list_from_gdc response={response}", flush=True)
        my_str: str = json.dumps(response.json(), indent=2)
        # print(f"get_file_list_from_gdc my_str={my_str}", flush=True)
        my_dict: Dict = json.loads(my_str)
        # print(f"get_file_list_from_gdc my_dict={my_dict}", flush=True)
        datafile_list: List[Dict] = my_dict['data']['hits']
        my_entry: Dict
        for my_entry in datafile_list:
            new_item: GdcApiDatafile = GdcApiDatafile(my_entry, True)
            if new_item.uuid in the_entries:
                # if existing entry found, use existing
                new_item = the_entries[new_item.uuid]
            results.append(new_item)
        len_datafile_list: int = len(datafile_list)
        print(f"get_file_list_from_gdc len(len_datafile_list)={len_datafile_list}", flush=True)
        if len_datafile_list < 1:
            continue_while = False
        # development short version
        if max_size > 0:
            if len(results) >= max_size:
                continue_while = False
        from_int += size
    # end of loop
    print("get_file_list_from_gdc done", flush=True)
    return results

# ########################################################
# writing files
# ########################################################


def write_headers_datafiles(the_out_file: io.TextIOWrapper) -> None:
    """
    Write the headers for datafile portion of object,
    after calling to write History portion.
    Writes newline at the end.
    :param the_out_file: stream to which to write headers
    :return: None
    """
    write_headers_history(the_out_file)
    write_tsv_list(the_out_file, DATAFILE_HEADERS, True, False)


def write_datafile_index(the_file_datafiles: str, the_sample_dir: str, the_entries: List[GdcApiDatafile]) -> None:
    """
    Write the headers and line items for index files.
    Replaces existing index file.
    :param the_file_datafiles: full path to datafile index file to write.
    :param the_sample_dir: full path to directory for index files
    :param the_entries: list of GdcApiDatafile to write to index.
    :return: None
    """
    print(f"update_file_datafiles the_file_datafiles={the_file_datafiles}", flush=True)
    the_entries.sort(key=lambda my_datafile: (my_datafile.program, my_datafile.project,
                                              my_datafile.experimental_strategy, my_datafile.data_type_display,
                                              my_datafile.workflow, my_datafile.history_version), reverse=False)
    out_file: io.TextIOWrapper
    with open(the_file_datafiles, 'w', encoding='utf-8') as out_file:
        write_headers_datafiles(out_file)
        my_entry: GdcApiDatafile
        for my_entry in the_entries:
            my_entry.write_datafile(out_file)
            my_entry.write_sample_index(the_sample_dir)

# ########################################################
# reading files
# ########################################################


def read_map_datafiles_by_program_projects(the_file_datafiles: str) -> Dict:
    """
    Read datafiles index and return dictionary, where key is a tuple
    of the Program and Project strings, and the value is List of GdcApiDatafile.
    :param the_file_datafiles: full path and filename for datafiles index file.
    :return: Dictionary with key as tuple of program and project. Value is List of GdcApiDatafile.
    """
    datafiles: Dict[str, GdcApiDatafile] = read_datafile_index(the_file_datafiles)
    mapping_datafiles: Dict[(str, str), List[GdcApiDatafile]] = {}
    datafile: GdcApiDatafile
    for datafile in datafiles.values():
        # build program/project tuple for key
        key: (str, str) = (datafile.program, datafile.project)
        list_of_datafiles: List[GdcApiDatafile] = mapping_datafiles.get(key, [])
        list_of_datafiles.append(datafile)
        mapping_datafiles[key] = list_of_datafiles
    return mapping_datafiles


def read_datafile_for_uuid_barcodes(the_datafile_list: List[GdcApiDatafile], the_sample_dir: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Read datafiles index, and associated samples, and populate two dictionaries.
    One dictionary maps UUIDs to Barcodes.
    The other dictionary maps Barcodes to UUIDs.
    Necessary for some conversion processes.
    :param the_datafile_list: List of GdcApiDatafile objects (REQUIRED for some conversions)
    :param the_sample_dir: full path to sample dir, which are associated with the dataset index (REQUIRED for some conversions)
    :return: Tuple of mapping dictionaries. First dictionary maps uuid to barcode. Second dictionary maps barcode to uuid.
    """
    mapping_uuid_to_barcode: Dict[str, str] = {}
    mapping_barcode_to_uuid: Dict[str, str] = {}
    datafile: GdcApiDatafile
    for datafile in the_datafile_list:
        # read datafile sample information for mapping (if needed)
        datafile.read_sample_file(the_sample_dir)
        # build program/project tuple for key
        sample: GdcApiSample
        for sample in datafile.sample_dict.values():
            mapping_uuid_to_barcode[sample.sample_uuid] = sample.sample_barcode
            mapping_barcode_to_uuid[sample.sample_barcode] = sample.sample_uuid
    return mapping_uuid_to_barcode, mapping_barcode_to_uuid


def read_datafile_index(the_file_datafiles: str) -> Dict[str, GdcApiDatafile]:
    """
    Read index file and create Dictionary of GdcApiDatafile object.
    Do not read sample info yet, as it takes too long to load.
    :param the_file_datafiles: full path and filename of datafile index file to read.
    :return: Dictionary with uuid as keys and GdcApiDatafile objects as values.
    """
    print(f"read_file_datafiles the_file_datafiles={the_file_datafiles}", flush=True)
    in_file: io.TextIOWrapper
    my_entries: Dict[str, GdcApiDatafile] = {}
    with open(the_file_datafiles, 'r', encoding='utf-8') as in_file:
        keys: List[str] = read_headers(in_file)
        line: str = in_file.readline().rstrip('\n')
        index: int = 1
        while '' != line:
            if 0 == index % 1000:
                print(f"read_file_datafiles datafile index={index}", flush=True)
            index += 1
            values: List[str] = line.split("\t")
            tsv_dict: Dict[str, str] = dict(zip(keys, values))
            entry: GdcApiDatafile = GdcApiDatafile(tsv_dict, False)
            # do not read sample file here, read it only when needed
            # too many files to pre-read
            my_entries[entry.uuid] = entry
            line = in_file.readline().rstrip('\n')
    return my_entries

# ########################################################
#
# ########################################################
