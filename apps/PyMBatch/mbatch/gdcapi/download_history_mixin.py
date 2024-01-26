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
import os
import json
import gzip
import shutil
import hashlib
from typing import List, Dict
import time
import zipfile
import requests
from mbatch.test.common import write_tsv_list
from mbatch.test.common import add_error

# ########################################################
# History Endpoint, getting history data for a data entry
# ########################################################

# GDC API endpoint for history API
HISTORY_ENDPOINT: str = 'https://api.gdc.cancer.gov/history/'

# GDC API list of fields to include in results
HISTORY_FIELDS: str = '' + \
    'uuid,' + \
    'version,' + \
    'file_change,' + \
    'release_date,' + \
    'data_release'

# GDC API parameters for HTTPS request.
# added expand so entries get populated even without end data
HISTORY_PARAMS: Dict[str, str] = {
    'pretty': 'true',
    'size': '50',
    'fields': HISTORY_FIELDS
}

HISTORY_HEADERS: List[str] = [
    'uuid',
    'md5sum',
    'file_name',
    'history_version',
    'history_status',
    'history_release_date',
    'history_data_release'
]


# pylint: disable=too-many-instance-attributes,too-many-arguments,too-few-public-methods
class GdcApiHistoryMixin:
    """
    Top level class that gets added to clinical, datafile, and biospecimen classes.
    GdcApiHistoryMixin implements shared functionality around file information (uuid, md5, and filename),
    and history data related to release and versioning status and numbers.
    uuid - GDC UUID for file endpoint
    md5sum - md5sum for downloaded file
    file_name - GDC name for this file, with extension. This needs to be retained, because
                it may be needed for processing some files (like RPPA).
    history_version - integer (stored as string) version of release for this file, with 1 as older and larger numbers being newer.
    history_status - values like released (current release), superseded (replaced by newer release). There may be others.
    history_release_date - YYYY-MM-DD release, corresponding to history_data_release
    history_data_release - "xx.0" release number string, with larger numbers being newer.
    In case of dataset wierdnesses, double check the data_release only labels new files in dataset (I think it does)
    This means that a dataset will have some release N (for newest) and some release N-X (for earlier release) data files in it
    """
    # declare but do not set member attributes
    uuid: str
    md5sum: str
    file_name: str
    history_version: str
    history_status: str
    history_release_date: str
    history_data_release: str

    def __init__(self: 'GdcApiHistoryMixin', the_dict: Dict, the_json_flag: bool) -> None:
        """
        Initialize GdcApiHistoryMixin from JSON or TSV built Dictionary.
        JSON data will be for non-History labeled attributes.
        :param the_dict: dictionary with attribute values
        :param the_json_flag: if True Dictionary is from JSON, otherwise from TSV
        """
        super().__init__()
        # declare instance variables
        self.uuid: str = ""
        self.md5sum: str = ""
        self.file_name: str = ""
        self.history_version: str = ""
        self.history_status: str = ""
        self.history_release_date: str = ""
        self.history_data_release: str = ""
        # set instance variables
        if the_json_flag:
            self.uuid = the_dict.get('file_id')
            self.md5sum = the_dict.get('md5sum')
            self.file_name = the_dict.get('file_name')
            self.history_version = ""
            self.history_status = ""
            self.history_release_date = ""
            self.history_data_release = ""
        else:
            entry: str
            for entry in HISTORY_HEADERS:
                setattr(self, entry, the_dict.get(entry))

    def update_history(self: 'GdcApiHistoryMixin', the_history: List[Dict]) -> bool:
        """
        Called from update_history_from_gdc - used to update History related attributes
        from JSON (using history endpoint).
        :param the_history: Dictionary build from /history JSON with history attribute data
        :return: True is data was changed or set
        """
        changed: bool = False
        found: bool = False
        my_dict: Dict[str, str]
        for my_dict in the_history:
            # only one dictionary should match UUID
            if (not changed) & (not found):
                if my_dict['uuid'] == self.uuid:
                    # only one dictionary should match UUID
                    found = False
                    if my_dict['version'] != self.history_version:
                        changed = True
                        self.history_version = my_dict['version']
                    if my_dict['file_change'] != self.history_status:
                        changed = True
                        self.history_status = my_dict['file_change']
                    if my_dict['data_release'] != self.history_data_release:
                        changed = True
                        self.history_data_release = my_dict['data_release']
                    if my_dict['release_date'] != self.history_release_date:
                        changed = True
                        self.history_release_date = my_dict['release_date']
        return changed

    def write_history(self: 'GdcApiHistoryMixin', the_out_file: io.TextIOWrapper) -> None:
        """
        Write out just the history portion of a file elment to a file.
        Adds a trailing tab and does NOT write a newline.
        :param the_out_file: stream to which to write strings
        :return: None
        """
        data_list: List[str] = []
        for attr in HISTORY_HEADERS:
            data_list.append(getattr(self, attr))
        write_tsv_list(the_out_file, data_list, False, True)

    def raise_error(self: 'GdcApiHistoryMixin', the_message: str) -> None:
        """
        Raise an error -- done from function to avoid no return value or
        unreachable code problems with "abstract" functions.
        :param the_message: message to raise
        :return: nothing
        """
        raise TypeError(f"Processing {the_message} as GdcApiHistoryMixin, which should not happen")

    def get_download_file(self: 'GdcApiHistoryMixin', the_base_dir: str) -> str:
        """
        Build the "history only" version of a file path. Should never be called,
        provided for override purposes.
        :param the_base_dir: full path to base directory to add specialized path and filename
        :return: (nothing, raise error) full path with filename to which to download the file
        """
        print("GdcApiHistoryMixin::get_download_file", flush=True)
        self.raise_error("Processing get_download_file as GdcApiHistoryMixin, which should not happen")
        return os.path.join(the_base_dir, self.history_version, self.file_name)

    def get_file_archive(self: 'GdcApiHistoryMixin', the_base_dir: str) -> str:
        """
        Build the "history only" version of a file path for the ZIP archive.
        Should never be called, provided for override purposes.
        :param the_base_dir: full path to base directory to add specialized path and filename
        :return: (nothing, raise error) full path with filename which is the ZIP archive
        """
        print("GdcApiHistoryMixin::get_download_file", flush=True)
        self.raise_error("Processing get_file_archive as GdcApiHistoryMixin, which should not happen")
        return os.path.join(the_base_dir, self.history_version, self.uuid + ".zip")

    def check_md5_sum(self: 'GdcApiHistoryMixin', the_file: str) -> bool:
        """
        Check the given file's md5 sum against the one provided
        :param the_file: full path to filename to check md5sum
        :return: True if md5 sum matches
        """
        passed: bool = False
        hash_md5: hashlib.md5 = hashlib.md5()
        sum_file: io.TextIOWrapper
        with open(the_file, "rb") as sum_file:
            for chunk in iter(lambda: sum_file.read(4096), b""):
                hash_md5.update(chunk)
        calc: str = hash_md5.hexdigest()
        if calc == self.md5sum:
            passed = True
        else:
            print(f"check md5sum failed for {self.uuid} file {the_file}", flush=True)
        return passed

    def download_gdcapi_file(self: 'GdcApiHistoryMixin', the_base_dir: str) -> bool:
        """
        Download self file, un-gzip if needed, then ZIP the file.
        Filename is GDC file name. ZIP archive name is uuid.zip.
        This makes finding the archive simpler, while sometimes filename
        is needed for processing the data. (RPPA data has sample id in filename.)
        :param the_base_dir: base directory to which to add more path and download the filename.
        :return: True is file was downloaded and MD5 matched.
        """
        success: bool = False
        my_file: str = self.get_download_file(the_base_dir)
        trimmed_file: str = my_file
        zip_file: str = os.path.join(os.path.dirname(trimmed_file), self.uuid + ".zip")
        if os.path.exists(zip_file):
            print(f"download_gdcapi_file ZIP exists {zip_file}", flush=True)
            success = True
        else:
            if trimmed_file.endswith(".gz"):
                trimmed_file: str = trimmed_file.removesuffix('.gz')
            if os.path.exists(trimmed_file):
                print(f"download_gdcapi_file exists {trimmed_file}", flush=True)
            else:
                print(f"download_gdcapi_file download {my_file}", flush=True)
                response: requests
                with requests.get('https://api.gdc.cancer.gov/data/' + self.uuid, stream=True, timeout=60, allow_redirects=True) as response:
                    if response.ok:
                        out_file: io.TextIOWrapper
                        success = True
                        new_dir: str = os.path.dirname(my_file)
                        if not os.path.exists(new_dir):
                            os.makedirs(new_dir, exist_ok=True)
                        with open(my_file, 'wb') as out_file:
                            shutil.copyfileobj(response.raw, out_file)
                if success:
                    # print(f"download_gdcapi_file check_md5_sum {my_file}", flush=True)
                    if not self.check_md5_sum(my_file):
                        add_error(f"download_gdcapi_file failed {my_file}")
                        os.unlink(my_file)
                        success = False
                    else:
                        # print(f"download_gdcapi_file passed {my_file}", flush=True)
                        if my_file.endswith(".gz"):
                            out_file: io.TextIOWrapper
                            in_file: io.TextIOWrapper
                            # print(f"download_gdcapi_file gz {trimmed_file}", flush=True)
                            with open(trimmed_file, 'wb') as out_file:
                                with gzip.open(my_file, 'rb') as in_file:
                                    shutil.copyfileobj(in_file, out_file)
                            os.unlink(my_file)
                        # zip check_file into zip_file
                        my_zip: zipfile.ZipFile
                        with zipfile.ZipFile(zip_file, 'w', compression=zipfile.ZIP_DEFLATED) as my_zip:
                            my_zip.write(trimmed_file, arcname=os.path.basename(trimmed_file))
                        os.unlink(trimmed_file)
        return success
# pylint: enable=too-many-instance-attributes,too-many-arguments,too-few-public-methods

# ########################################################
# retrieve data from GDC HTTP API
# ########################################################


# pylint: disable=too-many-nested-blocks,too-many-branches,broad-except
def update_history_from_gdc(the_entries: List[GdcApiHistoryMixin], the_depth: int) -> bool:
    """
    Loop through list of GdcApiHistoryMixin objects (inherited by GdcApiClinical, Datafile, and Biospecimen),
    and for entries without any history_X data or where the history_version is 'released',
    query the /history endpoint and update its data. (We update released for when it gets changed to
    superseded or some other value.
    :param the_entries: List of GdcApiHistoryMixin derived objects.
    :param the_depth: Number of times to retry failed updates
    :return: True if any elements were updated
    """
    print(f"update_history_from_gdc start the_depth={the_depth}", flush=True)
    changed: bool = False
    if the_depth > 0:
        # retry list for 404, etc
        retry_list: List[GdcApiHistoryMixin] = []
        # iterations
        total: int = len(the_entries)
        index: int
        my_entry: GdcApiHistoryMixin
        print("update_history_from_gdc iterate start", flush=True)
        for index, my_entry in enumerate(the_entries):
            if ('released' == my_entry.history_version) | ('' == my_entry.history_version):
                print(f"update_history_from_gdc UPDATE index={index+1}/{total} my_entry.uuid={my_entry.uuid}", flush=True)
                # GDC API filter for release open-access data. JSON format.
                # cannot filter history for just one UUID
                # does not use timeout by default (60 seconds)
                try:
                    # sleep for 10 seconds every 100 request to prevent overloading GDC
                    if 0 == index % 100:
                        time.sleep(10)
                        print(f"update_history_from_gdc index={index} from {len(the_entries)}", flush=True)
                    response: requests = requests.get(HISTORY_ENDPOINT + my_entry.uuid, params=HISTORY_PARAMS,
                                                      timeout=60, allow_redirects=True)
                    # print(f"update_history_from_gdc response={response}", flush=True)
                    if response.ok:
                        my_str: str = json.dumps(response.json(), indent=2)
                        # print(f"update_history_from_gdc my_str={my_str}", flush=True)
                        my_dict_list: List[Dict] = json.loads(my_str)
                        if my_entry.update_history(my_dict_list):
                            changed = True
                    else:
                        # print(f"update_history_from_gdc error={response.status_code} my_entry.uuid={my_entry.uuid}", flush=True)
                        retry_list.append(my_entry)
                except Exception as my_except:
                    print(f"update_history_from_gdc exception {my_except}", flush=True)
                    print(my_except, flush=True)
                    retry_list.append(my_entry)
            else:
                # do not need to update as it is old
                print(f"update_history_from_gdc {my_entry.history_version} does not need update index={index + 1}/{total} my_entry.uuid={my_entry.uuid}", flush=True)
        print(f"update_history_from_gdc iterate done - check for errors n={len(retry_list)}", flush=True)
        if len(retry_list) > 0:
            time.sleep(60)
            if update_history_from_gdc(retry_list, the_depth-1):
                changed = True
    else:
        for index, my_entry in enumerate(the_entries):
            add_error(f"update_history_from_gdc Getting history error for my_entry.uuid={my_entry.uuid}")
    print("update_history_from_gdc done", flush=True)
    return changed
# pylint: enable=too-many-nested-blocks,too-many-branches,broad-except

# ########################################################
# writing files
# ########################################################


def write_headers_history(the_out_file: io.TextIOWrapper) -> None:
    """
    Write the headers for history portion of object.
    Called by inheritors before they write their headers.
    Adds tab afterward, and does NOT add newline
    :param the_out_file: stream to which to write headers
    :return: None
    """
    write_tsv_list(the_out_file, HISTORY_HEADERS, False, True)


# ########################################################
#
# ########################################################
