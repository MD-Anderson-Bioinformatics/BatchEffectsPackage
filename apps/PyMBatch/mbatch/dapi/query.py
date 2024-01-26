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


from typing import List, Dict, Tuple, BinaryIO
import json
import requests


# pylint: disable=too-many-instance-attributes,too-few-public-methods
class DapiEntry:
    request_data: str
    request_results: str
    request_view: str
    dataset_id: str
    source: str
    program: str
    project: str
    category: str
    platform: str
    data: str
    details: str
    data_version: str
    test_version: str
    dataset_type: str

    def __init__(self: 'DapiQuery', the_data_entry: List[str]) -> None:
        actions: str = the_data_entry[0]
        splitted: List[str] = actions.split(" | ")
        self.request_data = splitted[0]
        self.request_results = splitted[1]
        self.request_view = splitted[2]
        self.dataset_id = the_data_entry[1]
        self.source = the_data_entry[2]
        self.program = the_data_entry[3]
        self.project = the_data_entry[4]
        self.category = the_data_entry[5]
        self.platform = the_data_entry[6]
        self.data = the_data_entry[7]
        self.details = the_data_entry[8]
        self.data_version = the_data_entry[9]
        self.test_version = the_data_entry[10]
        self.dataset_type = the_data_entry[11]

    def print_status(self: 'DapiEntry', the_summary: bool) -> None:
        if the_summary:
            print(f"DapiEntry {self.source} {self.project} {self.data} {self.details} {self.data_version} {self.test_version} {self.dataset_type}", flush=True)
        else:
            print(f"DapiEntry self.request_data={self.request_data}", flush=True)
            print(f"DapiEntry self.request_results={self.request_results}", flush=True)
            print(f"DapiEntry self.request_view={self.request_view}", flush=True)
            print(f"DapiEntry self.dataset_id={self.dataset_id}", flush=True)
            print(f"DapiEntry self.source={self.source}", flush=True)
            print(f"DapiEntry self.program={self.program}", flush=True)
            print(f"DapiEntry self.project={self.project}", flush=True)
            print(f"DapiEntry self.category={self.category}", flush=True)
            print(f"DapiEntry self.platform={self.platform}", flush=True)
            print(f"DapiEntry self.data={self.data}", flush=True)
            print(f"DapiEntry self.details={self.details}", flush=True)
            print(f"DapiEntry self.data_version={self.data_version}", flush=True)
            print(f"DapiEntry self.test_version={self.test_version}", flush=True)
            print(f"DapiEntry self.dataset_type={self.dataset_type}", flush=True)
# pylint: enable=too-many-instance-attributes,too-few-public-methods


# pylint: disable=too-many-arguments,too-many-instance-attributes,too-few-public-methods
class DapiQuery:
    """
    Store information on DAPI Query results
    """
    # do not set method variables, as they should be initialized in the init function
    # #####
    # URL
    api_url_base: str
    api_url_query: str
    api_url_index: str
    api_url_dsblob: str
    api_url_blobdata: str
    # AVAILABLE VALUES
    available_jobtype: List[str]
    available_platforms: List[str]
    available_categories: List[str]
    available_projects: List[str]
    available_details: List[str]
    available_dataversions: List[str]
    available_testversions: List[str]
    available_program: List[str]
    available_data: List[str]
    available_sources: List[str]
    # SELECTED VALUES
    selected_jobtype: List[str]
    selected_platforms: List[str]
    selected_categories: List[str]
    selected_projects: List[str]
    selected_details: List[str]
    selected_dataversions: List[str]
    selected_testversions: List[str]
    selected_program: List[str]
    selected_data: List[str]
    selected_sources: List[str]
    # DATASETS
    available_datasets: List[DapiEntry]

    def __init__(self: 'DapiQuery', the_base_url: str) -> None:
        """
        init and empty/nan values.
        Members described at class level
        """
        # URL
        self.api_url_base = the_base_url
        self.api_url_query = f'{the_base_url}/query'
        self.api_url_index = f'{the_base_url}/dsindex'
        self.api_url_dsblob = f'{the_base_url}/dsblob'
        self.api_url_blobdata = f'{the_base_url}/blobdata'
        # AVAILABLE
        self.available_jobtype = []
        self.available_platforms = []
        self.available_categories = []
        self.available_projects = []
        self.available_details = []
        self.available_dataversions = []
        self.available_testversions = []
        self.available_program = []
        self.available_data = []
        self.available_sources = []
        # SELECTED
        self.selected_jobtype = []
        self.selected_platforms = []
        self.selected_categories = []
        self.selected_projects = []
        self.selected_details = []
        self.selected_dataversions = []
        self.selected_testversions = []
        self.selected_program = []
        self.selected_data = []
        self.selected_sources = []
        # DATASETS
        self.available_datasets = []

    # pylint: disable=too-many-arguments,too-many-locals
    def update_from_selected(self: 'DapiQuery') -> None:
        print(f"update_from_selected self.api_url_query={self.api_url_query}", flush=True)
        search_args: Dict[str, List[str]] = {
            "mSources": self.selected_sources,
            "mProgram": self.selected_program,
            "mProjects": self.selected_projects,
            "mCategories": self.selected_categories,
            "mPlatforms": self.selected_platforms,
            "mData": self.selected_data,
            "mDetails": self.selected_details,
            "mDataVersions": self.selected_dataversions,
            "mTestVersions": self.selected_testversions,
            "mJobType": self.selected_jobtype
        }
        json_args: str = json.dumps(search_args)
        params: Dict[str, str] = {'search': json_args}
        print(f"update_from_selected params={params}", flush=True)
        response: requests = requests.get(self.api_url_query, params=params, timeout=60, allow_redirects=True)
        print(f"update_from_selected response={response}", flush=True)
        if response.ok:
            my_str: str = json.dumps(response.json(), indent=2)
            my_dict: Dict = json.loads(my_str)
            # headers (title), data
            self.available_datasets.clear()
            my_data_list: List[List[str]] = my_dict['data']
            my_dataset: List[str]
            for my_dataset in my_data_list:
                self.available_datasets.append(DapiEntry(my_dataset))
            # availableJobType, availablePlatforms, availableCategories, availableProjects, availableDetails,
            # availableDataVersions, availableTestVersions, availableProgram, availableData, availableSources
            self.available_jobtype.clear()
            self.available_jobtype.extend(my_dict['availableJobType'])
            self.available_platforms.clear()
            self.available_platforms.extend(my_dict['availablePlatforms'])
            self.available_categories.clear()
            self.available_categories.extend(my_dict['availableCategories'])
            self.available_projects.clear()
            self.available_projects.extend(my_dict['availableProjects'])
            self.available_details.clear()
            self.available_details.extend(my_dict['availableDetails'])
            self.available_dataversions.clear()
            self.available_dataversions.extend(my_dict['availableDataVersions'])
            self.available_testversions.clear()
            self.available_testversions.extend(my_dict['availableTestVersions'])
            self.available_program.clear()
            self.available_program.extend(my_dict['availableProgram'])
            self.available_data.clear()
            self.available_data.extend(my_dict['availableData'])
            self.available_sources.clear()
            self.available_sources.extend(my_dict['availableSources'])
        else:
            raise ValueError('Server returned bad response')
    # pylint: enable=too-many-arguments,too-many-locals

    def print_status(self: 'DapiQuery', the_entries: bool) -> None:
        print(f"self.api_url_base={self.api_url_base}", flush=True)
        print(f"self.available_jobtype={len(self.available_jobtype)}", flush=True)
        print(f"self.available_platforms={len(self.available_platforms)}", flush=True)
        print(f"self.available_categories={len(self.available_categories)}", flush=True)
        print(f"self.available_projects={len(self.available_projects)}", flush=True)
        print(f"self.available_details={len(self.available_details)}", flush=True)
        print(f"self.available_dataversions={len(self.available_dataversions)}", flush=True)
        print(f"self.available_testversions={len(self.available_testversions)}", flush=True)
        print(f"self.available_program={len(self.available_program)}", flush=True)
        print(f"self.available_data={len(self.available_data)}", flush=True)
        print(f"self.available_sources={len(self.available_sources)}", flush=True)
        print(f"self.selected_jobtype={len(self.selected_jobtype)}", flush=True)
        print(f"self.selected_platforms={len(self.selected_platforms)}", flush=True)
        print(f"self.selected_categories={len(self.selected_categories)}", flush=True)
        print(f"self.selected_projects={len(self.selected_projects)}", flush=True)
        print(f"self.selected_details={len(self.selected_details)}", flush=True)
        print(f"self.selected_dataversions={len(self.selected_dataversions)}", flush=True)
        print(f"self.selected_testversions={len(self.selected_testversions)}", flush=True)
        print(f"self.selected_program={len(self.selected_program)}", flush=True)
        print(f"self.selected_data={len(self.selected_data)}", flush=True)
        print(f"self.selected_sources={len(self.selected_sources)}", flush=True)
        print(f"self.available_datasets={len(self.available_datasets)}", flush=True)
        my_entry: DapiEntry
        if the_entries:
            for my_entry in self.available_datasets:
                my_entry.print_status(True)

    def find_data_versions(self: 'DapiQuery', the_dropdowns: List[Dict]) -> List[str]:
        versions: List[str] = []
        dropdown: Dict
        for dropdown in the_dropdowns:
            label: str = dropdown['entry_label']
            if label.startswith("DATA_"):
                versions.append(label)
            next_level: List[Dict] = dropdown['dropdown_entries']
            if len(next_level)>0:
                more_versions: List[str] = self.find_data_versions(next_level)
                versions.extend(more_versions)
                versions = list(set(versions))
                versions.sort()
        versions = list(set(versions))
        versions.sort()
        return versions

    def find_ngchm_files(self: 'DapiQuery', the_dropdowns: List[Dict]) -> List[str]:
        ngchms: List[str] = []
        dropdown: Dict
        for dropdown in the_dropdowns:
            label: str = dropdown['ngchm']
            if '' != label:
                ngchms.append(label)
            next_level: List[Dict] = dropdown['dropdown_entries']
            if len(next_level)>0:
                more_ngchms: List[str] = self.find_ngchm_files(next_level)
                ngchms.extend(more_ngchms)
                ngchms = list(set(ngchms))
                ngchms.sort()
        ngchms = list(set(ngchms))
        ngchms.sort()
        return ngchms

    def get_downloadable(self: 'DapiQuery', the_entry: DapiEntry) -> Tuple[List[str], List[str]]:
        data_versions: List[str] = []
        ngchm_files: List[str] = []
        print(f"get_index self.api_url_index={self.api_url_index}", flush=True)
        params: Dict[str, str] = {"id": the_entry.dataset_id}
        print(f"get_index params={params}", flush=True)
        response: requests = requests.get(self.api_url_index, params=params, timeout=60, allow_redirects=True)
        print(f"get_index response={response}", flush=True)
        if response.ok:
            # my_str: str = json.dumps(response.json(), indent=2)
            # print(f"get_index json={my_str}", flush=True)
            my_dict: Dict = response.json()
            dropdowns: List[Dict] = my_dict['mbatch']['dropdown_entries']
            data_versions = self.find_data_versions(dropdowns)
            print(f"get_index data_versions={data_versions}", flush=True)
            ngchm_files = self.find_ngchm_files(dropdowns)
            print("get_index ngchm_files=", flush=True, sep="\n")
            filename: str
            for filename in ngchm_files:
                print(f"  {filename}", flush=True, sep="\n")
        else:
            raise ValueError('Server returned bad response')
        return data_versions, ngchm_files

    def download_ngchm_html_to_file(self: 'DapiQuery', the_file_path: str, the_id: str, the_ngchm: str) -> None:
        print(f"download_ngchm_html_to_file self.api_url_dsblob={self.api_url_dsblob}", flush=True)
        if not the_ngchm.endswith(".html"):
            the_ngchm = f'{the_ngchm}.html'
        params: Dict[str, str] = {
            "id": the_id,
            "text": the_ngchm
        }
        print(f"download_ngchm_html_to_file params={params}", flush=True)
        response: requests = requests.get(self.api_url_dsblob, params=params, timeout=60, allow_redirects=True)
        print(f"download_ngchm_html_to_file response={response}", flush=True)
        if response.ok:
            file: BinaryIO
            chunk: bytes
            print(f"download_ngchm_html_to_file download to {the_file_path}", flush=True)
            with open(the_file_path, 'wb') as file:
                for chunk in response.iter_content(chunk_size=8192):
                    file.write(chunk)
        else:
            raise ValueError('Server returned bad response')

    def download_ngchm_ngchm_to_file(self: 'DapiQuery', the_file_path: str, the_id: str, the_ngchm: str) -> None:
        print(f"download_ngchm_ngchm_to_file self.api_url_dsblob={self.api_url_dsblob}", flush=True)
        if not the_ngchm.endswith(".ngchm"):
            if the_ngchm.endswith(".html"):
                the_ngchm = the_ngchm.replace(".html", "")
        params: Dict[str, str] = {
            "id": the_id,
            "text": the_ngchm
        }
        print(f"download_ngchm_ngchm_to_file params={params}", flush=True)
        response: requests = requests.get(self.api_url_dsblob, params=params, timeout=60, allow_redirects=True)
        print(f"download_ngchm_ngchm_to_file response={response}", flush=True)
        if response.ok:
            file: BinaryIO
            chunk: bytes
            print(f"download_ngchm_ngchm_to_file download to {the_file_path}", flush=True)
            with open(the_file_path, 'wb') as file:
                for chunk in response.iter_content(chunk_size=8192):
                    file.write(chunk)
        else:
            raise ValueError('Server returned bad response')

    def download_data_matrix_to_file(self: 'DapiQuery', the_file_path: str, the_id: str, the_version: str, the_original: bool) -> None:
        print(f"download_data_to_file self.api_url_blobdata={self.api_url_blobdata}", flush=True)
        # build path to data
        int_path: str = f"/versions/{the_version}/{'original' if the_original else 'pipeline'}/matrix.tsv"
        params: Dict[str, str] = {
            "id": the_id,
            "text": int_path
        }
        print(f"download_data_to_file params={params}", flush=True)
        response: requests = requests.get(self.api_url_blobdata, params=params, timeout=60, allow_redirects=True)
        print(f"download_data_to_file response={response}", flush=True)
        if response.ok:
            file: BinaryIO
            chunk: bytes
            print(f"download_data_to_file download to {the_file_path}", flush=True)
            with open(the_file_path, 'wb') as file:
                for chunk in response.iter_content(chunk_size=8192):
                    file.write(chunk)
        else:
            raise ValueError('Server returned bad response')

    def download_data_batches_to_file(self: 'DapiQuery', the_file_path: str, the_id: str, the_version: str, the_original: bool) -> None:
        print(f"download_data_batches_to_file self.api_url_blobdata={self.api_url_blobdata}", flush=True)
        # build path to data
        int_path: str = f"/versions/{the_version}/{'original' if the_original else 'pipeline'}/batches.tsv"
        params: Dict[str, str] = {
            "id": the_id,
            "text": int_path
        }
        print(f"download_data_batches_to_file params={params}", flush=True)
        response: requests = requests.get(self.api_url_blobdata, params=params, timeout=60, allow_redirects=True)
        print(f"download_data_batches_to_file response={response}", flush=True)
        if response.ok:
            file: BinaryIO
            chunk: bytes
            print(f"download_data_batches_to_file download to {the_file_path}", flush=True)
            with open(the_file_path, 'wb') as file:
                for chunk in response.iter_content(chunk_size=8192):
                    file.write(chunk)
        else:
            raise ValueError('Server returned bad response')
# pylint: too-many-arguments,enable=too-many-instance-attributes,too-few-public-methods
