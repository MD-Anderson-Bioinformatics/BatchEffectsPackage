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


from typing import Dict
from pathlib import Path
import json


# pylint: disable=too-many-instance-attributes,too-many-arguments,too-few-public-methods
class OriginalData:
    """
    source:       value from original_data.json, used to populate dropdown menus for dataset
    project:      value from original_data.json, used to populate dropdown menus for dataset
    category:     value from original_data.json, used to populate dropdown menus for dataset
    platform:     value from original_data.json, used to populate dropdown menus for dataset
    data:         value from original_data.json, used to populate dropdown menus for dataset
    details:      value from original_data.json, used to populate dropdown menus for dataset
    """
    # do not set method variables, as they should be initialized in the init function
    source: str
    program: str
    project: str
    category: str
    platform: str
    data: str
    details: str
    version: str

    def __init__(self: 'OriginalData', the_source: str, the_program: str, the_project: str,
                 the_category: str, the_platform: str, the_data: str,
                 the_details: str, the_version: str) -> None:
        """
        init and empty/nan values.
        Members described at class level
        """
        super().__init__()
        self.source: str = the_source
        self.program: str = the_program
        self.project: str = the_project
        self.category: str = the_category
        self.platform: str = the_platform
        self.data: str = the_data
        self.details: str = the_details
        self.version: str = the_version
# pylint: enable=too-many-instance-attributes,too-many-arguments,too-few-public-methods


# Note, unlike Java GSON solutions, this requires updating when member attributes change
# uses object_decoder, since JSON attribute names allow values invalid for Python
def object_decoder_from_convert(the_obj: Dict[str, str]) -> OriginalData:
    """
    covert the string-based dict into an object for easier type-checking/hinting
    Uses values in converted data repo, no version
    :param the_obj: dictionary string/string object representing file data
    :return: instance of OriginalData
    """
    # manually add source and version
    return OriginalData("source-tbd",
                        the_obj['program'], the_obj['project'], the_obj['category'], the_obj['platform'],
                        the_obj['data'], the_obj['details'], "version-tbd")


# Note, unlike Java GSON solutions, this requires updating when member attributes change
# uses object_decoder, since JSON attribute names allow values invalid for Python
def object_decoder(the_obj: Dict[str, str]) -> OriginalData:
    """
    covert the string-based dict into an object for easier type-checking/hinting
    :param the_obj: dictionary string/string object representing file data
    :return: instance of OriginalData
    """
    return OriginalData(the_obj['source'], the_obj['program'], the_obj['project'],
                        the_obj['category'], the_obj['platform'], the_obj['data'],
                        the_obj['details'], the_obj['version'])


def read_json_original_data(the_file: str) -> OriginalData:
    """
    Read the file and convert from JSON string to OriginalData instance
    :param the_file: full path to the original json data file
    :return: OriginalData instance
    """
    print("read_json_original_data start", flush=True)
    my_obj: OriginalData
    if Path(the_file).exists():
        with open(the_file, 'r', encoding='utf-8') as my_file:
            print(f"read_json_original_data the_file={the_file}", flush=True)
            my_obj = json.load(my_file, object_hook=object_decoder)
            print(f"read_json_original_data my_obj={my_obj}", flush=True)
    else:
        my_obj = OriginalData("user-defined", "user-defined", "user-defined", "user-defined",
                              "user-defined", "user-defined", "user-defined", "user-defined")
    print(f"read_json_original_data finished {my_obj}", flush=True)
    return my_obj
