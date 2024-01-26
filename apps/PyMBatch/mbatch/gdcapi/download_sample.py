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
from typing import Dict, List, Optional
from mbatch.test.common import write_tsv_list, read_headers


SAMPLE_HEADERS: List[str] = [
    'sample_uuid',
    'sample_barcode'
]


# pylint: disable=too-many-instance-attributes,too-many-arguments,too-few-public-methods
class GdcApiSample:
    """
    Class for mapping sample UUID to sample Barcode.
    """
    # declare but do not set member attributes
    sample_uuid: str
    sample_barcode: str

    def __init__(self: 'GdcApiSample', the_dict: Optional[Dict] = None,
                 the_barcode: Optional[str] = None, the_uuid: Optional[str] = None) -> None:
        """
        init and empty/nan values.
        Members described at class level
        """
        super().__init__()
        # declare instance variables
        self.sample_uuid: str = ""
        self.sample_barcode: str = ""
        # set  instance variables
        if the_dict is not None:
            for entry in SAMPLE_HEADERS:
                setattr(self, entry, the_dict.get(entry))
        else:
            self.sample_barcode = the_barcode
            self.sample_uuid = the_uuid

    def write_sample(self: 'GdcApiSample', the_out_file: io.TextIOWrapper) -> None:
        """
        Write sample row (line) to streaming wrapper using header array.
        :param the_out_file: TextIOWrapper writing stream
        :return: Nothing
        """
        my_samples: List[str] = []
        for attr in SAMPLE_HEADERS:
            my_samples.append(getattr(self, attr))
        write_tsv_list(the_out_file, my_samples, True, False)
# pylint: enable=too-many-instance-attributes,too-many-arguments,too-few-public-methods


def write_sample_index_file(the_index_file: str, the_samples: List[GdcApiSample]) -> None:
    """
    Sort and write sample index file.
    :param the_index_file: the full path and filename for sample index.
    :param the_samples: list of sample objects for this index file
    :return: Nothing
    """
    print(f"write_sample_index_file {len(the_samples)} to the_index_file={the_index_file}", flush=True)
    the_samples.sort(key=lambda my_sample: my_sample.sample_barcode, reverse=False)
    out_file: io.TextIOWrapper
    index: int
    with open(the_index_file, 'w', encoding='utf-8') as out_file:
        write_tsv_list(out_file, SAMPLE_HEADERS, True, False)
        my_entry: GdcApiSample
        for index, my_entry in enumerate(the_samples):
            if 0 == index % 10000:
                print(f"write_sample_index_file index={index}", flush=True)
            my_entry.write_sample(out_file)


def read_sample_index_file(the_index_file: str, the_sample_dict: Dict[str, GdcApiSample]) -> None:
    """
    Read sample index file and populate sample dictionary passed in.
    :param the_index_file: sample file index file to read
    :param the_sample_dict: sample dictionary to populate
    :return: Nothing
    """
    # print(f"read_sample_index_file the_index_file={the_index_file}", flush=True)
    in_file: io.TextIOWrapper
    with open(the_index_file, 'r', encoding='utf-8') as in_file:
        keys: List[str] = read_headers(in_file)
        line: str = in_file.readline().rstrip('\n')
        index: int = 1
        while '' != line:
            # if 0 == index % 10000:
            #     print(f"read_sample_index_file samples index={index}", flush=True)
            index += 1
            values: List[str] = line.split("\t")
            tsv_dict: Dict[str, str] = dict(zip(keys, values))
            entry: GdcApiSample = GdcApiSample(tsv_dict, None, None)
            the_sample_dict[entry.sample_uuid] = entry
            line = in_file.readline().rstrip('\n')
