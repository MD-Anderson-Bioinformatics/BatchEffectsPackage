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
from typing import Dict, List
from mbatch.test.common import read_headers


# pylint: disable=too-many-instance-attributes,too-few-public-methods
class Hg38Meth:
    """
    Hg38Meth represents the Meth450/27 (SeSAME) methylation probes.
    Probe information includes:
        self.probe_id: str = Pulled from column 'probe-id'. Unique probe id.
        chromosome: str = Pulled from column 'chromosome'. Strings, X, Y, 1-21. Mitochondrial and other "chromosomes" are already removed.
        start_loc: int = Pulled from column 'start'. Very large integer.
        end_loc: int = Pulled from column 'end'. Very large integer.
        strand: str = Pulled from column 'strand'. String + or -
        genes: List[str] = Pulled from column 'genes'. ; delimited list. Gene symbols.
    """
    # declare but do not set member attributes
    probe_id: str
    chromosome: str
    start_loc: int
    end_loc: int
    strand: str
    genes: List[str]

    def __init__(self: 'Hg38Meth', the_dict: Dict[str, str]) -> None:
        super().__init__()
        self.probe_id: str = the_dict['probe-id']
        self.chromosome: str = the_dict['chromosome']
        self.start_loc: int = int(the_dict['start'])
        self.end_loc: int = int(the_dict['stop'])
        self.strand: str = the_dict['strand']
        gene_string: str = the_dict['genes']
        self.genes: List[str] = gene_string.split(";")
# pylint: enable=too-many-instance-attributes,too-few-public-methods


# read HG38_Meth file into a Dict with probe as the key.
def read_hg38_meth_file(the_file: str) -> Dict[str, Hg38Meth]:
    """
    Read the HG38 meth file and populate dictionary of methylation probes.
    :param the_file: Full path to file to read.
    :return: Dictionary of Hg38Gene objects with probe_id as dict key.
    """
    in_file: io.TextIOWrapper
    my_entries: Dict[str, Hg38Meth] = {}
    with open(the_file, 'r', encoding='utf-8') as in_file:
        keys: List[str] = read_headers(in_file)
        line: str = in_file.readline().rstrip('\n')
        index: int = 1
        while '' != line:
            if 0 == index % 1000:
                print(f"read_hg38_meth_file datafile index={index}", flush=True)
            index += 1
            values: List[str] = line.split("\t")
            tsv_dict: Dict[str, str] = dict(zip(keys, values))
            entry: Hg38Meth = Hg38Meth(tsv_dict)
            my_entries[entry.probe_id] = entry
            line = in_file.readline().rstrip('\n')
    return my_entries
