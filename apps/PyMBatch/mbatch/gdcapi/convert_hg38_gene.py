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
class Hg38Gene:
    """
    Hg38Gene is a class that holds the gene mapping information,
    including:
        unique: str = Pulled from column 'unique'. Unique id from gene symbol and Ensembl id for gene with pipe between.
        gene_symbol: str = Pulled from column 'gene-symbol'.
        ensembl_id: str = Pulled from column 'ensembl-id'.
        chromosome: str = Pulled from column 'chromosome'. Strings, X, Y, 1-21. Mitochondrial and other "chromosomes" are already removed.
        start_loc: int = Pulled from column 'start-loc'. Very large integer.
        end_loc: int = Pulled from column 'end-loc'. Very large integer.
        strand: str = Pulled from column 'strand'. String + or -
        gene_type: str = Pulled from column 'gene-type'.
        calling_center: str = Pulled from column 'calling-center'.
    """
    # declare but do not set member attributes
    unique: str
    gene_symbol: str
    ensembl_id: str
    chromosome: str
    start_loc: int
    end_loc: int
    strand: str
    gene_type: str
    calling_center: str

    def __init__(self: 'Hg38Gene', the_dict: Dict[str, str]) -> None:
        super().__init__()
        self.unique: str = the_dict['unique']
        self.gene_symbol: str = the_dict['gene-symbol']
        self.ensembl_id: str = the_dict['ensembl-id']
        self.chromosome: str = the_dict['chromosome']
        self.start_loc: int = int(the_dict['start-loc'])
        self.end_loc: int = int(the_dict['end-loc'])
        self.strand: str = the_dict['strand']
        self.gene_type: str = the_dict['gene-type']
        self.calling_center: str = the_dict['calling-center']
# pylint: enable=too-many-instance-attributes,too-few-public-methods


# read HG38_Genes file into a Dict with gene_id as the key.
def read_hg38_gene_file(the_file: str, the_key_field: str) -> Dict[str, Hg38Gene]:
    """
    Read the HG38 gene file and populate dictionary of genes.
    :param the_file: Full path to file to read.
    :param the_key_field: What column to use for dictionary key? (Usually 'unique'.)
    :return: Dictionary of Hg38Gene objects with given dict key.
    """
    in_file: io.TextIOWrapper
    my_entries: Dict[str, Hg38Gene] = {}
    with open(the_file, 'r', encoding='utf-8') as in_file:
        keys: List[str] = read_headers(in_file)
        line: str = in_file.readline().rstrip('\n')
        index: int = 1
        while '' != line:
            if 0 == index % 1000:
                print(f"read_hg38_gene_file datafile index={index}", flush=True)
            index += 1
            values: List[str] = line.split("\t")
            tsv_dict: Dict[str, str] = dict(zip(keys, values))
            entry: Hg38Gene = Hg38Gene(tsv_dict)
            my_entries[getattr(entry, the_key_field)] = entry
            line = in_file.readline().rstrip('\n')
    return my_entries
