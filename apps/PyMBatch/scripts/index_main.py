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

import os
from mbatch.test.test_index import create_index_archive
from mbatch.test.common import delete_directory_contents

result_dir: str = "/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/index/ZIP-RESULTS"
info_dir: str = "/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/index/RESULTS/info"
data_dir: str = "/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/index/ZIP-DATA"
zip_dir: str = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/index/"

if __name__ == '__main__':
    delete_directory_contents(zip_dir)
    os.makedirs(zip_dir, exist_ok=True)
    create_index_archive(result_dir, data_dir, zip_dir, info_dir, None, None)
