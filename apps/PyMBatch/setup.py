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


from setuptools import setup, find_packages

setup(name='mbatch',
      version='1.0',
      # list folders, not files
      packages=find_packages(),
      # script files
      scripts=['scripts/legend_main.py', 'scripts/dsc_main.py'],
      # No package data (yet)
      # package_data={'capitalize': ['data/cap_data.txt']},
      install_requires=['matplotlib', 'pandas', 'numpy', 'pillow', 'jsonpickle', 'requests', 'xmltodict'],
      )
