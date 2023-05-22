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

from mbatch.test.test_legend import test_legend
from mbatch.legend.legend import combine_legends


dyn_leg_java: str = "/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_java.png"
dyn_leg_pyth: str = "/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_python.png"
dyn_leg_comb: str = "/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_combined.png"
sta_leg_java: str = "/BatchEffectsPackage_data/testing_static/PyMBatch/Legend/legend_java.png"
sta_leg_pyth: str = "/BatchEffectsPackage_data/testing_static/PyMBatch/Legend/legend_python.png"
sta_leg_comb: str = "/BatchEffectsPackage_data/testing_static/PyMBatch/Legend/legend_combined.png"


if __name__ == '__main__':
    test_legend(True, "/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_python.png", False, False, False)
    test_legend(True, "/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_python_no_symbol.png", True, False, False)
    test_legend(True, "/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_python_no_color.png", False, True, False)
    test_legend(True, "/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_python_no_color_symbol.png", True, True, False)
    test_legend(False, "/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_java.png", False, False, False)
    test_legend(False, "/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_java_no_symbol.png", True, False, False)
    test_legend(False, "/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_java_no_color.png", False, True, False)
    test_legend(False, "/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_java_no_color_symbol.png", True, True, False)
    combine_legends([sta_leg_java, sta_leg_pyth, sta_leg_java, sta_leg_pyth, sta_leg_java], dyn_leg_comb)
