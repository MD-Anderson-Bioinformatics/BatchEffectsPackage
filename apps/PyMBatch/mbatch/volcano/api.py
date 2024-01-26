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

from typing import List, Optional
import pandas
import numpy
from mbatch.volcano.calc import volcano_calc, VolcanoData
from mbatch.volcano.plot import volcano_plot, volcano_data


# pylint: disable=too-many-arguments
def volcano_calc_plot(the_title: str, the_sample_id_col: str, the_matrix: pandas.DataFrame, the_batches: pandas.DataFrame,
                      the_batch_types: List[str], the_output_dir: str, the_log_frame_flag: bool,
                      the_sub_dir_a: Optional[str] = None, the_sub_dir_b: Optional[str] = None) -> None:
    print("volcano_calc_plot perform calculations", flush=True)
    volcano_list: List[VolcanoData] = volcano_calc(the_sample_id_col, the_matrix, the_batches, the_log_frame_flag)
    print("volcano_calc_plot write data", flush=True)
    volcano_data(the_output_dir, the_batch_types, volcano_list, the_log_frame_flag, the_sub_dir_a, the_sub_dir_b)
    print("volcano_calc_plot plot data", flush=True)
    volcano_plot(the_output_dir, volcano_list, the_title, the_log_frame_flag, the_sub_dir_a, the_sub_dir_b)
    print("volcano_calc_plot done", flush=True)
# pylint: enable=too-many-arguments


# pylint: disable=too-many-arguments
def volcano_calc_plot_from_r(the_title: str, the_sample_id_col: str,
                             the_matrix: numpy.ndarray, the_features: List[str], the_samples: List[str],
                             the_batches: numpy.ndarray,
                             the_batch_types: List[str], the_output_dir: str, the_log_frame_flag: bool,
                             the_sub_dir_a: Optional[str] = None, the_sub_dir_b: Optional[str] = None) -> None:
    print("volcano_calc_plot_from_r make matrix", flush=True)
    my_matrix: pandas.DataFrame = pandas.DataFrame(the_matrix, index=the_features, columns=the_samples)
    print("volcano_calc_plot_from_r make batches", flush=True)
    print(the_batches)
    my_batches: pandas.DataFrame = pandas.DataFrame(data=the_batches, dtype=str)
    print("volcano_calc_plot_from_r call plot", flush=True)
    volcano_calc_plot(the_title, the_sample_id_col, my_matrix, my_batches, the_batch_types,
                      the_output_dir, the_log_frame_flag, the_sub_dir_a, the_sub_dir_b)
# pylint: enable=too-many-arguments
