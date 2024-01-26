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

import os
from typing import List, Optional
import io
import matplotlib.pylab
import matplotlib.pyplot
from mbatch.test.common import convert_to_filename
from mbatch.volcano.calc import VolcanoData


# pylint: disable=too-many-arguments
def volcano_data(the_file_path: str, the_batch_type_list: List[str], the_data_list: List[VolcanoData],
                 the_log_frame_flag: bool, the_sub_dir_a: Optional[str] = None, the_sub_dir_b: Optional[str] = None) -> None:
    print(f"volcano_data the_file_path={the_file_path}", flush=True)
    batch_type: str
    for batch_type in the_batch_type_list:
        written: bool = False
        sub_dir: str = os.path.join(the_file_path, batch_type)
        if the_sub_dir_a is not None:
            sub_dir = os.path.join(sub_dir, the_sub_dir_a)
        if the_sub_dir_b is not None:
            sub_dir = os.path.join(sub_dir, the_sub_dir_b)
        if not os.path.exists(sub_dir):
            os.makedirs(sub_dir)
        print(f"volcano_data batch type dir={sub_dir}", flush=True)
        data_found: bool = False
        my_entry: VolcanoData
        for my_entry in the_data_list:
            if batch_type == my_entry.m_batch_type:
                data_found = True
                my_batch: str = my_entry.m_batch_a
                batch_data_file: str = os.path.join(sub_dir, f'Volcano-Data-{convert_to_filename(my_batch)}.json')
                print(f"volcano_data batch_data_file={batch_data_file}", flush=True)
                out_file: io.TextIOWrapper
                with open(batch_data_file, 'w', encoding='utf-8') as out_file:
                    my_entry.write_json(out_file)
                if not written:
                    calc_file: str = os.path.join(sub_dir, 'Volcano-Calc.txt')
                    with open(calc_file, 'w', encoding='utf-8') as out_file:
                        if the_log_frame_flag:
                            out_file.write("log-frame-data")
                        else:
                            out_file.write("linear-data")
                written = True
        if not data_found:
            # write error.log Unable to calculate volcano plot results\n
            error_file: str = os.path.join(sub_dir, 'error.log')
            print(f"volcano_data error_file={error_file}", flush=True)
            out_file: io.TextIOWrapper
            with open(error_file, 'w', encoding='utf-8') as out_file:
                out_file.write('Unable to calculate volcano plot results\n')
# pylint: enable=too-many-arguments


# pylint: disable=too-many-arguments
def volcano_plot(the_out_dir: str, the_data_list: List[VolcanoData], the_title: str, the_log_frame_flag: bool,
                 the_sub_dir_a: Optional[str] = None, the_sub_dir_b: Optional[str] = None) -> None:
    print(f"volcano_plot the_out_dir={the_out_dir}", flush=True)
    my_entry: VolcanoData
    for my_entry in the_data_list:
        my_batch_type: str = my_entry.m_batch_type
        base_path: str = os.path.join(the_out_dir, my_batch_type)
        if the_sub_dir_a is not None:
            base_path = os.path.join(base_path, the_sub_dir_a)
        if the_sub_dir_b is not None:
            base_path = os.path.join(base_path, the_sub_dir_b)
        if not os.path.exists(base_path):
            os.makedirs(base_path)
        png_path: str = os.path.join(base_path, f"Volcano-Diagram-{convert_to_filename(my_entry.m_batch_a)}.png")
        title_path: str = os.path.join(base_path, f"Volcano-Title-{convert_to_filename(my_entry.m_batch_a)}.txt")
        print(f"volcano_plot png_path={png_path}", flush=True)
        print(f"volcano_plot title_path={title_path}", flush=True)
        volcano_plot_obj(png_path, title_path, my_entry, the_title, the_log_frame_flag)
# pylint: enable=too-many-arguments


def volcano_plot_obj(the_png_path: str, the_title_path: str, the_data: VolcanoData, the_title: str, the_log_frame_flag: bool) -> None:
    print(f"volcano_plot_obj the_png_path={the_png_path}", flush=True)
    print(f"volcano_plot_obj the_title_path={the_title_path}", flush=True)
    my_title: str = f"{the_title} - {the_data.m_batch_type} - {the_data.m_batch_a}"
    print(f"volcano_plot_obj my_title={my_title}", flush=True)
    out_file: io.TextIOWrapper
    with open(the_title_path, 'w', encoding='utf-8') as out_file:
        out_file.write(my_title)
    # x = logged fold change
    # y = logged adjusted p-values
    figure: matplotlib.pylab.Figure
    axes: matplotlib.pylab.Axes
    figure, axes = matplotlib.pyplot.subplots()
    axes.scatter(x=the_data.m_fold_change, y=the_data.m_trans_pvalues, s=1)
    print(f"volcano_plot the_data.m_fold_change={min(the_data.m_fold_change)} <-> {max(the_data.m_fold_change)}", flush=True)
    # pylint: disable=nested-min-max
    axes.set_xlim(min(-3.0, min(the_data.m_fold_change)), max(3.0, max(the_data.m_fold_change)))
    print(f"volcano_plot the_data.m_trans_pvalues={min(the_data.m_trans_pvalues)} <-> {max(the_data.m_trans_pvalues)}", flush=True)
    axes.set_ylim(min(-3.0, min(the_data.m_trans_pvalues)), max(3.0, max(the_data.m_trans_pvalues)))
    # pylint: enable=nested-min-max
    if the_log_frame_flag:
        axes.set_title(f"{my_title} (log-frame-data)")
    else:
        axes.set_title(f"{my_title} (linear-data)")
    axes.set_xlabel("Log2 Fold Change")
    axes.set_ylabel("-Log10 Adjusted P-Value")
    axes.axvline(-2, color="grey", linestyle="--")
    axes.axvline(2, color="grey", linestyle="--")
    axes.axhline(2, color="grey", linestyle="--")
    figure.savefig(the_png_path)
    matplotlib.pyplot.close(figure)
    print("volcano_plot done", flush=True)
