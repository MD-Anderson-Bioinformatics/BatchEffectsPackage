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

from typing import List, Union, Dict
import itertools
import math
import matplotlib.transforms as mtrans
import matplotlib.pyplot as mplot
import matplotlib.lines as mlines
from PIL import Image

map_j2p_shapes: Dict[str, str] = \
    {
        '0': 's',
        '1': 'o',
        '2': '^',
        '3': '+',
        '4': 'x',
        '5': 'D',
        '6': 'v',
        '7': 'xs',
        '8': '+x',
        '9': '+D',
        '10': '+o',
        '11': '^v',
        '12': '+s',
        '13': 'xo',
        '14': 'vs'
    }


def convert_java_to_python_shape(the_shape: str) -> str:
    """
    Map Java/R shape code to Python shape code

    Shape					Java Code	Python Code
    square					0			s
    circle					1			o
    pt up triangle			2			^
    plus-cross				3			+
    x-cross					4			x
    diamond					5			D
    pt dn triangle			6			v
    square & x-cross		7			xs    was xs
    plus-cross & x-cross	8			+x
    plus-cross & square		9			+D
    plus-cross & circle		10			+o
    pt up & pt dn triangle	11			^v
    plus-cross & square		12			+s
    x-cross & circle		13			xo
    pt dn triangle & square	14			vs

    Note: Dict throws exception (as desired) if no mapping found

    :param the_shape: Java shape string
    :return: Python shape string
    """
    # print(f"convert_java_to_python_shape the_shape={the_shape}", flush=True)
    return_shape: str = map_j2p_shapes[the_shape]
    return return_shape


def export_legend(the_colors: List[str], the_shapes: List[str], the_labels: List[str], the_title: str, the_file_path: str) -> str:
    """
    Build a stand-along legend for the MBatch package, with color, shape, and labels

    original concept from
    https://stackoverflow.com/questions/4534480/get-legend-as-a-separate-picture-in-matplotlib

    :param the_colors: list of colors acceptable to Line2D
    :param the_shapes: list of marker shapes acceptable to Line2D
    :param the_labels: list of string labels acceptable to Line2D
    :param the_title: title for legend
    :param the_file_path: full path and name for PNG file
    :return: the_file_path, the full path and name for PNG file
    """
    # print(f"export_legend the_colors={the_colors}", flush=True)
    # print(f"export_legend type(the_colors)={type(the_colors)}", flush=True)
    # print(f"export_legend the_shapes={the_shapes}", flush=True)
    # print(f"export_legend type(the_shapes)={type(the_shapes)}", flush=True)
    # print(f"export_legend the_labels={the_labels}", flush=True)
    # print(f"export_legend type(the_labels)={type(the_labels)}", flush=True)
    # print(f"export_legend the_title={the_title}", flush=True)
    # print(f"export_legend type(the_title)={type(the_title)}", flush=True)
    print(f"export_legend the_file_path={the_file_path}", flush=True)
    print(f"export_legend type(the_file_path)={type(the_file_path)}", flush=True)
    handles: List[Union[tuple[mlines.Line2D, mlines.Line2D], mlines.Line2D]] = []
    color: str
    shape: str
    for (color, shape) in itertools.zip_longest(the_colors, the_shapes):
        # print(f"export_legend color={color}", flush=True)
        # print(f"export_legend shape={shape}", flush=True)
        # print(f"export_legend type(shape)={type(shape)}", flush=True)
        # is either mlines.Line2D or a tuple of Line2D
        if (color is None) & (shape is None):
            # print(f"0 export_legend color={color} shape={shape} type(shape)={type(shape)}", flush=True)
            my_line: mlines.Line2D = mlines.Line2D([], [], color='#000000', linestyle='none', fillstyle='none', marker='', markersize=10)
            handles.append(my_line)
        elif shape is None:
            # print(f"2 export_legend color={color} shape={shape} type(shape)={type(shape)}", flush=True)
            my_line: mlines.Line2D = mlines.Line2D([], [], color=color, linestyle='none', fillstyle='full', marker='s', markersize=10)
            handles.append(my_line)
        elif (color is None) & (2 == len(shape)):
            # print(f"1A export_legend color={color} shape={shape} type(shape)={type(shape)}", flush=True)
            my_line: tuple[mlines.Line2D, mlines.Line2D] = \
                (mlines.Line2D([], [], color='#000000', fillstyle='none', linestyle='none', marker=shape[0], markersize=10),
                 mlines.Line2D([], [], color='#000000', fillstyle='none', linestyle='none', marker=shape[1], markersize=10))
            handles.append(my_line)
        elif color is None:
            # print(f"1B export_legend color={color} shape={shape} type(shape)={type(shape)}", flush=True)
            my_line: mlines.Line2D = mlines.Line2D([], [], color='#000000', linestyle='none', fillstyle='none', marker=shape, markersize=10)
            handles.append(my_line)
        elif 2 == len(shape):
            # print(f"3 export_legend color={color} shape={shape} type(shape)={type(shape)}", flush=True)
            my_line: tuple[mlines.Line2D, mlines.Line2D] = \
                (mlines.Line2D([], [], color=color, fillstyle='none', linestyle='none', marker=shape[0], markersize=10),
                 mlines.Line2D([], [], color=color, fillstyle='none', linestyle='none', marker=shape[1], markersize=10))
            handles.append(my_line)
        else:
            # print(f"4 export_legend color={color} shape={shape} type(shape)={type(shape)}", flush=True)
            my_line: mlines.Line2D = mlines.Line2D([], [], color=color, linestyle='none', fillstyle='none', marker=shape, markersize=10)
            handles.append(my_line)
    # create a figure for the legend (allows us to turn off axes).
    # otherwise, axes appear to be border of legend, because legend takes up entire plot area.
    # print("export_legend plot figure", flush=True)
    fig_legend: mplot = mplot.figure(frameon=False, edgecolor='none', layout='tight')
    # print("export_legend plot legend", flush=True)
    legend: mplot = fig_legend.legend(handles, the_labels, loc=3, framealpha=0.0, frameon=False, edgecolor='none', title=the_title)
    # print("export_legend draw", flush=True)
    legend.figure.canvas.draw()
    # print("export_legend box", flush=True)
    bbox: mtrans.BboxBase = legend.get_window_extent().transformed(legend.figure.dpi_scale_trans.inverted())
    print(f"export_legend initial height '{bbox.height}'", flush=True)
    print(f"export_legend initial width '{bbox.width}'", flush=True)
    if bbox.height > 20.0:
        bbox = mtrans.Bbox.from_bounds(bbox.xmin, bbox.ymin, bbox.width, 20.0)
    print(f"export_legend final height '{bbox.height}'", flush=True)
    print(f"export_legend final width '{bbox.width}'", flush=True)
    print(f"export_legend save '{the_file_path}'", flush=True)
    legend.figure.savefig(the_file_path, dpi=300.0, bbox_inches=bbox)
    print("export_legend done", flush=True)
    return the_file_path


def export_legend_convert_shapes(the_colors: List[str], the_shapes: List[str], the_labels: List[str], the_title: str, the_file_path: str) -> str:
    """
    Build a stand-along legend for the MBatch package, with color, shape, and labels
    Convert shapes from Java/R to Python shapes

    :param the_colors: list of colors acceptable to Line2D
    :param the_shapes: list of marker shapes as used in older Java/R code
    :param the_labels: list of string labels acceptable to Line2D
    :param the_title: title for legend
    :param the_file_path: full path and name for PNG file
    :return: the_file_path, the full path and name for PNG file
    """
    # print(f"export_legend_convert_shapes the_colors={the_colors}", flush=True)
    # print(f"export_legend_convert_shapes type(the_colors)={type(the_colors)}", flush=True)
    # print(f"export_legend_convert_shapes the_shapes={the_shapes}", flush=True)
    # print(f"export_legend_convert_shapes type(the_shapes)={type(the_shapes)}", flush=True)
    # print(f"export_legend_convert_shapes the_labels={the_labels}", flush=True)
    # print(f"export_legend_convert_shapes type(the_labels)={type(the_labels)}", flush=True)
    # print(f"export_legend_convert_shapes the_title={the_title}", flush=True)
    # print(f"export_legend_convert_shapes type(the_title)={type(the_title)}", flush=True)
    print(f"export_legend_convert_shapes the_file_path={the_file_path}", flush=True)
    # print(f"export_legend_convert_shapes type(the_file_path)={type(the_file_path)}", flush=True)
    python_shapes: List[str] = []
    if len(the_shapes) > 0:
        my_shape: str
        for my_shape in the_shapes:
            # print(f"export_legend_convert_shapes my_shape={my_shape}", flush=True)
            python_shapes.append(convert_java_to_python_shape(my_shape))
    # print("export_legend_convert_shapes call export_legend", flush=True)
    return export_legend(the_colors, python_shapes, the_labels, the_title, the_file_path)


def combine_legends(the_list_of_files: List[str], the_filename_path: str) -> str:
    """
    Combined the images in the list into a new image written to the passed in filename
    :param the_list_of_files: list of full paths to file to combine
    :param the_filename_path: full path to file to write
    :return: the_filename_path
    """
    image_list: List[Image] = []
    image_file: str
    # print(f"combine_legends the_list_of_files={the_list_of_files}", flush=True)
    print(f"combine_legends the_filename_path={the_filename_path}", flush=True)
    for image_file in the_list_of_files:
        # print(f"combine_legends image_file={image_file}", flush=True)
        open_image: Image = Image.open(image_file)
        image_list.append(open_image)
    # determine largest legend sizes
    max_width: int = 0
    max_height: int = 0
    my_image: Image
    for my_image in image_list:
        if my_image.width > max_width:
            max_width = my_image.width
        if my_image.height > max_height:
            max_height = my_image.height
    # determine layout
    num_rows: int
    num_cols: int
    if len(the_list_of_files) <= 3:
        num_rows = 1
        num_cols = len(the_list_of_files)
    else:
        num_rows = 2
        num_cols = math.ceil(len(the_list_of_files) / 2)
    # calculate size
    image_height: int = num_rows * max_height
    # add padding based on columns
    image_width: int = (num_cols * max_width) + (num_cols * 5)
    row_index: int = 0
    col_index: int = 0
    # print(f"combine_legends max_width={max_width}", flush=True)
    # print(f"combine_legends max_height={max_height}", flush=True)
    # print(f"combine_legends num_rows={num_rows}", flush=True)
    # print(f"combine_legends num_cols={num_cols}", flush=True)
    # print(f"combine_legends image_width={image_width}", flush=True)
    # print(f"combine_legends image_height={image_height}", flush=True)
    #  2-tuple: (width, height)
    new_image: Image = Image.new("RGB", (image_width, image_height), "white")
    for my_image in image_list:
        # print(f"combine_legends loc indexes={col_index} {row_index}", flush=True)
        # 2-tuple: (column-location-width, row-location-height) for upper left corner of paste
        new_image.paste(my_image, ((col_index * max_width), (row_index * max_height)))
        row_index += 1
        if row_index >= num_rows:
            row_index = 0
            col_index += 1
    new_image.save(the_filename_path)
    return the_filename_path
