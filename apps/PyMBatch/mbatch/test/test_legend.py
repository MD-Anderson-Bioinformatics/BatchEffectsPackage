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


import unittest
import os
import shutil
from typing import List
from PIL import Image, ImageChops
from mbatch.legend.legend import export_legend, export_legend_convert_shapes
from mbatch.legend.legend import convert_java_to_python_shape, combine_legends


def test_legend(the_shape_flag: bool, the_file_path: str, the_empty_symbols: bool, the_empty_colors: bool,
                the_huge_legend: bool) -> str:
    """
    create a test legend
    :param the_shape_flag: if true use python shapes, otherwise java
    :param the_file_path: full path to filename
    :param the_empty_symbols: if true use symbols, otherwise send empty list
    :param the_empty_colors: if true use colors, otherwise send empty list
    :param the_huge_legend: if true use make huge lengend for test
    :return: path to filename
    """
    print(f"test_legend the_file_path={the_file_path}", flush=True)
    colors: List[str] = ["#FF0000", "#00FF00", "#000000", "#FA0000",
                         "#AD0000", "#00FA00", "#00AF00", "#0000FA",
                         "#0000AF", "#000000", "#AFAFAF", "#FF00FF",
                         "#00FFFF", "#FFFF00", "#FAFFAF"]
    shapes: List[str] = ['s', 'o', '^', '+',
                         'x', 'D', 'v', 'xs',
                         '+x', '+D', '+o', '^v',
                         '+s', 'xo', 'vs']
    java_shapes: List[str] = ['0', '1', '2', '3',
                              '4', '5', '6', '7',
                              '8', '9', '10', '11',
                              '12', '13', '14']
    labels: List[str] = ["00 A1 - UCSF (13)",
                         "01 GM - MD Anderson (15)",
                         "02 D8 - Greater Poland Cancer Center (63)",
                         "03 HN - Ontario Institute for Cancer Research (OICR) (1)",
                         "04 Lorem ipsum dolor sit amet, consectetur adipiscing elit",
                         "05 Curabitur non eros eget nibh egestas mollis",
                         "06 Ut ipsum sem, rutrum eget aliquam luctus, posuere at metus",
                         "07 Nam nisi dui, pulvinar vitae dignissim ullamcorper, dignissim id metus",
                         "08 Lorem ipsum dolor sit amet, consectetur adipiscing elit",
                         "09 Mauris posuere arcu et leo fringilla rhoncus",
                         "10 Donec nec neque tellus, sed accumsan magna",
                         "11 In massa sem, pretium vitae aliquet vehicula, varius id risus",
                         "12 Aenean nec cursus enim",
                         "13 Suspendisse congue elementum facilisis",
                         "14 Nulla non faucibus libero"]
    if the_empty_symbols:
        shapes = []
        java_shapes = []
    if the_empty_colors:
        colors = []
    if the_huge_legend:
        for _ in range(100000):
            colors.append("#FAFFAF")
            shapes.append('vs')
            labels.append("14 Nulla non faucibus libero")
        export_legend(colors, shapes, labels, "Huge Legend Title", the_file_path)
    elif the_shape_flag:
        export_legend(colors, shapes, labels, "Test Legend Title", the_file_path)
    else:
        export_legend_convert_shapes(colors, java_shapes, labels, "Test Legend Title Convert", the_file_path)
    return the_file_path


def two_png_same_p(the_file_a: str, the_file_b: str) -> bool:
    """
    Use ImageChops to see if images are indentical
    :param the_file_a: full path to image file a
    :param the_file_b: full path to image file b
    :return: true if they are the same
    """
    im1: Image
    im2: Image
    with Image.open(the_file_a) as im1:
        with Image.open(the_file_b) as im2:
            diff: Image = ImageChops.difference(im2, im1)
            if None is diff.getbbox():
                return True
    return False


# pylint: disable=too-many-instance-attributes
class TestLegend(unittest.TestCase):
    """
    Class for setting up Legend testing - clear/make directory for output
    """
    # declare but do not set member attributes
    dyn_leg_java: str
    dyn_leg_java_no_symbol: str
    dyn_leg_java_no_color: str
    dyn_leg_java_no_color_symbol: str
    dyn_leg_pyth: str
    dyn_leg_pyth_no_symbol: str
    dyn_leg_pyth_no_color: str
    dyn_leg_pyth_no_color_symbol: str
    dyn_leg_comb: str
    dyn_leg_comb_single: str
    sta_leg_java: str
    sta_leg_pyth: str
    sta_leg_comb: str

    def setUp(self: 'TestLegend'):
        # called from init
        self.dyn_leg_java: str = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_java.png"
        self.dyn_leg_java_no_symbol: str = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_java_no_symbol.png"
        self.dyn_leg_java_no_color: str = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_java_no_color.png"
        self.dyn_leg_java_no_color_symbol: str = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_java_no_color_symbol.png"
        self.dyn_leg_pyth: str = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_python.png"
        self.dyn_leg_pyth_no_symbol: str = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_python_no_symbol.png"
        self.dyn_leg_pyth_no_color: str = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_python_no_color.png"
        self.dyn_leg_pyth_no_color_symbol: str = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_python_no_color_symbol.png"
        self.dyn_leg_comb: str = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_combined.png"
        self.dyn_leg_comb_single: str = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_combined_single.png"
        self.dyn_leg_huge: str = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/Legend/legend_huge.png"
        self.sta_leg_java: str = "/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/Legend/legend_java.png"
        self.sta_leg_pyth: str = "/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/Legend/legend_python.png"
        self.sta_leg_comb: str = "/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/Legend/legend_combined.png"
        if os.path.exists(os.path.dirname(self.dyn_leg_java)):
            shutil.rmtree(os.path.dirname(self.dyn_leg_java))
        os.makedirs(os.path.dirname(self.dyn_leg_java))

    def test_convert_java_to_python_shape(self: 'TestLegend') -> None:
        """
        test the convert_java_to_python_shape function
        :return: nothing
        """
        print("test_convert_java_to_python_shape", flush=True)
        self.assertEqual(convert_java_to_python_shape('0'), 's')
        self.assertEqual(convert_java_to_python_shape('1'), 'o')
        self.assertEqual(convert_java_to_python_shape('2'), '^')
        self.assertEqual(convert_java_to_python_shape('3'), '+')
        self.assertEqual(convert_java_to_python_shape('4'), 'x')
        self.assertEqual(convert_java_to_python_shape('5'), 'D')
        self.assertEqual(convert_java_to_python_shape('6'), 'v')
        self.assertEqual(convert_java_to_python_shape('7'), 'xs')
        self.assertEqual(convert_java_to_python_shape('8'), '+x')
        self.assertEqual(convert_java_to_python_shape('9'), '+D')
        self.assertEqual(convert_java_to_python_shape('10'), '+o')
        self.assertEqual(convert_java_to_python_shape('11'), '^v')
        self.assertEqual(convert_java_to_python_shape('12'), '+s')
        self.assertEqual(convert_java_to_python_shape('13'), 'xo')
        self.assertEqual(convert_java_to_python_shape('14'), 'vs')

    def test_export_legend(self: 'TestLegend') -> None:
        """
        test the export_legend function - uses Python shape codes
        :return: nothing
        """
        print("test_export_legend", flush=True)
        test_legend(True, self.dyn_leg_pyth, False, False, False)
        test_legend(True, self.dyn_leg_pyth_no_symbol, True, False, False)
        test_legend(True, self.dyn_leg_pyth_no_color, False, True, False)
        test_legend(True, self.dyn_leg_pyth_no_color_symbol, True, True, False)
        # compare self.dyn_leg_pyth to self.sta_leg_pyth
        self.assertTrue(two_png_same_p(self.dyn_leg_pyth, self.sta_leg_pyth))

    def test_export_huge_legend(self: 'TestLegend') -> None:
        """
        test the export_legend function - uses Python shape codes
        :return: nothing
        """
        print("test_export_huge_legend start", flush=True)
        test_legend(True, self.dyn_leg_huge, False, False, True)
        print("test_export_huge_legend done", flush=True)

    def test_export_legend_convert_shapes(self: 'TestLegend') -> None:
        """
        test the export_legend_convert_shapes function - uses Java shape codes
        :return: nothing
        """
        print("test_export_legend_convert_shapes", flush=True)
        test_legend(False, self.dyn_leg_java, False, False, False)
        test_legend(False, self.dyn_leg_java_no_symbol, True, False, False)
        test_legend(False, self.dyn_leg_java_no_color, False, True, False)
        test_legend(False, self.dyn_leg_java_no_color_symbol, True, True, False)
        # compare self.dyn_leg_java to self.sta_leg_java
        self.assertTrue(two_png_same_p(self.dyn_leg_java, self.sta_leg_java))

    def test_combine_legends_multiple(self: 'TestLegend') -> None:
        """
        test the combine_legends function
        :return: nothing
        """
        print("test_combine_legends_multiple", flush=True)
        combine_legends([self.sta_leg_java, self.sta_leg_pyth, self.sta_leg_java, self.sta_leg_pyth, self.sta_leg_java], self.dyn_leg_comb)
        # compare self.dyn_leg_java to self.sta_leg_java
        self.assertTrue(two_png_same_p(self.dyn_leg_comb, self.sta_leg_comb))

    def test_combine_legends_single(self: 'TestLegend') -> None:
        """
        test the combine_legends function
        :return: nothing
        """
        print("test_combine_legends_single", flush=True)
        combine_legends([self.sta_leg_java], self.dyn_leg_comb_single)
        # only tests for not failing with one legend
# pylint: enable=too-many-instance-attributes


if __name__ == '__main__':
    unittest.main()
