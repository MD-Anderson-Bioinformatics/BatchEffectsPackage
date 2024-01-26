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


from typing import List
import unittest
import os
import shutil
import csv
import pandas
import numpy
import scanpy
from mbatch.volcano.api import volcano_calc_plot


dynamic_test_volcano_dir: str = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/volcano"
dynamic_test_volcano_onebatch_dir: str = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/volcano_onebatch"
dynamic_test_volcano_ov_dir: str = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/volcano_ov"
dynamic_test_volcano_tumor_dir_linear: str = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/volcano_tumor_linear"
dynamic_test_volcano_tumor_dir_logframe: str = "/BEA/BatchEffectsPackage_data/testing_dynamic/PyMBatch/volcano_tumor_logframe"


# pylint: disable=too-many-instance-attributes
class TestVolcano(unittest.TestCase):
    """
    Class for setting up volcano plot testing - clear/make directory for output
    """
    # do not set method variables, as they should be initialized in the init function
    # No local method variables

    def setUp(self: 'TestVolcano') -> None:
        """
        setup script to clear and re-populate test directory
        :return:
        """
        #############################
        # files for configout testing
        if self._testMethodName == 'test_volcano_calc':
            print(f"TestJob::setUp dynamic_test_dir={dynamic_test_volcano_dir}", flush=True)
            if os.path.exists(dynamic_test_volcano_dir):
                shutil.rmtree(dynamic_test_volcano_dir)
            os.makedirs(dynamic_test_volcano_dir)
        if self._testMethodName == 'test_volcano_onebatch':
            print(f"TestJob::setUp dynamic_test_dir={dynamic_test_volcano_onebatch_dir}", flush=True)
            if os.path.exists(dynamic_test_volcano_onebatch_dir):
                shutil.rmtree(dynamic_test_volcano_onebatch_dir)
            os.makedirs(dynamic_test_volcano_onebatch_dir)
        if self._testMethodName == 'test_volcano_ov_calc':
            print(f"TestJob::setUp dynamic_test_dir={dynamic_test_volcano_ov_dir}", flush=True)
            if os.path.exists(dynamic_test_volcano_ov_dir):
                shutil.rmtree(dynamic_test_volcano_ov_dir)
            os.makedirs(dynamic_test_volcano_ov_dir)
        if self._testMethodName == 'test_volcano_tumor_calc_linear':
            print(f"TestJob::setUp dynamic_test_dir_linear={dynamic_test_volcano_tumor_dir_linear}", flush=True)
            if os.path.exists(dynamic_test_volcano_tumor_dir_linear):
                shutil.rmtree(dynamic_test_volcano_tumor_dir_linear)
            os.makedirs(dynamic_test_volcano_tumor_dir_linear)
        if self._testMethodName == 'test_volcano_tumor_calc_logframe':
            print(f"TestJob::setUp dynamic_test_dir_logframe={dynamic_test_volcano_tumor_dir_logframe}", flush=True)
            if os.path.exists(dynamic_test_volcano_tumor_dir_logframe):
                shutil.rmtree(dynamic_test_volcano_tumor_dir_logframe)
            os.makedirs(dynamic_test_volcano_tumor_dir_logframe)
        #############################
        print("TestVolcano:setUp done", flush=True)

    # TODO: flag to convert NA to 0?

    def test_volcano_calc(self: 'TestVolcano') -> None:
        """
        test the volcano.calc function
        :return: nothing
        """
        print("TestVolcano:test_volcano_calc start", flush=True)
        sample_id_col: str = 'Sample'
        print("TestVolcano:test_volcano_calc read matrix", flush=True)
        # read Standardized Data matrix
        my_matrix: pandas.DataFrame = pandas.read_csv("/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/StdData/TCGA/matrix_data.tsv", sep='\t', quoting=csv.QUOTE_NONE, encoding='utf-8', index_col=0)
        # log normalize data with a +1
        my_matrix = numpy.log10(my_matrix + 1)
        # make AnnData object from Standardized Data Matrix
        # need to transpose matrix, obs is samples, var is features
        anndata: scanpy.AnnData = scanpy.AnnData(X=my_matrix.values.transpose(), obs=my_matrix.columns.to_frame(), var=my_matrix.index.to_frame())
        # calculate highly variable
        scanpy.pp.highly_variable_genes(anndata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=500)
        # list of features to keep
        keep_features: List[str] = anndata.var.highly_variable.index.values[anndata.var.highly_variable.values]
        my_matrix = my_matrix.loc[keep_features]
        print("TestVolcano:test_volcano_calc read batches", flush=True)
        my_batches: pandas.DataFrame = pandas.read_csv("/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/StdData/TCGA/batches.tsv",
                                                       sep='\t', quoting=csv.QUOTE_NONE, encoding='utf-8', index_col=False, dtype=str)
        # reduce number of Batch Types
        my_batches = my_batches[[sample_id_col, 'ShipDate', 'TSS', 'PlateId']]
        print("TestVolcano:test_volcano_ov_calc call volcano_calc_plot", flush=True)
        volcano_calc_plot("Example Title", sample_id_col, my_matrix, my_batches, ['ShipDate', 'TSS', 'PlateId'],
                          dynamic_test_volcano_dir, True)
        print("TestVolcano:test_volcano_calc done", flush=True)

    def test_volcano_onebatch(self: 'TestVolcano') -> None:
        """
        test the volcano.onebatch function
        :return: nothing
        """
        print("TestVolcano:test_volcano_onebatch start", flush=True)
        sample_id_col: str = 'Sample'
        print("TestVolcano:test_volcano_onebatch read matrix", flush=True)
        # read Standardized Data matrix
        my_matrix: pandas.DataFrame = pandas.read_csv("/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/StdData/ONE_BATCH/matrix.tsv", sep='\t', quoting=csv.QUOTE_NONE, encoding='utf-8', index_col=0)
        # log normalize data with a +1
        my_matrix = numpy.log10(my_matrix + 1)
        # make AnnData object from Standardized Data Matrix
        # need to transpose matrix, obs is samples, var is features
        anndata: scanpy.AnnData = scanpy.AnnData(X=my_matrix.values.transpose(), obs=my_matrix.columns.to_frame(), var=my_matrix.index.to_frame())
        # calculate highly variable
        scanpy.pp.highly_variable_genes(anndata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=500)
        # list of features to keep
        keep_features: List[str] = anndata.var.highly_variable.index.values[anndata.var.highly_variable.values]
        my_matrix = my_matrix.loc[keep_features]
        print("TestVolcano:test_volcano_onebatch read batches", flush=True)
        my_batches: pandas.DataFrame = pandas.read_csv("/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/StdData/ONE_BATCH/batches.tsv",
                                                       sep='\t', quoting=csv.QUOTE_NONE, encoding='utf-8', index_col=False, dtype=str)
        # reduce number of Batch Types
        my_batches = my_batches[[sample_id_col, 'example']]
        print("TestVolcano:test_volcano_onebatch call volcano_calc_plot", flush=True)
        volcano_calc_plot("Example Title", sample_id_col, my_matrix, my_batches, ['example'],
                          dynamic_test_volcano_onebatch_dir, True)
        print("TestVolcano:test_volcano_onebatch done", flush=True)

    def test_volcano_ov_calc(self: 'TestVolcano') -> None:
        """
        test the volcano.calc function
        :return: nothing
        """
        print("TestVolcano:test_volcano_ov_calc start", flush=True)
        sample_id_col: str = 'Sample'
        print("TestVolcano:test_volcano_ov_calc read matrix", flush=True)
        # read Standardized Data matrix
        my_matrix: pandas.DataFrame = pandas.read_csv("/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/StdData/OV/matrix.tsv", sep='\t', quoting=csv.QUOTE_NONE, encoding='utf-8', index_col=0)
        # DO NOT log normalize data with a +1
        # rows with all zeros already removed
        # do not use scanpy.pp.highly_variable_genes, it fails on this data
        print("TestVolcano:test_volcano_calc read batches", flush=True)
        my_batches: pandas.DataFrame = pandas.read_csv("/BEA/BatchEffectsPackage_data/testing_static/PyMBatch/StdData/OV/batches.tsv",
                                                       sep='\t', quoting=csv.QUOTE_NONE, encoding='utf-8', index_col=False, dtype=str)
        # reduce number of Batch Types
        my_batches = my_batches[[sample_id_col, 'batch_id', 'sample_type_name']]
        print("TestVolcano:test_volcano_ov_calc call volcano_calc_plot", flush=True)
        volcano_calc_plot("Example Title", sample_id_col, my_matrix, my_batches, ['batch_id', 'sample_type_name'],
                          dynamic_test_volcano_ov_dir, True)
        print("TestVolcano:test_volcano_ov_calc done", flush=True)

    def test_volcano_tumor_calc_linear(self: 'TestVolcano') -> None:
        """
        test the volcano.calc function
        :return: nothing
        """
        print("TestVolcano:test_volcano_tumor_calc_linear start", flush=True)
        sample_id_col: str = 'Sample'
        print("TestVolcano:test_volcano_tumor_calc_linear read matrix", flush=True)
        my_matrix: pandas.DataFrame = pandas.read_csv("/BEA/BatchEffectsPackage_data/testing_static/MATRIX_DATA/matrix_data-Tumor.tsv", sep='\t', quoting=csv.QUOTE_NONE, encoding='utf-8', index_col=0)
        print("TestVolcano:test_volcano_tumor_calc_linear read batches", flush=True)
        my_batches: pandas.DataFrame = pandas.read_csv("/BEA/BatchEffectsPackage_data/testing_static/MATRIX_DATA/batches-Tumor.tsv",
                                                       sep='\t', quoting=csv.QUOTE_NONE, encoding='utf-8', index_col=False, dtype=str)
        # reduce number of Batch Types
        my_batches = my_batches[[sample_id_col, 'TSS', 'BatchId']]
        print("TestVolcano:test_volcano_tumor_calc_linear call volcano_calc_plot", flush=True)
        # not clear whether this data is log frame or not. Use not here.
        volcano_calc_plot("Example Title", sample_id_col, my_matrix, my_batches, ['TSS', 'BatchId'],
                          dynamic_test_volcano_tumor_dir_linear, False)
        print("TestVolcano:test_volcano_tumor_calc_linear done", flush=True)

    def test_volcano_tumor_calc_logframe(self: 'TestVolcano') -> None:
        """
        test the volcano.calc function
        :return: nothing
        """
        print("TestVolcano:test_volcano_tumor_calc_logframe start", flush=True)
        sample_id_col: str = 'Sample'
        print("TestVolcano:test_volcano_tumor_calc_logframe read matrix", flush=True)
        my_matrix: pandas.DataFrame = pandas.read_csv("/BEA/BatchEffectsPackage_data/testing_static/MATRIX_DATA/matrix_data-Tumor.tsv", sep='\t', quoting=csv.QUOTE_NONE, encoding='utf-8', index_col=0)
        print("TestVolcano:test_volcano_tumor_calc_logframe read batches", flush=True)
        my_batches: pandas.DataFrame = pandas.read_csv("/BEA/BatchEffectsPackage_data/testing_static/MATRIX_DATA/batches-Tumor.tsv",
                                                       sep='\t', quoting=csv.QUOTE_NONE, encoding='utf-8', index_col=False, dtype=str)
        # reduce number of Batch Types
        my_batches = my_batches[[sample_id_col, 'TSS', 'BatchId']]
        print("TestVolcano:test_volcano_tumor_calc_logframe call volcano_calc_plot", flush=True)
        # not clear whether this data is log frame or not. Use not here.
        volcano_calc_plot("Example Title", sample_id_col, my_matrix, my_batches, ['TSS', 'BatchId'],
                          dynamic_test_volcano_tumor_dir_logframe, True)
        print("TestVolcano:test_volcano_tumor_calc_logframe done", flush=True)
# pylint: enable=too-many-instance-attributes


if __name__ == '__main__':
    unittest.main()
