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
import io
import json
import pandas
import numpy
import scipy


# pylint: disable=too-many-arguments,too-many-instance-attributes,too-few-public-methods
class VolcanoData:
    """
    Store information on Volcano Calc results
    """
    # do not set method variables, as they should be initialized in the init function
    # #####
    # URL
    m_batch_type: str
    m_batch_a: str
    m_features: List[str]
    m_fold_change: List[float]
    m_pvalues: List[float]
    m_trans_pvalues: List[float]

    def __init__(self: 'VolcanoData', the_batch_type: str,
                 the_batch_a: str, the_features: List[str],
                 the_fold_change: List[float],
                 the_pvalues: List[float], the_trans_pvalues: List[float]) -> None:
        """
        init and empty/nan values.
        Members described at class level
        """
        self.m_batch_type = the_batch_type
        self.m_batch_a = the_batch_a
        self.m_features = the_features
        self.m_fold_change = the_fold_change
        self.m_trans_pvalues = the_trans_pvalues
        self.m_pvalues = the_pvalues

    def write_json(self: 'VolcanoData', the_out: io.TextIOWrapper) -> bool:
        written: bool = False
        # ####
        # convert to DICT, with fields
        # batch_type: str
        # batch_name: str
        # values: List[DICT]
        # ####
        # the values DICT have fields
        # feature: str
        # fold_change: float
        # trans_pvalue: float
        # pvalue: float
        # ####
        my_values: List[dict] = []
        feature: str
        fold_change: float
        trans_pvalue: float
        pvalue: float
        for feature, fold_change, trans_pvalue, pvalue in zip(self.m_features, self.m_fold_change, self.m_trans_pvalues, self.m_pvalues):
            sub_dict: dict = {
                'feature': feature,
                'fold_change': round(float(fold_change),4),
                'trans_pvalue': round(float(trans_pvalue),4),
                'pvalue': round(float(pvalue),4)
            }
            my_values.append(sub_dict)
        batch_dict: dict = {
            'batch_type': self.m_batch_type,
            'batch_name': self.m_batch_a,
            'values': my_values
        }
        # convert dict to JSON string
        my_str: str = json.dumps(batch_dict, indent=2)
        the_out.write(my_str)
        written = True
        return written
# pylint: enable=too-many-arguments,too-many-instance-attributes,too-few-public-methods


# pylint: disable=too-many-locals
def volcano_calc(the_sample_col: str, the_data: pandas.DataFrame, the_batches: pandas.DataFrame, the_log_frame_flag: bool) -> List[VolcanoData]:
    print(f"volcano_calc {the_sample_col}", flush=True)
    print(f"volcano_calc size = {the_data.shape}", flush=True)
    # remove rows that are all zero
    the_data = the_data[~(the_data == 0.0).all(axis=1)]
    print(f"volcano_calc removed zeros = {the_data.shape}", flush=True)
    # remove rows that are all +/- inf
    # replace infinite values with NaN
    the_data.replace([numpy.inf, -numpy.inf], numpy.nan, inplace=True)
    # drop rows with all NaN values
    the_data.dropna(axis=0, how='all', inplace=True)
    print(f"volcano_calc removed +/- inf = {the_data.shape}", flush=True)
    results: List['VolcanoData'] = []
    column_names: List[str] = list(the_batches.columns.values)
    column_names.remove(the_sample_col)
    feature_list: List[str] = the_data.index.values.tolist()
    name_group: str
    for name_group in column_names:
        # get list of batch values
        batch_values: List[str] = list(set(the_batches[name_group]))
        values_length: int = len(batch_values)
        print(f"Batch Type {name_group} has {values_length} values {batch_values}", flush=True)
        if values_length > 1:
            # process batch types with more than one batch
            index_a: int = 0
            while index_a < values_length:
                # get the batch to compare
                batch_a: str = batch_values[index_a]
                # get list of samples for each batch
                samples_a: List[str] = (the_batches[the_batches[name_group] == batch_a][the_sample_col]).to_list()
                samples_b: List[str] = the_data.columns.difference(samples_a).tolist()
                # get batch dataframes
                batch_df_a: pandas.DataFrame = the_data[samples_a]
                batch_df_b: pandas.DataFrame = the_data[samples_b]
                print(f"Batch Type {name_group} - {batch_a}", flush=True)
                # calculate fold changes
                means_a: numpy.ndarray = batch_df_a.to_numpy().mean(axis=1)
                means_b: numpy.ndarray = batch_df_b.to_numpy().mean(axis=1)
                # reviewed by rehan!
                if not the_log_frame_flag:
                    # out and where makes sure that NaNs and infinities become zero
                    means_a = numpy.log2(means_a, out=numpy.zeros_like(means_a), where=means_a > 0)
                    means_b = numpy.log2(means_b, out=numpy.zeros_like(means_b), where=means_b > 0)
                # since everything is log2 can now subtract
                foldchanges: List[float] = list(means_b - means_a)
                # calculate p-values
                # pvals: List[float] = []
                num_values: int = means_a.shape[0]
                row_int: int
                pvalue_list: List[float] = []
                for row_int in range(0, num_values):
                    # print(f"Batch Type {name_group}: {batch_a} vs {batch_b} for {batch_df_a.index[row_int]}", flush=True)
                    row_a: numpy.ndarray = batch_df_a.to_numpy()[row_int, :]
                    row_b: numpy.ndarray = batch_df_b.to_numpy()[row_int, :]
                    # print(f"min={row_a.min()}|{row_b.min()} max={row_a.max()}|{row_b.max()} ")
                    # TODO: any checks before T-Test?
                    ttest_result = scipy.stats.ttest_ind(row_a, row_b, nan_policy='omit')
                    pvalue: float = ttest_result.pvalue
                    if numpy.isnan(pvalue):
                        # use 1 as default to mean "ignore"
                        # prevents divide by zero error during log10
                        pvalue = 1
                    if numpy.ma.is_masked(pvalue):
                        pvalue = 1
                    pvalue_list.append(pvalue)
                transformed_pvals: List[float] = list(-1 * numpy.log10(num_values * numpy.array(pvalue_list)))
                results.append(VolcanoData(name_group, batch_a, feature_list, foldchanges, pvalue_list, transformed_pvals))
                # go to next batch a
                index_a += 1
                # print("--------------------------------------------------", flush=True)
    return results
# pylint: enable=too-many-locals
