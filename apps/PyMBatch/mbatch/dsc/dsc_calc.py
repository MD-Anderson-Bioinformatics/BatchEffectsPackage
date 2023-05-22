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

from typing import List
import math
import time
import numpy
from mbatch.dsc.dsc_info import DscInfo, epsilon_zero_check_value


# pylint: disable=too-many-locals,too-many-statements
def dsc_calc(the_matrix: numpy.ndarray, the_batches: numpy.ndarray, the_time_flag: bool = False) -> DscInfo:
    """
    Perform DSC calculations for the given data.
    For more details see: https://bioinformatics.mdanderson.org/public-software/tcga-batch-effects/#the-dsc-metric
    DSC is a ratio of between batch dispersion vs. within batch dispersion.
    Dw can roughly be viewed as the average “distance” between samples within a batch and the batch mean, or centroid.
    Db can roughly be viewed as the average distance between batch centroids and global mean.
    Same three values are also calculated and returned for each feature.
    :param the_matrix: samples across the top, features down the side, Decimal values
    :param the_batches: list of strings with batch ids for samples
    :param the_time_flag: true to write time string (defaults to False)
    :return: DscInfo object contain results of DSC calculation
    """
    # print("*******************start******************", flush=True)
    start: float = time.time()
    # check if the_batches is a dataframe -- the way R passes value
    # not needed for numpy?
    # if isinstance(the_batches, pandas.DataFrame):
    #    # print("Squeeze dataframe to series", flush=True)
    #    the_batches = the_batches.squeeze()
    # past_squeeze: float = time.time()
    # print(f"dsc_calc time to squeeze={(past_squeeze-start)} seconds", flush=True)
    # print(f"dsc_calc the_matrix={the_matrix}", flush=True)
    # print(f"dsc_calc the_batches={the_batches}", flush=True)
    # print(f"dsc_calc type(the_matrix)={type(the_matrix)}", flush=True)
    # print(f"dsc_calc type(the_batches)={type(the_batches)}", flush=True)
    # print(f"dsc_calc the_matrix.shape={the_matrix.shape}", flush=True)
    batches_cnt: int = the_batches.size
    # print(f"dsc_calc batches_cnt={batches_cnt}", flush=True)
    # print(f"dsc_calc type(batches_cnt)={type(batches_cnt)}", flush=True)
    # feature_cnt: int = the_matrix.shape[0]
    # print(f"dsc_calc feature_cnt={feature_cnt}", flush=True)
    # print(f"dsc_calc type(feature_cnt)={type(feature_cnt)}", flush=True)
    sample_cnt: int = the_matrix.shape[1]
    # print(f"dsc_calc sample_cnt={sample_cnt}", flush=True)
    # NOTE: assumes samples in StdData and batches are in the same order
    # print("Number of batches should match number of samples (columns)", flush=True)
    assert sample_cnt == batches_cnt, "Number of batches should match number of samples (columns)"
    # get unique list of batches
    # print("get unique list of batches", flush=True)
    unique_batches: numpy.ndarray = numpy.unique(the_batches)
    # print(f"len(unique_batches)={len(unique_batches)}", flush=True)
    # to calculate overall dw and db values
    dw_value: float = 0.0
    db_value: float = 0.0
    # to calculate feature dw and db values
    dsc_feature_list: List[float] = []
    dw_feature_list: List[float] = []
    db_feature_list: List[float] = []
    # iterate over features
    feature_counter: int = 0
    # print("iterate over features", flush=True)
    # use indexes to iterate over features
    feature_row: numpy.ndarray
    for feature_row in the_matrix:
        # print("====ROW=====", flush=True)
        feature_counter += 1
        # print(f"feature_name={feature_name}", flush=True)
        # print(f"feature_row.shape={feature_row.shape}", flush=True)
        # if 0 == math.remainder(feature_counter, 1000):
        #     print(f"feature_counter={feature_counter}", flush=True)
        dw_feature: float = 0.0
        db_feature: float = 0.0
        batch_name: str
        # print("get feature_row", flush=True)
        # feature_row: pandas.Series = the_matrix.loc[feature_name]
        # print("get feature_row_no_nan", flush=True)
        # remove_nan_pre: float = time.time()
        # do not do in separate function as python function calls are expensive
        feature_row_no_nan: numpy.ndarray = feature_row[~numpy.isnan(feature_row)]
        # remove_nan_post: float = time.time()
        # print(f"dsc_calc time to remove_nan={(remove_nan_post - remove_nan_pre)} seconds", flush=True)
        # print("get feature_mean", flush=True)
        # mean_pre: float = time.time()
        # use numpy math -- faster
        feature_mean: float = numpy.ndarray.mean(feature_row_no_nan)
        # print(f"feature_mean={feature_mean}", flush=True)
        # mean_post: float = time.time()
        # print(f"dsc_calc time to mean={(mean_post - mean_pre)} seconds", flush=True)
        # print("iterate over batch names", flush=True)
        # batch_variance_pre: float = time.time()
        for batch_name in unique_batches:
            # print("-----------------------------", flush=True)
            # print(f"batch_name={batch_name}", flush=True)
            # do without function call -- save expense
            type_count_tuple = numpy.unique(the_batches, return_counts=True)
            batch_sample_count: int = type_count_tuple[1][numpy.where(type_count_tuple[0] == batch_name)[0][0]]
            indices_for_batch: numpy.ndarray = numpy.where(the_batches == batch_name)
            batch_values: numpy.ndarray = feature_row[indices_for_batch]
            # batch_mean: float = statistics.mean(batch_values)
            batch_mean: float = numpy.ndarray.mean(batch_values)
            batch_variance: float = 0.0
            if batch_values.size > 1:
                # this gives us population variance
                batch_variance = numpy.ndarray.var(batch_values)
                # sample variance is pop-var * (batch_values.size/((batch_values.size - 1)*1.0))
                batch_variance = batch_variance * (batch_values.size / ((batch_values.size - 1) * 1.0))
                # batch_variance = statistics.variance(batch_values)
            # print(f"pre add dw_feature={dw_feature}", flush=True)
            # print(f"pre add db_feature={db_feature}", flush=True)
            # print(f"delta dw_feature={((float(batch_sample_count-1) / float(sample_cnt)) * batch_variance)}", flush=True)
            # print(f"delta db_feature={((float(batch_sample_count-1) / float(sample_cnt)) * batch_variance)}", flush=True)
            dw_feature = dw_feature + ((float(batch_sample_count-1) / float(sample_cnt)) * batch_variance)
            db_feature = db_feature + ((float(batch_sample_count) / float(sample_cnt)) *
                                       (batch_mean - feature_mean) * (batch_mean - feature_mean))
            # print(f"batch_sample_count={batch_sample_count}", flush=True)
            # print(f"batch_mean={batch_mean}", flush=True)
            # print(f"batch_variance={batch_variance}", flush=True)
            # print(f"dw_feature={dw_feature}", flush=True)
            # print(f"db_feature={db_feature}", flush=True)
        # print("------------done-------------", flush=True)
        # batch_variance_post: float = time.time()
        # print(f"dsc_calc time to batch_variance={(batch_variance_post - batch_variance_pre)} seconds", flush=True)
        # calculate feature values
        # print(f"pre sum dw_value={dw_value}", flush=True)
        # print(f"pre sum db_value={db_value}", flush=True)
        dw_value = dw_value + dw_feature
        db_value = db_value + db_feature
        # print(f"post sum dw_value={dw_value}", flush=True)
        # print(f"post sum db_value={db_value}", flush=True)
        # print(f"pre sqrt dw_feature={dw_feature}", flush=True)
        # print(f"pre sqrt db_feature={db_feature}", flush=True)
        dw_feature = math.sqrt(dw_feature)
        # print("between", flush=True)
        db_feature = math.sqrt(db_feature)
        # print(f"post sqrt dw_feature={dw_feature}", flush=True)
        # print(f"post sqrt db_feature={db_feature}", flush=True)
        # determine feature DSC
        # feature_dsc_pre: float = time.time()
        dsc_feature: float = float('nan')
        dw_feature = epsilon_zero_check_value(dw_feature)
        if 0.0 == dw_feature:
            dsc_feature = float('Inf')
        elif 0.0 == epsilon_zero_check_value(db_feature):
            dsc_feature = 0.0
        elif 0.0 != dw_feature:
            dsc_feature = db_feature / dw_feature
        # add to feature lists
        dsc_feature_list.append(dsc_feature)
        dw_feature_list.append(dw_feature)
        db_feature_list.append(db_feature)
        # feature_dsc_post: float = time.time()
        # print(f"dsc_calc time to feature_dsc={(feature_dsc_post - feature_dsc_pre)} seconds", flush=True)
    # calculate overall values
    # overall_pre: float = time.time()
    dw_value = math.sqrt(dw_value)
    db_value = math.sqrt(db_value)
    dsc_value = float('nan')
    if dw_value > 0:
        dsc_value = db_value / dw_value
    # handle too few samples case
    if sample_cnt < 2:
        # makes a list of 0.0 that is sample_cnt long
        dsc_feature_list = [0.0] * sample_cnt
        db_feature_list = [0.0] * sample_cnt
        dsc_value = 0.0
        db_value = 0.0
    # populate and return result
    result_info: DscInfo = DscInfo()
    result_info.m_dsc = dsc_value
    result_info.m_db = db_value
    result_info.m_dw = dw_value
    result_info.m_list_of_feature_dsc = dsc_feature_list
    result_info.m_list_of_feature_db = db_feature_list
    result_info.m_list_of_feature_dw = dw_feature_list
    result_info.epsilon_zero_check()
    # overall_post: float = time.time()
    # print(f"dsc_calc time to overall={(overall_post - overall_pre)} seconds", flush=True)
    finish: float = time.time()
    if the_time_flag:
        print(f"dsc_calc time to run={(finish-start)} seconds", flush=True)
    return result_info
# pylint: enable=too-many-locals,too-many-statements
