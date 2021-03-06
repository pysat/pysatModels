#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2019, AGB & pysat team
# Full license can be found in License.md
#-----------------------------------------------------------------------------
"""
Routines to align and work with pairs of modelled and observational data

Routines
--------
compare_model_and_inst

"""

from __future__ import absolute_import
from __future__ import unicode_literals

import numpy as np

import verify  # PyForecastTools
import pysat

import pysatModelUtils as pysat_mu


def compare_model_and_inst(pairs=None, inst_name=[], mod_name=[],
                           methods=['all']):
    """Compare modelled and measured data

    Parameters
    ------------
    pairs : xarray.Dataset instance
        Dataset containing only the desired observation-model data pairs
    inst_name : list of strings
        ordered list of instrument measurements to compare to modelled data
    mod_name : list of strings
        ordered list of modelled data to compare to instrument measurements
    methods : list of strings
        statistics to calculate.  See Notes for accecpted inputs

    Returns
    ----------
    stat_dict : dict of dicts
        Dictionary where the first layer of keys denotes the instrument data
        name and the second layer provides the desired statistics
    data_units : dict
        Dictionary containing the units for the data

    Notes
    -----
    Statistics are calculated using PyForecastTools (imported as verify).
    See notes there for more details.

    all - all statistics
    all_bias - bias, meanPercentageError, medianLogAccuracy,
               symmetricSignedBias
    accuracy - returns dict with mean squared error, root mean squared error,
               mean absolute error, and median absolute error
    scaledAccuracy - returns dict with normaled root mean squared error, mean
                     absolute scaled error, mean absolute percentage error,
                     median absolute percentage error, median symmetric
                     accuracy
    bias - scale-dependent bias as measured by the mean error
    meanPercentageError - mean percentage error
    medianLogAccuracy - median of the log accuracy ratio
    symmetricSignedBias - Symmetric signed bias, as a percentage
    meanSquaredError - mean squared error
    RMSE - root mean squared error
    meanAbsError - mean absolute error
    medAbsError - median absolute error
    nRMSE - normaized root mean squared error
    scaledError - scaled error (see PyForecastTools for references)
    MASE - mean absolute scaled error
    forecastError - forecast error (see PyForecastTools for references)
    percError - percentage error
    absPercError - absolute percentage error
    logAccuracy - log accuracy ratio
    medSymAccuracy - Scaled measure of accuracy
    meanAPE - mean absolute percentage error

    """
    method_rout = {"bias": verify.bias, "accuracy": verify.accuracy,
                   "meanPercentageError": verify.meanPercentageError,
                   "medianLogAccuracy": verify.medianLogAccuracy,
                   "symmetricSignedBias": verify.symmetricSignedBias,
                   "meanSquaredError": verify.meanSquaredError,
                   "RMSE": verify.RMSE, "meanAbsError": verify.meanAbsError,
                   "medAbsError": verify.medAbsError, "MASE": verify.MASE,
                   "scaledAccuracy": verify.scaledAccuracy,
                   "nRMSE": verify.nRMSE, "scaledError": verify.scaledError,
                   "forecastError": verify.forecastError,
                   "percError": verify.percError, "meanAPE": verify.meanAPE,
                   "absPercError": verify.absPercError,
                   "logAccuracy": verify.logAccuracy,
                   "medSymAccuracy": verify.medSymAccuracy}

    replace_keys = {'MSE': 'meanSquaredError', 'MAE': 'meanAbsError',
                    'MdAE': 'medAbsError', 'MAPE': 'meanAPE',
                    'MdSymAcc': 'medSymAccuracy'}

    # Grouped methods for things that don't have convenience functions
    grouped_methods = {"all_bias": ["bias", "meanPercentageError",
                                    "medianLogAccuracy",
                                    "symmetricSignedBias"],
                       "all": list(method_rout.keys())}

    # Replace any group method keys with the grouped methods
    for gg in [(i, mm) for i, mm in enumerate(methods)
               if mm in list(grouped_methods.keys())]:
        # Extend the methods list to include all the grouped methods
        methods.extend(grouped_methods[gg[1]])
        # Remove the grouped method key
        methods.pop(gg[0])

    # Ensure there are no duplicate methods
    methods = list(set(methods))

    # Test the input
    if pairs is None:
        raise ValueError('must provide Dataset of paired observations')

    if len(inst_name) != len(mod_name):
        raise ValueError('must provide equal number of instrument and model ' +
                         'data names for comparison')

    if not np.all([iname in pairs.data_vars.keys() for iname in inst_name]):
        raise ValueError('unknown instrument data value supplied')

    if not np.all([iname in pairs.data_vars.keys() for iname in mod_name]):
        raise ValueError('unknown model data value supplied')

    if not np.all([mm in list(method_rout.keys()) for mm in methods]):
        known_methods = list(method_rout.keys())
        known_methods.extend(list(grouped_methods.keys()))
        unknown_methods = [mm for mm in methods
                           if mm not in list(method_rout.keys())]
        raise ValueError('unknown statistical method(s) requested:\n' +
                         '{:}\nuse only:\n{:}'.format(unknown_methods,
                                                      known_methods))

    # Initialize the output
    stat_dict = {iname: dict() for iname in inst_name}
    data_units = {iname: pairs.data_vars[iname].units for iname in inst_name}

    # Cycle through all of the data types
    for i, iname in enumerate(inst_name):
        # Determine whether the model data needs to be scaled
        iscale = pysat.utils.scale_units(pairs.data_vars[iname].units,
                                         pairs.data_vars[mod_name[i]].units)
        mod_scaled = pairs.data_vars[mod_name[i]].values.flatten() * iscale

        # Flatten both data sets, since accuracy routines require 1D arrays
        inst_dat = pairs.data_vars[iname].values.flatten()

        # Ensure no NaN are used in statistics
        inum = np.where(np.isfinite(mod_scaled) & np.isfinite(inst_dat))[0]

        # Calculate all of the desired statistics
        for mm in methods:
            try:
                stat_dict[iname][mm] = method_rout[mm](mod_scaled[inum],
                                                       inst_dat[inum])

                # Convenience functions add layers to the output, remove
                # these layers
                if hasattr(stat_dict[iname][mm], "keys"):
                    for nn in stat_dict[iname][mm].keys():
                        new = replace_keys[nn] if nn in replace_keys.keys() \
                            else nn
                        stat_dict[iname][new] = stat_dict[iname][mm][nn]
                    del stat_dict[iname][mm]
            except ValueError or NotImplementedError as err:
                # Not all data types can use all statistics.  Inform the user
                # instead of stopping processing.  Only valid statistics will
                # be included in output
                pysat_mu.logger.info("{:s} can't use {:s}: {:}".format(iname,
                                                                       mm, err))

    return stat_dict, data_units
