#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2022, pysat development team
# Full license can be found in License.md
# -----------------------------------------------------------------------------
"""Routines to align and work with pairs of modelled and observational data."""

import numpy as np

import pysat
import verify  # PyForecastTools

import pysatModels as ps_mod


def compare_model_and_inst(pairs, inst_name, mod_name, methods=['all'],
                           unit_label='units'):
    """Compare modelled and measured data.

    Parameters
    ----------
    pairs : xarray.Dataset
        Dataset containing only the desired observation-model data pairs
    inst_name : list
        Ordered list of strings indicating whicch instrument measurements to
         compare to modelled data
    mod_name : list
        Ordered list of strings indicating which modelled data to compare to
        instrument measurements
    methods : list
        Statistics to calculate. See Notes for accecpted inputs.
        (default=['all'])
    unit_label : str
        Unit attribute for data in `pairs` (default='units')

    Returns
    -------
    stat_dict : dict
        Dict of dicts where the first layer of keys denotes the instrument data
        name and the second layer provides the desired statistics
    data_units : dict
        Dict containing the units for the data

    Raises
    ------
    ValueError
        If input parameters are improperly formatted

    See Also
    --------
    PyForecastTools

    Notes
    -----
    Statistics are calculated using PyForecastTools (imported as verify).

    1. all: all statistics
    2. all_bias: bias, meanPercentageError, medianLogAccuracy,
       symmetricSignedBias
    3. accuracy: returns dict with mean squared error, root mean squared error,
       mean absolute error, and median absolute error
    4. scaledAccuracy: returns dict with normaled root mean squared error, mean
       absolute scaled error, mean absolute percentage error,  median absolute
       percentage error, median symmetric accuracy
    5. bias: scale-dependent bias as measured by the mean error
    6. meanPercentageError: mean percentage error
    7. medianLogAccuracy: median of the log accuracy ratio
    8. symmetricSignedBias: Symmetric signed bias, as a percentage
    9. meanSquaredError: mean squared error
    10. RMSE: root mean squared error
    11. meanAbsError: mean absolute error
    12. medAbsError: median absolute error
    13. nRMSE: normaized root mean squared error
    14. scaledError: scaled error (see PyForecastTools for references)
    15. MASE: mean absolute scaled error
    16. forecastError: forecast error (see PyForecastTools for references)
    17. percError: percentage error
    18. absPercError: absolute percentage error
    19. logAccuracy: log accuracy ratio
    20. medSymAccuracy: Scaled measure of accuracy
    21. meanAPE: mean absolute percentage error
    22. medAPE: median absolute perceentage error

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
                    'MdAE': 'medAbsError', 'MAPE': 'meanAPE', 'MdAPE': 'medAPE',
                    'MdSymAcc': 'medSymAccuracy'}

    # Grouped methods for things that don't have convenience functions
    grouped_methods = {"all_bias": ["bias", "meanPercentageError",
                                    "medianLogAccuracy",
                                    "symmetricSignedBias"],
                       "all": list(method_rout.keys())}

    # Ensure the list inputs are lists
    inst_name = pysat.utils.listify(inst_name)
    mod_name = pysat.utils.listify(mod_name)
    methods = pysat.utils.listify(methods)

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
        raise ValueError(''.join(['must provide equal number of instrument ',
                                  'and model data names for comparison']))

    if not np.all([iname in pairs.data_vars.keys() for iname in inst_name]):
        if not np.any([iname in pairs.data_vars.keys() for iname in inst_name]):
            raise ValueError('unknown instrument data value supplied')

        # Remove the missing instrument/model data variable names
        bad_ind = list()
        for i, iname in enumerate(inst_name):
            if iname not in pairs.data_vars.keys():
                bad_ind.append(i)

        ps_mod.logger.warning(
            'removed {:d} of the provided instrument variable names'.format(
                len(inst_name) - len(bad_ind)))

        for i in bad_ind:
            inst_name.pop(i)
            mod_name.pop(i)

    if not np.all([iname in pairs.data_vars.keys() for iname in mod_name]):
        if not np.any([iname in pairs.data_vars.keys() for iname in mod_name]):
            raise ValueError('unknown model data value supplied')

        # Remove the missing instrument/model data variable names
        bad_ind = list()
        for i, iname in enumerate(mod_name):
            if iname not in pairs.data_vars.keys():
                bad_ind.append(i)

        ps_mod.logger.warning(
            'removed {:d} of the provided model variable names'.format(
                len(mod_name) - len(bad_ind)))

        for i in bad_ind:
            inst_name.pop(i)
            mod_name.pop(i)

    if not np.all([mm in list(method_rout.keys()) for mm in methods]):
        known_methods = list(method_rout.keys())
        known_methods.extend(list(grouped_methods.keys()))
        unknown_methods = [mm for mm in methods
                           if mm not in list(method_rout.keys())]
        raise ValueError(''.join(['unknown statistical method(s) requested: ',
                                  '{:}\nuse only: '.format(unknown_methods),
                                  '{:}'.format(known_methods)]))

    # Initialize the output
    stat_dict = {iname: dict() for iname in inst_name}
    data_units = {iname: pairs.data_vars[iname].attrs[unit_label]
                  for iname in inst_name}

    # Cycle through all of the data types
    for i, iname in enumerate(inst_name):
        # Determine whether the model data needs to be scaled
        iscale = pysat.utils.scale_units(
            pairs.data_vars[iname].attrs[unit_label],
            pairs.data_vars[mod_name[i]].attrs[unit_label])
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
            except (ValueError, NotImplementedError, ZeroDivisionError) as err:
                # Not all data types can use all statistics.  Inform the user
                # instead of stopping processing.  Only valid statistics will
                # be included in output.
                # New version of pyForecastTools doesn't trigger this, using
                # Inf or masked output instead
                ps_mod.logger.info("{:s} can't use {:s}: {:}".format(iname,
                                                                     mm, err))

    return stat_dict, data_units
