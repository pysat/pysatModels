#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2019, AGB & pysat team
# Full license can be found in License.md
# -----------------------------------------------------------------------------
"""Routines to extract observational-style data from model output."""

import numpy as np
import re
import scipy.interpolate as interpolate

import pandas as pds
import pysat.utils as pyutils

import pysatModels as ps_mod


def instrument_altitude_to_model_pressure(inst, model, inst_name, mod_name,
                                          mod_datetime_name, mod_time_name,
                                          mod_units, inst_alt, mod_alt,
                                          mod_alt_units, scale=100.,
                                          inst_out_alt='model_altitude',
                                          inst_out_pres='model_pressure',
                                          tol=1.0):
    """Interpolates altitude values onto model pressure levels.

    Parameters
    ----------
    inst : pysat.Instrument
        Instrument object with observational data
    model : xarray.Dataset
        Model data in xarray format
    inst_name : array-like
        List of variable names containing the observational data coordinates
        at which the model data will be interpolated. Must be in the same order
        as mod_name.
    mod_name : array-like
        list of names of the coordinates to use for determining model locations.
        Must be in the same order as `mod_alt` is stored within xarray.
        The coordinates must make up a regular grid.
    mod_datetime_name : str
        Name of the data series in the model Dataset containing datetime info
    mod_time_name : str
        Name of the time coordinate in the model Dataset
    mod_units : list
        units for each of the mod_name location attributes.  Currently
        supports: rad/radian(s), deg/degree(s), h/hr(s)/hour(s), m, km, and cm
    inst_alt : str
        String identifier used in inst for the altitude variable
    mod_alt : str
        Variable identifier for altitude data in the model
        e.g. 'ZG' in standard TIEGCM files.
    mod_alt_units : str
        units for the altitude variable. Currently supports: m, km, and cm
    scale : float
        Scalar used to roughly translate a change in altitude with a
        change in pressure level, the scale height. Same units as used by inst.
        (default=100.)
    inst_out_alt : str
        Label assigned to the model altitude data when attached to inst
        (default='model_altitude').
    inst_out_pres : str
        Label assigned to the model pressure level when attached to inst
        (default='model_pressure').
    tol : float
        Allowed difference between observed and modelled altitudes.
        Interpreted to have the same units as `inst_alt` (default=1.0).

    Returns
    -------
    [inst_out_alt, inst_out_pres] : list
        List of keys corresponding to the modelled data that was added to the
        instrument.

    Raises
    ------
    ValueError
        For incorrect input arguments

    Notes
    -----
    Uses an iterative regular grid interpolation to find the
    appropriate pressure level for the given input locations.

    """

    # Ensure the coordinate and data variable names are array-like
    inst_name = pyutils.listify(inst_name)
    mod_name = pyutils.listify(mod_name)

    # Test input
    if len(inst_name) == 0:
        estr = 'Must provide inst_name as a list of strings.'
        raise ValueError(estr)

    if len(mod_name) == 0:
        estr = 'Must provide mod_name as a list of strings.'
        raise ValueError(estr)

    if len(inst_name) != len(mod_name):
        estr = 'Must provide the same number of instrument and model '
        estr += 'location attribute names as a list'
        raise ValueError(estr)

    if len(mod_name) != len(mod_units):
        raise ValueError(''.join(['Must provide units for each model location',
                                  ' attribute']))

    if mod_time_name not in model.coords:
        raise ValueError("Unknown model time coordinate key name")

    if inst_alt not in inst.variables:
        raise ValueError("Unknown Instrument altitude key name")

    if mod_alt not in model.data_vars:
        raise ValueError("Unknown Model altitude key name")

    # Determine the scaling between model and instrument data
    inst_scale = np.ones(shape=len(inst_name), dtype=float)
    for i, iname in enumerate(inst_name):
        if iname not in inst.data.keys():
            raise ValueError(''.join(['Unknown instrument location index ',
                                      '{:}'.format(iname)]))
        # Altitude variable check
        if iname == inst_alt:
            alt_scale = pyutils.scale_units(
                mod_alt_units, inst.meta[iname, inst.meta.labels.units])
            inst_scale[i] = 1.
        else:
            inst_scale[i] = pyutils.scale_units(
                mod_units[i], inst.meta[iname, inst.meta.labels.units])

    # Pull out datetime data
    if mod_datetime_name in model.data_vars:
        mod_datetime = model.data_vars[mod_datetime_name]
        mod_datetime = mod_datetime.values.astype(np.int64)
    elif mod_datetime_name in model.coords:
        mod_datetime = model.coords[mod_datetime_name].values.astype(np.int64)
    else:
        raise ValueError("".join(["unknown model name for datetime: ",
                                  mod_datetime_name]))

    # Create initial fake regular grid index in inst
    inst_model_coord = inst[inst_name[0]] * 0

    # We need to create altitude index from the model.
    # First, handle model locations for interpolation.
    # We use the dimensions associated with model altitude
    # in the order provided.
    points = [model[dim].values / temp_scale
              for dim, temp_scale in zip(mod_name, inst_scale)]
    # Add time in first
    points.insert(0, mod_datetime)

    # Create interpolator
    interp = interpolate.RegularGridInterpolator(points,
                                                 np.log(model[mod_alt].values
                                                        / alt_scale),
                                                 bounds_error=False,
                                                 fill_value=None)
    # Use this interpolator to figure out what altitudes we are at.
    # Each loop, use the current estimated path through model expressed in
    # pressure and interp above to get the equivalent altitude of this path.
    # Compare this altitude to the actual instrument altitude and
    # shift the equivalent pressure for the instrument up/down
    # until difference between altitudes is small.

    # Log of instrument altitude
    log_ialt = np.log(inst[inst_alt])

    # Initial difference signal
    diff = log_ialt * 0 + 2.0 * tol
    while np.any(np.abs(diff) > tol):
        # Create input array using satellite time/position.
        # Replace the altitude coord with the fake tiegcm one
        coords = []
        for iscale, coord in zip(inst_scale, inst_name):
            if coord == inst_alt:
                # Don't scale altitude-like model coordinate
                coords.append(inst_model_coord)
            else:
                # Scale other dimensions to the model
                coords.append(inst[coord] * iscale)

        coords.insert(0, inst.index.values.astype(int))
        # To peform the interpolation we need points
        # like (x1, y1, z1), (x2, y2, z2)
        # but we currently have a list like
        # [(x1, x2, x3), (y1, y2, y3), ...].
        # To convert to the form we need the * below breaks each main list out,
        # and zip repeatedly outputs a tuple (xn, yn, zn).
        sat_pts = [inp for inp in zip(*coords)]

        # Altitude pulled out from model
        orbit_alt = interp(sat_pts)
        # Difference in altitude
        diff = np.e ** orbit_alt - np.e ** log_ialt
        # Shift index in inst for model pressure level in the opposite direction
        # to diff, reduced in value by scale, the 'scale height'
        inst_model_coord -= diff / scale

    # Store achieved model altitude
    inst[inst_out_alt] = np.e ** orbit_alt
    notes_str = ''.join(('Interpolated Model altitude corresponding to `',
                         inst_out_pres, '` pressure level.',
                         inst.meta[inst_alt, inst.meta.labels.notes]))
    inst.meta[inst_out_alt] = {inst.meta.labels.units:
                               inst.meta[inst_alt, inst.meta.labels.units],
                               inst.meta.labels.notes: notes_str}

    # Add pressure level that goes with altitude to Instrument
    inst[inst_out_pres] = inst_model_coord
    notes_str = ''.join(('Interpolated Model pressure corresponding to `',
                         inst_out_alt, '` altitude level.'))

    # TODO(#95) add units from Model when that input is upgraded to be
    # a pysat Custom compatible method
    inst.meta[inst_out_pres] = {inst.meta.labels.notes: notes_str}

    return [inst_out_alt, inst_out_pres]


def instrument_view_through_model(inst, model, inst_name, mod_name,
                                  mod_datetime_name, mod_time_name,
                                  mod_units, sel_name=None,
                                  methods=['linear'], model_label='model'):
    """Interpolates model values onto instrument locations.

    Parameters
    ----------
    inst : pysat.Instrument
        Instrument object with observational data
    model : xarray
        Modelled data
    inst_name : array-like
        List of variable names containing the observational data coordinates
        at which the model data will be interpolated. Do not include 'time',
        only spatial coordinates.
    mod_name : array-like
        List of model dimension names used for organizing model data
        in the same order as `inst_name`.  These must make up a regular grid.
        Do not include 'time', only spatial dimensions.
    mod_datetime_name : str
        Name of the data series in the model Dataset containing datetime info.
    mod_time_name : str
        Name of the time coordinate in the model Dataset.
    mod_units : list
        Units for each of the mod_name location dimensions.  Currently
        supports: rad/radian(s), deg/degree(s), h/hr(s)/hour(s), m, km, and cm
    sel_name : array-like or NoneType
        List of names of modelled data indices to append to Instrument object,
        or None to append all modelled data. (default=None)
    methods : str
        'linear' interpolation or 'nearest' neighbor options for RegularGrid.
        Must supply an option for each variable. (default=['linear'])
    model_label : str
        Name of model, used to identify interpolated data values in instrument
        (default="model")

    Returns
    -------
    interp_data.keys() : Keys
        Keys of modelled data added to the instrument

    Raises
    ------
    ValueError
        For incorrect input arguments

    Notes
    -----
    Updates the inst Instrument with interpolated data from the model
    Instrument. The interpolation is performed via the RegularGridInterpolator
    for quick performance.

    This method may require the use of a pre-processor on coordinate
    dimensions to ensure that a regular interpolation may actually be
    performed.

    Models, such as TIEGCM, have a regular grid in pressure, not in altitude.
    To use this routine for TIEGCM please use
    instrument_altitude_to_model_pressure first to transform instrument
    altitudes to pressure levels suitable for this method.

    Variables that vary exponentially in height
    may be approximated by taking a log before interpolating, though
    this does also yield an exponential variation along the horizontal
    directions as well.

    Expects units strings to have the units as the first word, if a long
    description is provided (e.g., 'degrees', 'degrees North', or 'deg_N' and
    not 'geographic North degrees')

    See Also
    --------
    pysat.utils.scale_units

    """

    # Ensure the coordinate and data variable names are array-like
    inst_name = pyutils.listify(inst_name)
    mod_name = pyutils.listify(mod_name)
    methods = pyutils.listify(methods)

    # Interpolate over all variables if None provided
    if sel_name is None:
        sel_name = pyutils.listify(model.data_vars.keys())
    else:
        sel_name = pyutils.listify(sel_name)

    if len(methods) != len(sel_name):
        estr = ' '.join(('Must provide interpolation selection',
                        'for each variable via methods keyword.'))
        raise ValueError(estr)

    for method in methods:
        if method not in ['linear', 'nearest']:
            estr = ''.join(('Methods only supports "linear" or "nearest".',
                            ' Not ', method))
            raise ValueError(estr)

    # Test input
    if len(inst_name) == 0:
        estr = 'Must provide inst_name as a list of strings.'
        raise ValueError(estr)

    if len(mod_name) == 0:
        estr = 'Must provide mod_name as a list of strings.'
        raise ValueError(estr)

    if len(sel_name) == 0:
        estr = 'Must provide sel_name as a list of strings.'
        raise ValueError(estr)

    if len(inst_name) != len(mod_name):
        estr = 'Must provide the same number of instrument and model '
        estr += 'location attribute names as a list'
        raise ValueError(estr)

    if len(mod_name) != len(mod_units):
        raise ValueError(''.join(['Must provide units for each model location',
                                  ' attribute']))

    if mod_time_name not in model.coords:
        raise ValueError("Unknown model time coordinate key name")

    for name in sel_name:
        if name not in model:
            estr = ''.join((name, ' is not a valid model variable.'))
            raise ValueError(estr)

    if mod_datetime_name in model.data_vars:
        mod_datetime = model.data_vars[mod_datetime_name]
        mod_datetime = mod_datetime.values.astype(np.int64)
    elif mod_datetime_name in model.coords:
        mod_datetime = model.coords[mod_datetime_name].values.astype(np.int64)
    else:
        raise ValueError("".join(["unknown model name for datetime: ",
                                  mod_datetime_name]))

    # Determine the scaling between model and instrument data
    inst_scale = np.ones(shape=len(inst_name), dtype=np.float64)
    for i, iname in enumerate(inst_name):
        if iname not in inst.data.keys():
            raise ValueError(''.join(['Unknown instrument location index ',
                                      '{:}'.format(iname)]))

        # Some units may have extra information (e.g., 'degrees North').
        # Use only the actual units in the scaling function.  These are assumed
        # to come first.
        long_units = re.split(r"\W+|_",
                              inst.meta[iname, inst.meta.labels.units])[0]
        inst_scale[i] = pyutils.scale_units(mod_units[i], long_units)

    del_list = list()
    keep_list = list()
    for mdat in sel_name:
        mname = '_'.join((model_label, mdat))
        if mname in inst.data.keys():
            ps_mod.logger.warning("".join(["model data already interpolated:",
                                           " {:}".format(mname)]))
            del_list.append(mdat)
        else:
            keep_list.append(mdat)

    if len(sel_name) == len(del_list):
        raise ValueError("instrument object already contains all model data")
    elif len(del_list) > 0:
        sel_name = keep_list

    # Create inst input based upon provided dimension names
    coords = [inst[coord_name] for coord_name in inst_name]

    # Time goes first
    coords.insert(0, inst.index.values.astype(np.int64))

    # Move from a list of lists [ [x1, x2, ...], [y1, y2, ...]]
    # to a list of tuples
    # [(x1, y1, ...), (x2, y2, ...)]
    # required for regulargrid interpolator
    inst_pts = [inp for inp in zip(*coords)]

    # Perform the interpolation
    interp = {}
    output_names = []
    for label, method in zip(sel_name, methods):
        # Determine the unit scaling between model and instrument data
        inst_scale = np.ones(shape=len(inst_name), dtype=float)

        # Collect model grid points together
        points = []

        # Time dim first
        points.append(mod_datetime)

        # Now spatial
        for iscale, var in zip(inst_scale, mod_name):
            points.append(model[var].values / iscale)

        # Create the interpolator
        interp[label] = interpolate.RegularGridInterpolator(
            points, model[label].values, bounds_error=False, fill_value=None,
            method=method)

        # Apply it at observed locations and store result
        output_names.append('_'.join((model_label, label)))
        inst[output_names[-1]] = interp[label](inst_pts)

    return output_names


def interp_inst_w_irregular_model_coord(inst, model, inst_name, mod_name,
                                        mod_datetime_name, mod_units,
                                        mod_reg_dim, mod_irreg_var,
                                        mod_var_delta, sel_name=None,
                                        model_label='model'):
    """Interpolate model data onto Instrument path using an irregular coordinate.

    Parameters
    ----------
    inst : pysat.Instrument
        pysat object that will receive interpolated data based upon position.
    model : pysat.Instrument
        Xarray pysat Instrument with model data that will be interpolated onto
        the `inst` locations.
    inst_name : list
        List of variable names containing the instrument data coordinates
        at which the model data will be interpolated. Do not include 'time',
        only spatial coordinates. Same ordering as used by mod_name.
    mod_name : list
        List of names of the data dimensions used to organize model data,
        in the same order as `inst_name`. These dimensions must make up a
        regular grid. Values from `mod_irreg_var` will be used to replace
        one of these regular dimensions, `mod_reg_dim`, with irregular values.
    mod_datetime_name : str
        Name of the data series in the model Dataset containing datetime info.
    mod_units : list
        Units for each of the `mod_name` dimensions. Users must provide units
        for `mod_irreg_var' instead of the units for `mod_reg_dim`.
        Currently supports: rad/radian(s), deg/degree(s), h/hr(s)/hour(s), m,
        km, and cm.
    mod_reg_dim : str
        Existing regular dimension name (must be in `mod_name`) used to
        organize model data that will be replaced with values from
        `mod_irreg_var` before performing interpolation.
    mod_irreg_var : str
        Variable name in model used to define the irregular grid value
        locations along `mod_reg_dim`. Must have same dimensions as `mod_name`.
    mod_var_delta : list
        List of delta values to be used when downselecting model values
        before interpolation, max(min(inst) - delta, min(model)) <= val <=
        min(max(inst) + delta, max(model)). Interpreted in the same order as
        `mod_name`.
    sel_name : list
        List of strings denoting model variable names that will be
        interpolated onto inst. The coordinate dimensions for these variables
        must correspond to those in `mod_irreg_var`.
    model_label : str
        name of model, used to identify interpolated data values in instrument
        (default="model")

    Returns
    -------
    output_names : list
        Keys of interpolated model data added to the instrument

    Raises
    ------
    ValueError
        For incorrect input arguments

    Notes
    -----
    Expects units strings to have the units as the first word, if a long
    description is provided (e.g., 'degrees', 'degrees North', or 'deg_N' and
    not 'geographic North degrees').

    See Also
    --------
    pysat.utils.scale_units

    """

    # Ensure the inputs are array-like
    inst_name = pyutils.listify(inst_name)

    # interp over all vars if None provided
    if sel_name is None:
        sel_name = list(model.data_vars.keys())
    else:
        sel_name = pyutils.listify(sel_name)

    # Test input
    if len(inst_name) == 0:
        estr = 'Must provide inst_name as a list of strings.'
        raise ValueError(estr)

    if len(mod_name) == 0:
        estr = 'Must provide mod_name as a list of strings.'
        raise ValueError(estr)

    if len(sel_name) == 0:
        estr = 'Must provide sel_name as a list of strings.'
        raise ValueError(estr)

    if len(mod_var_delta) == 0:
        estr = 'Must provide mod_var_delta as a list of floats.'
        raise ValueError(estr)

    if len(mod_var_delta) != len(mod_name):
        estr = ' '.join(['Must provide the same number of model deltas',
                         'and dimension names.'])
        raise ValueError(estr)

    if len(inst_name) != len(mod_name):
        estr = ' '.join(['Must provide the same number of instrument and model',
                         'location attribute names as a list'])
        raise ValueError(estr)

    if len(mod_name) != len(mod_units):
        raise ValueError(''.join(['Must provide units for each model location',
                                  ' attribute']))

    # Ensure all variables in sel_name exist
    for iname in sel_name:
        if iname not in model:
            raise ValueError(''.join(['Unknown model variable index ',
                                      '{:}'.format(iname)]))

    # Ensure that `mod_irreg_var` exists
    if mod_irreg_var not in model:
        raise ValueError(''.join(['Unknown irregular model variable index ',
                                  '{:}'.format(mod_irreg_var)]))

    # Ensure coordinate dimensions match
    for var in sel_name:
        if model[var].dims != model[mod_irreg_var].dims:
            estr = ' '.join(('Coordinate dimensions must match for',
                             '"mod_irreg_var" and', model[var].name))
            raise ValueError(estr)

    # Ensure mod_reg_dim in mod_irreg_var
    if mod_reg_dim not in model[mod_irreg_var].dims:
        estr = 'mod_reg_dim must be a coordinate dimension for mod_irreg_var.'
        raise ValueError(estr)

    for mname in mod_name:
        if mname not in model[mod_irreg_var].dims:
            estr = ' '.join(('mod_name must contain coordinate dimension',
                             'labels for mod_irreg_var.'))
            raise ValueError(estr)

    # Determine the scaling between model and instrument data
    inst_scale = np.ones(shape=len(inst_name), dtype=float)
    for i, iname in enumerate(inst_name):
        if iname not in inst.data.keys():
            raise ValueError(''.join(['Unknown instrument location index ',
                                      '{:}'.format(iname)]))
        if mod_name[i] != mod_reg_dim:
            # Some units may have extra information (e.g., 'degrees North').
            # Use only the actual units in the scaling function. These are
            # assumed to come first.
            long_units = re.split(r"\W+|_",
                                  inst.meta[iname, inst.meta.labels.units])[0]
            inst_scale[i] = pyutils.scale_units(mod_units[i], long_units)
        else:
            # TODO(#95) update this else to use metadata from
            #  `mod_irreg_var` when made a pysat compatible fuction.
            long_units = re.split(r"\W+|_",
                                  inst.meta[iname, inst.meta.labels.units])[0]
            inst_scale[i] = pyutils.scale_units(mod_units[i], long_units)

            # Store scaling for mod_irreg_var
            dvar_scale = inst_scale[i]

    # Pull out datetime data
    if mod_datetime_name in model.data_vars:
        mod_datetime = model.data_vars[mod_datetime_name]
        mod_datetime = mod_datetime.values.astype(np.int64)
    elif mod_datetime_name in model.coords:
        mod_datetime = model.coords[mod_datetime_name].values.astype(np.int64)
    else:
        raise ValueError("".join(["unknown model name for datetime: ",
                                  mod_datetime_name]))

    # First, model locations for interpolation (regulargrid)
    coords = [model[dim].values / temp_scale
              for dim, temp_scale in zip(mod_name, inst_scale)]

    # Time first
    coords.insert(0, mod_datetime)

    # Translate regular locations to equivalent irregular ones
    # pull out irregular grid locations for variables that will be interpolated
    dvar = model[mod_irreg_var]

    # Make a mesh of data coordinate location values. Start by getting the
    # total number of elements and find which dimension will be updated
    num_pts = 1
    update_dim = -1000
    for i, dim in enumerate(dvar.dims):
        num_pts *= len(dvar.coords[dim])
        if dim == mod_reg_dim:
            update_dim = i

    # Array to store irregular locations of measurements
    points = np.zeros((num_pts, 4))

    # Create mesh corresponding to coordinate values and store
    pts = np.meshgrid(*coords, indexing='ij')
    for i, pt in enumerate(pts):
        points[:, i] = np.ravel(pt)

    # Replace existing regular dimension with irregular data. Model values
    # need to account for any change in units wrt the Instrument.
    points[:, update_dim] = np.ravel(dvar) / dvar_scale

    # Pull out model data to be interpolated. Needed so same downselection
    # procedure may be applied to data as to dimensions.
    model_vars = []
    for iname in sel_name:
        model_vars.append(np.ravel(model[iname].values))

    # Downselect all model points to those within some range of the Instrument
    # values. For each of the model dimensions, each with a corresponding
    # instrument variable, determine the selection criteria and store the
    # limits.
    for i, iname in enumerate(inst_name):
        # Range from Instrument values
        min_val = inst[iname].min()
        max_val = inst[iname].max()

        # Range from model values
        min_mval = np.nanmin(points[:, i + 1])
        max_mval = np.nanmax(points[:, i + 1])

        # Net range
        min_val = np.max([min_val - mod_var_delta[i], min_mval])
        max_val = np.min([max_val + mod_var_delta[i], max_mval])

        # Determine which points are within the specified tolerance range.
        idx, = np.where((points[:, i + 1] >= min_val)
                        & (points[:, i + 1] <= max_val))

        # Downselect model dimension values
        points = points[idx, :]

        # Perform downselection of model variables as well
        for j, ivar in enumerate(model_vars):
            model_vars[j] = ivar[idx]

    ps_mod.logger.debug('Remaining points after downselection {:d}'.format(
        len(idx)))

    # Create input array using inst time/position
    coords = [inst[coord] for coord in inst_name]
    coords.insert(0, inst.index.values.astype(np.int64))
    sat_pts = [inp for inp in zip(*coords)]

    # Perform interpolation of user desired variables
    output_names = []
    for var, idata in zip(sel_name, model_vars):
        ps_mod.logger.debug(''.join(['Creating interpolation object for ',
                                     var]))
        output_names.append('_'.join((model_label, var)))
        inst[output_names[-1]] = \
            interpolate.griddata(points, idata, sat_pts, rescale=True)
        ps_mod.logger.debug(' '.join([var, 'complete.']))

    return output_names


def extract_modelled_observations(inst, model, inst_name, mod_name,
                                  mod_datetime_name, mod_time_name, mod_units,
                                  sel_name=None, time_method='min',
                                  pair_method='closest', method='linear',
                                  model_label='model',
                                  model_units_attr='units'):
    """Extract instrument-aligned data from a modelled data set.

    Parameters
    ----------
    inst : pysat.Instrument
        instrument object for which modelled data will be extracted
    model : xarray.Dataset
        modelled data set
    inst_name : array-like
        list of names of the data series to use for determining instrument
        location
    mod_name : array-like
        list of names of the data series to use for determining model locations
        in the same order as inst_name.  These must make up a regular grid.
    mod_datetime_name : str
        Name of the data series in the model Dataset containing datetime info
    mod_time_name : str
        Name of the time coordinate in the model Dataset
    mod_units : list of strings
        units for each of the mod_name location attributes.  Currently
        supports: rad/radian(s), deg/degree(s), h/hr(s)/hour(s), m, km, and cm
    sel_name : array-like or NoneType
        list of names of modelled data indices to append to instrument object,
        or None to append all modelled data (default=None)
    time_method : str
        Pair data using larger (max) or smaller (min) of the smallest
        instrument/model time increments (default='min')
    pair_method : str
        Find all relevent pairs ('all') or just the closest pairs ('closest').
        (default='closest')
    method : str
        Interpolation method.  Supported are 'linear', 'nearest', and
        'splinef2d'.  The last is only supported for 2D data and is not
        recommended here.  (default='linear')
    model_label : str
        name of model, used to identify interpolated data values in instrument
        (default="model")
    model_units_attr : str
        Attribute for model xarray values that contains units (default='units')

    Returns
    -------
    interp_data.keys() : dict_keys
        Keys of modelled data added to the instrument

    Raises
    ------
    ValueError
        For incorrect input arguments

    Notes
    -----
    For best results, select clean instrument data after alignment with model

    """
    # Ensure the array-like inputs are arrays
    inst_name = np.asarray(inst_name)
    mod_name = np.asarray(mod_name)

    if sel_name is None:
        sel_name = np.asarray(list(model.data_vars.keys()))
    else:
        sel_name = np.asarray(sel_name)
        if sel_name.shape == ():
            sel_name = np.asarray([sel_name])

    # Ensure the method flags are all lower-case for easy testing
    time_method = time_method.lower()
    pair_method = pair_method.lower()

    # Test input
    if len(inst_name) == 0:
        estr = 'Must provide instrument location attribute names as a list'
        raise ValueError(estr)

    if len(inst_name) != len(mod_name):
        estr = 'Must provide the same number of instrument and model '
        estr += 'location attribute names as a list'
        raise ValueError(estr)

    if len(mod_name) != len(mod_units):
        raise ValueError(''.join(['Must provide units for each model location',
                                  ' attribute']))

    if mod_time_name not in model.coords:
        raise ValueError("Unknown model time coordinate key name")

    if time_method not in ['min', 'max']:
        raise ValueError("unknown time method, expects 'min' or 'max'")

    if pair_method not in ['all', 'closest']:
        raise ValueError("unknown pairing method, expects 'all' or 'closest'")

    # Ensure mod_name is a list
    mod_name = list(mod_name)

    # Remove any model coordinates from the modelled data to interpolate
    sel_name = sel_name[[mdat not in mod_name for mdat in sel_name]]

    if len(sel_name) == 0:
        raise ValueError('No model data keys to interpolate')

    # Determine the scaling between model and instrument data
    inst_scale = np.ones(shape=len(inst_name), dtype=float)
    for i, iname in enumerate(inst_name):
        if iname not in inst.data.keys():
            raise ValueError(''.join(['Unknown instrument location index ',
                                      '{:} '.format(iname),
                                      '(should not be epoch time)']))

        # Some units may have extra information (e.g., 'degrees North').
        # Use only the actual units in the scaling function.  These are assumed
        # to come first.
        long_units = re.split(r"\W+|_",
                              inst.meta[iname, inst.meta.labels.units])[0]
        inst_scale[i] = pyutils.scale_units(mod_units[i], long_units)

    # Determine the model time resolution
    if mod_datetime_name in model.data_vars.keys():
        mod_datetime = model.data_vars[mod_datetime_name].values
    elif mod_datetime_name in model.coords:
        mod_datetime = model.coords[mod_datetime_name].values
    else:
        raise ValueError("".join(["unknown model name for datetime: ",
                                  mod_datetime_name]))

    # Determine the appropriate time difference in seconds from the instrument
    # and model data.  If only one time value is present, assume anything is
    # close enough
    dtime = mod_datetime[1:] - mod_datetime[:-1]
    if len(dtime) == 0:
        tm_sec = np.inf
    else:
        tm_sec = dtime.min()
        tm_sec /= np.timedelta64(1, 's')

    # Casting as values creates an array of numpy.timedelta64 objects in ns
    dtime = inst.index.values[1:] - inst.index.values[:-1]
    ti_sec = np.inf if len(dtime) == 0 else dtime.min().astype(float) * 1.0e-9

    # This will still work if infinite, since it will cause all data to be
    # accepted as close enough.
    if time_method == 'max':
        min_del = tm_sec if tm_sec > ti_sec else ti_sec
    else:
        min_del = tm_sec if tm_sec < ti_sec else ti_sec

    # Determine which instrument observations are within the model time
    # resolution of a model run
    mind = list()
    iind = list()
    del_sec = abs(mod_datetime
                  - inst.index.values[:, np.newaxis]).astype(float) * 1.0e-9
    for inst_ind, mod_ind in enumerate(del_sec.argmin(axis=1)):
        if del_sec[inst_ind, mod_ind] <= min_del:
            if mod_ind in mind and pair_method == 'closest':
                # Test to see if this model observation has multiple pairings
                old_ind = mind.index(mod_ind)
                if(del_sec[inst_ind, mod_ind]
                   < del_sec[iind[old_ind], mind[old_ind]]):
                    # If this one is closer, keep it
                    iind[old_ind] = inst_ind
                    mind[old_ind] = mod_ind
            else:
                # If this is a new point, keep it
                iind.append(inst_ind)
                mind.append(mod_ind)

    # Determine the model coordinates closest to the satellite track
    interp_shape = inst.index.shape if inst.pandas_format else \
        [inst.data.sizes[ss] for ss in inst.data.coords.keys()]
    inst_coord = {kk: getattr(inst.data, inst_name[i]).values * inst_scale[i]
                  for i, kk in enumerate(mod_name)}

    # Initalize the interpolated data dictionary and test to ensure that the
    # instrument data doesn't already have the interpolated data
    interp_data = {"{:s}_{:s}".format(model_label, mdat):
                   np.full(shape=interp_shape, fill_value=np.nan)
                   for mdat in sel_name}

    del_list = list()
    for mdat in interp_data.keys():
        if mdat in inst.data.keys():
            ps_mod.logger.warning("".join(["model data already interpolated:",
                                           " {:}".format(mdat)]))
            del_list.append(mdat)

    if len(interp_data.keys()) == len(del_list):
        raise ValueError("instrument object already contains all model data")
    elif len(del_list) > 0:
        for mdat in del_list:
            del interp_data[mdat]

    dim_warned = False
    for i, ii in enumerate(iind):
        # Cycle through each model data type, since it may not depend on
        # all the dimensions
        for mdat in sel_name:
            # Define the output key
            attr_name = "{:s}_{:s}".format(model_label, mdat)

            # Skip if this is already included in the output instrument
            if attr_name not in interp_data.keys():
                break

            # Determine the dimension values
            dims = [kk for kk in model.data_vars[mdat].dims]
            indices = {mod_time_name: mind[i]}

            # Construct the data needed for interpolation, ensuring that
            # the types are appropriate
            values = model[indices][mdat].data
            points = [model.coords[kk].data for kk in dims if kk in mod_name]
            get_coords = True if len(points) > 0 else False
            idims = 0

            if not get_coords and not dim_warned:
                # Only warn the user once
                ps_mod.logger.warning("".join(["model dimensions ",
                                               "{:} don't match".format(dims),
                                               " interpolation parameters ",
                                               "{:}".format(mod_name)]))
                dim_warned = True

            while get_coords:
                if inst.pandas_format:
                    # This data iterates only by time
                    xout = ii
                    xi = [inst_coord[kk][ii] for kk in dims if kk in mod_name]
                    get_coords = False
                else:
                    # This data may have additional dimensions
                    if idims == 0:
                        # Determine the number of dimensions
                        idims = len(inst.data.coords)
                        idim_names = [ckey for i, ckey in
                                      enumerate(inst.data.coords.keys())
                                      if i > 0]

                        # Find relevent dimensions for cycling and slicing
                        ind_dims = [k for k, iname in enumerate(inst_name)
                                    if iname in idim_names]
                        imod_dims = [k for k in ind_dims
                                     if mod_name[k] in dims]
                        coord_keys = [ckey for ckey in inst.data.coords.keys()]
                        ind_dims = [coord_keys.index(inst_name[k])
                                    for k in imod_dims]

                        # Set the number of cycles
                        icycles = 0
                        ncycles = sum([len(inst.data.coords[inst_name[k]])
                                       for k in imod_dims])
                        cinds = np.zeros(shape=len(imod_dims), dtype=int)

                    # Get the instrument coordinate for this cycle
                    if icycles < ncycles or icycles == 0:
                        ss = [ii if k == 0 else 0 for k in range(idims)]
                        se = [ii + 1 if k == 0 else
                              len(inst.data.coords[idim_names[k - 1]])
                              for k in range(idims)]
                        xout = [cinds[ind_dims.index(k)] if k in ind_dims
                                else slice(ss[k], se[k]) for k in range(idims)]
                        xind = [cinds[ind_dims.index(k)] if k in ind_dims
                                else ss[k] for k in range(idims)]
                        xout = tuple(xout)
                        xind = tuple(xind)

                        xi = list()
                        for kk in dims:
                            if kk in mod_name:
                                # This is the next instrument coordinate
                                k = mod_name.index(kk)
                                if k in imod_dims:
                                    # This is an xarray coordiante
                                    cind = imod_dims.index(k)
                                    xi.append(inst_coord[kk][cinds[cind]])
                                else:
                                    # This is an xarray variable
                                    xi.append(inst_coord[kk][xind])

                        # Cycle the indices
                        if len(cinds) > 0:
                            k = 0
                            cinds[k] += 1

                            while cinds[k] > inst.data.coords.dims[
                                    inst_name[imod_dims[k]]]:
                                k += 1
                                if k < len(cinds):
                                    cinds[k - 1] = 0
                                    cinds[k] += 1
                                else:
                                    break
                        icycles += 1

                    # If we have cycled through all the coordinates for this
                    # time, move onto the next time
                    if icycles >= ncycles:
                        get_coords = False

                # Interpolate the desired value
                try:
                    yi = interpolate.interpn(points, values, xi, method=method)
                except ValueError as verr:
                    if str(verr).find("requested xi is out of bounds") > 0:
                        # This is acceptable, pad the interpolated data with
                        # NaN
                        ps_mod.logger.warning(
                            "{:} for {:s} data at {:}".format(verr, mdat, xi))
                        yi = [np.nan]
                    else:
                        raise ValueError(verr)

                # Save the output
                try:
                    interp_data[attr_name][xout] = yi[0]
                except TypeError as terr:
                    ps_mod.logger.error(''.join(['check mod_name input ',
                                                 'against dimensions for ',
                                                 mdat, ' in model data']))
                    raise TypeError(terr)

    # Update the instrument object and attach units to the metadata
    for mdat in interp_data.keys():
        # Assign the units, if keyword is present in xarray data
        attr_name = mdat.split("{:s}_".format(model_label))[-1]
        mod_units = "missing"
        if hasattr(model.data_vars[attr_name], model_units_attr):
            mod_units = getattr(model.data_vars[attr_name], model_units_attr)
        inst.meta[mdat] = {inst.meta.labels.units: mod_units}

        if inst.pandas_format:
            inst[mdat] = pds.Series(interp_data[mdat], index=inst.index)
        else:
            inst.data = inst.data.assign(interp_key=(inst.data.coords.keys(),
                                                     interp_data[mdat]))
            inst.data = inst.data.rename({"interp_key": mdat})

    return interp_data.keys()
