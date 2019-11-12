#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2019, AGB & pysat team
# Full license can be found in License.md
#-----------------------------------------------------------------------------
"""
Routines to extract observational-style data from model output

Routines
--------
satellite_view_through_model
extract_modelled_observations

"""

from __future__ import absolute_import
from __future__ import unicode_literals

import numpy as np
import pandas as pds
import scipy.interpolate as interpolate

import pysat.utils as pyutils

import pysatModelUtils as pysat_mu


def instrument_altitude_to_model_pressure(obs, mod, obs_coords,
                                          obs_alt, mod_alt, mod_coord,
                                          scale=100.,
                                          obs_out='model_altitude',
                                          tol=1.,
                                          autoscale=False):
    """Interpolates altitude values onto model pressure levels.

    Uses a recursive regular grid interpolation to find the
    appropriate pressure level for the given input locations.

    Parameters
    ----------
    obs : pysat.Instrument object
        Instrument object with observational data
    mod : pysat.Instrument object
        Instrument object with modelled data
    obs_coords : array-like
        List of variable names containing the observational data coordinates
        at which the model data will be interpolated. Must be in the same order
        as the model dimensions.
    obs_alt : string
        String identifier used in inst for the altitude variable
    mod_alt : string
        Variable identifier for altitude data in the model
        e.g. 'ZG' in standard TIEGCM files.
    mod_coord : string
        Coordinate dimension identifier for the closest pressure
        level in model, equivalent in a sense to altitude
        e.g. 'ilev', 'lev' for TIEGCM
    scale : float
        Scalar used to roughly translate a change in altitude with a
        change in pressure level, the scale height.
    obs_out : string
        Label used to add interpolated altitudes into pbs Instrument
        object.
    tol : float
        Allowed difference between observed and modeled altitudes.
    autoscale : bool
        if True, use 'units' in meta to determine relative scaling
        between obs and mod

    """

    # Ensure the coordinate and data variable names are array-like
    obs_coords = np.asarray(obs_coords)

    # create initial fake regular grid index in inst
    obs[mod_coord] = 0

    # we need to create altitude index from model
    # collect relevant inputs
    # First, model locations for interpolation
    # we use the dimensions associated with model altitude
    # in the order provided
    points = [mod.data.coords[dim].values if dim != 'time' else
              mod.data.coords[dim].values.astype(int)
              for dim in mod[mod_alt].dims]
    mod_units = [mod.meta[dim]['units'] for dim in mod[mod_alt].dims]
    # drop time units
    mod_units = mod_units[1:]

    # Determine the scaling between model and instrument data
    inst_scale = np.ones(shape=len(obs_coords), dtype=float)
    if autoscale:
        for i, iname in enumerate(obs_coords):
            if iname not in obs.variables:
                raise ValueError(''.join(['Unknown instrument location index ',
                                        '{:}'.format(iname)]))
            if iname == obs_alt:
                # pull out altitude units
                inst_scale[i] = pyutils.scale_units(obs.meta[iname]['units'],
                                                    mod.meta[mod_alt]['units'])

                # pull out altitude unit scalar
                alt_scale = inst_scale[i]
            else:
                inst_scale[i] = pyutils.scale_units(obs.meta[iname]['units'],
                                        mod_units[i])
    else:
        alt_scale = 1

    # create interpolator
    interp = interpolate.RegularGridInterpolator(points,
                                                 np.log(mod[mod_alt].values*alt_scale),
                                                 bounds_error=False,
                                                 fill_value=None)
    # use this interpolator to figure out what altitudes we are at
    # each loop, use the current estimated path through model expressed in pressure
    # and interp above to get the equivalent altitude of this path
    # compare this altitude to the actual instrument altitude
    # shift the equivalent pressure for the instrument up/down
    # until difference ibetween altitudes is small
    # currently fixed to 10 iterations

    # log of instrumeent altitude
    log_ialt = np.log(obs[obs_alt])
    diff = log_ialt*0 + 1000.
    while np.any(np.abs(diff) > tol):
        # create input array using satellite time/position
        # replace the altitude coord with the fake tiegcm one
        coords = []
        for iscale, coord in zip(inst_scale, obs_coords):
            if coord == obs_alt:
                # don't scale altitude-like model coordinate
                coords.append(obs[mod_coord])
            else:
                # scale other dimensions to the model
                coords.append(obs[coord]/iscale)

        coords.insert(0, obs.index.values.astype(int))
        # to peform the interpolation we need points
        # like (x1, y1, z1), (x2, y2, z2)
        # but we currently have a list like
        # [(x1, x2, x3), (y1, y2, y3), ...]
        # convert to the form we need
        # the * below breaks each main list out, and zip
        # repeatedly outputs a tuple (xn, yn, zn)
        sat_pts = [inp for inp in zip(*coords)]

        # altitude pulled out from model
        orbit_alt = interp(sat_pts)
        # difference in altitude
        diff = np.e**orbit_alt - np.e**log_ialt
        # shift index in inst for model pressure level
        # in the opposite direction to diff
        obs[mod_coord] -= diff/scale
    obs[obs_out] = np.e**orbit_alt


def instrument_view_through_model(obs, mod, obs_coords, mod_data_names):
    """Interpolates model values onto instrument locations.

    A RegularGrid interpolation is used for quick performance. This
    method may require the use of a pre-processor on coordinate
    dimensions to ensure that a regular interpolation may actually be
    performed. Exponential interpolation along the height direction
    is optionally supported.

    Parameters
    ----------
    obs : pysat.Instrument object
        Instrument object with observational data
    mod : pysat.Instrument object
        Instrument object with modelled data
    obs_coords : array-like
        List of variable names containing the observational data coordinates
        at which the model data will be interpolated. Do not include 'time',
        only spatial coordinates.
    mod_data_names : array-like
        List of model data output variable names to interpolate. All
        variables require the same dimensions.

    Notes
    -----
    Updates the obs Instrument with interpolated data from the mod Instrument.
    The interpolation is performed via the RegularGridInterpolator.

    Models, such as TIEGCM, have a regular grid in pressure, not in altitude.
    To use this routine for TIEGCM please use
        instrument_altitude_to_model_pressure
    first to transform instrument altitudes to pressure levels suitable for
    use in this method.

    """

    # Ensure the coordinate and data variable names are array-like
    obs_coords = np.asarray(obs_coords)
    mod_data_names = np.asarray(mod_data_names)

    # create obs input based upon provided dimension names
    coords = [obs[coord_name] for coord_name in obs_coords]
    # time goes first
    coords.insert(0, obs.index.values.astype(int))
    # move from a list of lists [ [x1, x2, ...], [y1, y2, ...]]
    # to a list of tuples
    # [(x1, y1, ...), (x2, y2, ...)]
    # required for regulargrid interpolator
    obs_pts = [inp for inp in zip(*coords)]

    # collect model grid points together
    points = []
    for dim in mod[mod_data_names[0]].dims:
        if dim == 'time':
            points.append(mod.data.coords[dim].values.astype(int) )
        else:
            points.append(mod.data.coords[dim].values)

    # perform the interpolation
    interp = {}
    for label in mod_data_names:
        # create the interpolator
        interp[label] = interpolate.RegularGridInterpolator(points,
                                                            mod[label].values,
                                                            bounds_error=False,
                                                            fill_value=None)
        # apply it at oberved locations
        obs[''.join(('model_', label))] = interp[label](obs_pts)
        # get meta data as well
        obs.meta[''.join(('model_', label))] = mod.meta[label]


def instrument_view_irregular_model(sat, model, dim1, dim_var, scoords, new_vars):
    """Interpolate model from irregular to regular sampling.

    Parameters
    ----------
    sat : pysat.Instrument object
        Satellite object that will recieve interpolated data based upon position
    model : pysat.Instrument object
        Model object that will be interpolated onto satellite path
    dim1 : string identifier
        Existing regular dimension to be replaced
    dim2 : string identifier
        Existing irregular variable used to define regular grid

    """

    import scipy.interpolate

    # create inputs for interpolation
    dvar = model[dim_var]

    # make a mesh of data location values using intrinsic
    # regular grid
    num_pts = 1
    coords = []
    update_dim = -1000
    for i, dim in enumerate(dvar.dims):
        num_pts *= len(dvar.coords[dim])
        if dim == 'time':
            coords.append(model.data.coords[dim].values.astype(int))
        else:
            coords.append(model.data.coords[dim].values)
        if dim == dim1:
            update_dim = i

    # locations of measurements
    points = np.zeros((num_pts, 4))

    pts = np.meshgrid(*coords, indexing='ij')
    for i, pt in enumerate(pts):
        points[:,i] = np.ravel(pt)

    # replace existing regular dimension with irregular data
    points[:, update_dim] = np.ravel(dvar)

    # downselect points to those in altitude range of satellite
    # print(points[:, update_dim])
    # print(sat['altitude'].min(), sat['altitude'].max())

    min_val = sat['altitude'].min() - 20. if sat['altitude'].min() < np.nanmax(points[:, update_dim]) else np.nanmax(points[:, update_dim] - 20.)
    max_val = sat['altitude'].max() + 20. if sat['altitude'].max() < np.nanmax(points[:, update_dim]) else np.nanmax(points[:, update_dim])

    idx, = np.where((points[:, update_dim] >= min_val) &
                    (points[:, update_dim] <= max_val))
    points = points[idx, :]
    print ('Remaining points after downselection', len(idx))
    print (points[:, update_dim], np.nanmin(points[:, update_dim]), np.nanmax(points[:, update_dim]))
    # create input array using satellite time/position
    coords = [sat[coord] for coord in scoords]
    coords.insert(0, sat.index.values.astype(int))
    sat_pts = [inp for inp in zip(*coords)]

    interp = {}
    for var in new_vars:
        print('Creating interpolation object for', var, '.')
        # interp[var] = scipy.interpolate.LinearNDInterpolator(points,
        #                                                      np.ravel(model[var].values),
        #                                                      rescale=True)
        sat[''.join(('model_', var))] = scipy.interpolate.griddata(points,
                                                 np.ravel(model[var].values)[idx],
                                                 sat_pts,
                                                 rescale=True)
        print('Interpolating', var, 'onto satellite.')
        # sat[''.join(('model_', var))] = interp[var](sat_pts)
    print('Complete.')



# Needs a better name, is this being used anywhere?
def instrument_view_through_model(obs, mod, obs_coords, mod_dat_names):
    """Interpolate model values onto satellite orbital path.

    Parameters
    ----------
    obs : pysat.Instrument object
        Instrument object with observational data
    mod : pysat.Instrument object
        Instrument object with modelled data
    obs_coords : array-like
        List of variable names containing the observational data coordinates
        at which the model data will be interpolated
    mod_dat_names : array-like
        List of model data output variable names  to interpolate

    Notes
    -----
    Updates the obs Instrument with interpolated data from the mod Instrument

    """
    # Ensure the coordinate and data variable names are array-like
    obs_coords = np.asarray(obs_coords)
    mod_dat_names = np.asarray(mod_dat_names)

    # Create input array using observational data's time/position
    # This needs to be changed, pretty sure it doesn't work for xarray data
    pysat_mu.logger.debug("the coordinate data section needs to be fixed")
    coords = [obs.data[cc] for cc in obs_coords]
    coords.insert(0, obs.index.values.astype(int))
    obs_pts = [inp for inp in zip(*coords)] # what is this doing?

    # Add optional scaling?

    # Interpolate each model data value onto the observations time and location
    for label in mod_dat_names:
        points = [mod.data.coords[dim].values if dim != 'time' else
                  mod.data.coords[dim].values.astype(int)
                  for dim in mod[label].dims]
        interp_val = interpolate.RegularGridInterpolator(points,
                                                         mod[label].values,
                                                         bounds_error=False,
                                                         fill_value=None)
        obs[''.join(('model_', label))] = interp_val(obs_pts)

    # Update the observation's meta data
    pysat_mu.logger.debug("Missing meta data update")

    return



def extract_modelled_observations(inst, model, inst_name, mod_name,
                                  mod_datetime_name, mod_time_name, mod_units,
                                  sel_name=None, method='linear',
                                  model_label='model'):
    """Extracts instrument-aligned data from a modelled data set

    Parameters
    ----------
    inst : pysat.Instrument instance
        instrument object for which modelled data will be extracted
    model : xarray Dataset
        modelled data set
    inst_name : array-like
        list of names of the data series to use for determing instrument
        location
    mod_name : array-like
        list of names of the data series to use for determing model locations
        in the same order as inst_name.  These must make up a regular grid.
    mod_datetime_name : string
        Name of the data series in the model Dataset containing datetime info
    mod_time_name : string
        Name of the time coordinate in the model Dataset
    mod_units : list of strings
        units for each of the mod_name location attributes.  Currently
        supports: rad/radian(s), deg/degree(s), h/hr(s)/hour(s), m, km, and cm
    sel_name : array-like or NoneType
        list of names of modelled data indices to append to instrument object,
        or None to append all modelled data (default=None)
    method : string
        Interpolation method.  Supported are 'linear', 'nearest', and
        'splinef2d'.  The last is only supported for 2D data and is not
        recommended here.  (default='linear')
    model_label : string
        name of model, used to identify interpolated data values in instrument
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
    --------
    For best results, select clean instrument data after alignment with model

    """

    # Ensure the array-like inputs are arrays
    inst_name = np.asarray(inst_name)
    mod_name = np.asarray(mod_name)

    if sel_name is None:
        sel_name = np.asarray(list(model.data_vars.keys()))
    else:
        sel_name = np.asarray(sel_name)

    # Test input
    if len(inst_name) == 0:
        estr = 'Must provide instrument location attribute names as a list'
        raise ValueError(estr)

    if len(inst_name) != len(mod_name):
        estr = 'Must provide the same number of instrument and model '
        estr += 'location attribute names as a list'
        raise ValueError(estr)

    if len(mod_name) != len(mod_units):
        raise ValueError('Must provide units for each model location ' +
                         'attribute')

    if mod_time_name not in model.coords:
        raise ValueError("Unknown model time coordinate key name")

    # Remove any model coordinates from the modelled data to interpolate
    sel_name = sel_name[[mdat not in mod_name for mdat in sel_name]]

    if len(sel_name) == 0:
        raise ValueError('No model data keys to interpolate')

    # Determine the scaling between model and instrument data
    inst_scale = np.ones(shape=len(inst_name), dtype=float)
    for i, iname in enumerate(inst_name):
        if iname not in inst.data.keys():
            raise ValueError(''.join(['Unknown instrument location index ',
                                      '{:}'.format(iname)]))
        inst_scale[i] = pyutils.scale_units(mod_units[i],
                                            inst.meta.data.units[iname])

    # Determine the model time resolution
    if mod_datetime_name in model.data_vars:
        mod_datetime = model.data_vars[mod_datetime_name].values
    elif mod_datetime_name in model.coords:
        mod_datetime = model.coords[mod_datetime_name].values
    else:
        raise ValueError("".join(["unknown model name for datetime: ",
                                  mod_datetime_name]))

    tm_sec = (mod_datetime[1:] - mod_datetime[:-1]).min()
    tm_sec /= np.timedelta64(1, 's')
    ti_sec = (inst.index[1:] - inst.index[:-1]).min().total_seconds()
    min_del = tm_sec if tm_sec < ti_sec else ti_sec

    # Determine which instrument observations are within the model time
    # resolution of a model run
    mind = list()
    iind = list()
    for i, tt in enumerate(mod_datetime):
        del_sec = abs(tt - inst.index).total_seconds()
        if del_sec.min() <= min_del:
            iind.append(del_sec.argmin())
            mind.append(i)

    # Determine the model coordinates closest to the satellite track
    interp_shape = inst.index.shape if inst.pandas_format else \
        inst.data.data_vars.items()[0][1].shape
    inst_coord = {kk: getattr(inst.data, inst_name[i]).values * inst_scale[i]
                  for i, kk in enumerate(mod_name)}

    # Initalize the interpolated data dictionary and test to ensure that the
    # instrument data doesn't already have the interpolated data
    interp_data = {"{:s}_{:s}".format(model_label, mdat):
                   np.full(shape=interp_shape, fill_value=np.nan)
                   for mdat in sel_name}

    for mdat in interp_data.keys():
        if mdat in inst.data.keys():
            pysat_mu.logger.warning("".join(["model data already interpolated:",
                                             " {:}".format(mdat)]))
            del interp_data[mdat]

    if len(interp_data.keys()) == 0:
        raise ValueError("instrument object already contains all model data")

    for i, ii in enumerate(iind):
        # Cycle through each model data type, since it may not depend on
        # all the dimensions
        for mdat in sel_name:
            # Define the output key
            attr_name = "{:s}_{:s}".format(model_label, mdat)

            # Determine the dimension values
            dims = list(model.data_vars[mdat].dims)
            ndim = model.data_vars[mdat].data.shape
            indices = {mod_time_name: mind[i]}

            # Construct the data needed for interpolation
            values = model[indices][mdat].data
            points = [model.coords[kk].data for kk in dims if kk in mod_name]
            get_coords = True if len(points) > 0 else False
            idims = 0

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
                        idim_names = inst.data.coords.keys()[1:]

                        # Find relevent dimensions for cycling and slicing
                        ind_dims = [k for k, kk in enumerate(inst_name)
                                    if kk in idim_names]
                        imod_dims = [k for k in ind_dims
                                     if mod_name[k] in dims]
                        ind_dims = [inst.data.coords.keys().index(inst_name[k])
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
                              len(inst.data.coords[idim_names[k-1]])
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
                                    xi.append(inst_coord[kk][cinds[k]])
                                else:
                                    # This is an xarray variable
                                    xi.append(inst_coord[kk][xind])

                        # Cycle the indices
                        if len(cinds) > 0:
                            k = 0
                            cinds[k] += 1

                            while cinds[k] > \
                                inst.data.coords.dims[inst_name[imod_dims[k]]]:
                                k += 1
                                if k < len(cinds):
                                    cinds[k-1] = 0
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
                        print("Warning: {:} for ".format(verr) +
                              "{:s} data at {:}".format(mdat, xi))
                        yi = [np.nan]
                    else:
                        raise ValueError(verr)

                # Save the output
                interp_data[attr_name][xout] = yi[0]

    # Update the instrument object and attach units to the metadata
    for mdat in interp_data.keys():
        attr_name = mdat.split("{:s}_".format(model_label))[-1]
        inst.meta[mdat] = {inst.units_label: model.data_vars[attr_name].units}

        if inst.pandas_format:
            inst[mdat] = pds.Series(interp_data[mdat], index=inst.index)
        else:
            inst.data = inst.data.assign(interp_key=(inst.data.coords.keys(),
                                                     interp_data[mdat]))
            inst.data.rename({"interp_key": mdat}, inplace=True)

    return interp_data.keys()
