#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2022, pysat development team
# Full license can be found in License.md
# -----------------------------------------------------------------------------
"""Routines to match modelled and observational data."""

import datetime as dt
import numpy as np
import pandas as pds

import pysat

import pysatModels
from pysatModels.utils.convert import load_model_xarray
from pysatModels.utils import extract


def collect_inst_model_pairs(start, stop, tinc, inst, inst_download_kwargs=None,
                             model_load_rout=load_model_xarray,
                             model_load_kwargs=None, inst_clean_rout=None,
                             inst_lon_name=None, mod_lon_name=None,
                             lon_pos='end', inst_name=None, mod_name=None,
                             mod_datetime_name=None, mod_time_name=None,
                             mod_units=None, sel_name=None, time_method='min',
                             pair_method='closest', method='linear',
                             model_label='model', comp_clean='clean'):
    """Pair instrument and model data.

    Parameters
    ----------
    start : dt.datetime
        Starting datetime
    stop : dt.datetime
        Ending datetime
    tinc : dt.timedelta
        Time incriment for model files
    inst : pysat.Instrument
        Instrument object for which modelled data will be extracted
    inst_download_kwargs : dict or NoneType
        Optional keyword arguments for downloading instrument data
        (default=None)
    model_load_rout : func
        Routine to load model data into an xarray using datetime as argument
        input input and other necessary data as keyword arguments.  If the
        routine requires a time-dependent filename, ensure that the load
        routine uses the datetime input to construct the correct filename, as
        done in load_model_xarray. (default=load_model_xarray)
    model_load_kwargs : dict or NoneType
        Keyword arguments for the model loading routine. (default=None)
    inst_clean_rout : func
        Routine to clean the instrument data
    inst_lon_name : str
        variable name for instrument longitude
    mod_lon_name : str
        variable name for model longitude
    lon_pos : str or int
        Accepts zero-offset integer for list order or 'end' (default='end')
    inst_name : list or NoneType
        List of names of the data series to use for determing instrument
        location. (default=None)
    mod_name : list or NoneType
        List of names of the data series to use for determing model locations
        in the same order as inst_name.  These must make up a regular grid.
        (default=None)
    mod_datetime_name : str
        Name of the data series in the model Dataset containing datetime info
    mod_time_name : str
        Name of the time coordinate in the model Dataset
    mod_units : list or NoneType
        units for each of the mod_name location attributes.  Currently
        supports: rad/radian(s), deg/degree(s), h/hr(s)/hour(s), m, km, and cm.
        (default=None)
    sel_name : list or NoneType
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
    comp_clean : str
        Clean level for the comparison data ('clean', 'dusty', 'dirty', 'none')
        (default='clean')

    Returns
    -------
    matched_inst : pysat.Instrument
        Instrument object with observational data from `inst` and paired
         modelled data.

    Raises
    ------
    ValueError
        If input is incorrect

    Notes
    -----
    Perform the data cleaning after finding the times and locations where the
    observations and model align.

    """

    # Initialize the output
    matched_inst = None

    # Test the input
    if inst_lon_name is None:
        raise ValueError('Need longitude name for instrument data')

    if mod_lon_name is None:
        raise ValueError('Need longitude name for model data')

    if mod_datetime_name is None:
        raise ValueError('Need datetime coordinate name for model data')

    if mod_time_name is None:
        raise ValueError('Need time coordinate name for model data')

    if inst_name is None or len(inst_name) == 0:
        estr = 'Must provide instrument location attribute names as a list'
        raise ValueError(estr)

    if mod_name is None:
        estr = 'Must provide model location attribute names as a list'
        raise ValueError(estr)

    if mod_units is None:
        raise ValueError('Must provide model units as a list')

    if len(inst_name) != len(mod_name):
        estr = ''.join(['Must provide the same number of instrument and ',
                       'model location attribute names as a list'])
        raise ValueError(estr)

    if len(mod_name) != len(mod_units):
        raise ValueError(''.join(['Must provide units for each model location',
                                  ' attribute']))

    if inst_clean_rout is None:
        raise ValueError('Need routine to clean the instrument data')

    if inst_download_kwargs is None:
        inst_download_kwargs = {}

    if model_load_kwargs is None:
        model_load_kwargs = {}

    skip_download = False
    if "skip_download" in inst_download_kwargs.keys():
        skip_download = inst_download_kwargs['skip_download']
        del inst_download_kwargs['skip_download']

    # Download the instrument data, if needed and wanted
    if not skip_download and (stop
                              - start).days != len(inst.files[start:stop]):
        missing_times = [tt for tt in pds.date_range(start, stop, freq='1D',
                                                     closed='left')
                         if tt not in inst.files[start:stop].index]
        for tt in missing_times:
            inst.download(start=tt, stop=tt + pds.DateOffset(days=1),
                          **inst_download_kwargs)

    # Cycle through the times, loading the model and instrument data as needed
    istart = start
    inst_lon_adjust = True
    inst_dims = []
    while start < stop:
        # Load the model data for each time
        try:
            mdata = model_load_rout(start, **model_load_kwargs)
        except (IOError, ValueError) as err:
            pysatModels.logger.info(
                'unable to load model data at {:}\n{:}'.format(start, err))
            mdata = None

        if mdata is not None:
            # Get the range for model longitude, if it has not already been set
            if inst_lon_adjust:
                if mod_lon_name in mdata.coords:
                    lon_high = float(mdata.coords[mod_lon_name].max())
                    lon_low = float(mdata.coords[mod_lon_name].min())
                elif mod_lon_name in mdata.data_vars:
                    lon_high = float(np.nanmax(mdata.data_vars[mod_lon_name]))
                    lon_low = float(np.nanmin(mdata.data_vars[mod_lon_name]))
                else:
                    raise ValueError("".join(["unknown name for model ",
                                              "longitude: ", mod_lon_name]))

                if lon_high > 180.0 and lon_low < 0.0:
                    raise ValueError("unexpected longitude range")
                elif lon_high > 180.0 or lon_low >= 0.0:
                    lon_low = 0.0
                    lon_high = 360.0
                else:
                    lon_low = -180.0
                    lon_high = 180.0

                # Set the range of the instrument longitude
                inst.custom_attach(pysat.utils.coords.update_longitude,
                                   kwargs={'low': lon_low,
                                           'lon_name': inst_lon_name,
                                           'high': lon_high})
                inst.load(date=istart)

                # Set flag to false now that the range has been set
                inst_lon_adjust = False

            # Load the instrument data, if needed
            if inst.empty or inst.index[-1] < istart:
                inst.load(date=istart)

            if not inst.empty and np.any(inst.index >= istart):
                added_names = extract.extract_modelled_observations(
                    inst=inst, model=mdata, inst_name=inst_name,
                    mod_name=mod_name, mod_datetime_name=mod_datetime_name,
                    mod_time_name=mod_time_name, mod_units=mod_units,
                    sel_name=sel_name, time_method=time_method, method=method,
                    pair_method=pair_method, model_label=model_label)

                if len(added_names) > 0:
                    # Clean the instrument data
                    inst.clean_level = comp_clean
                    inst_clean_rout(inst)

                    im = list()
                    for aname in added_names:
                        # Determine the number of good points
                        imnew = np.where(np.isfinite(inst[aname].values))

                        # Some data types are higher dimensions than others,
                        # make sure we end up choosing a high dimension one
                        # so that we don't accidently throw away paired data
                        if len(im) == 0 or len(im[0]) < len(imnew[0]):
                            im = imnew

                    # If the data is 1D, save it as a list instead of a tuple
                    if len(im) == 1:
                        im = im[0]
                    else:
                        # If the dimension data hasn't been set yet, do it here
                        if len(inst_dims) == 0:
                            inst_dims = [inst.index.name]
                            inst_dims.extend([dd for dd in inst.data.dims.keys()
                                              if dd != inst.index.name])

                        im = {kk: np.unique(im[i])
                              for i, kk in enumerate(inst_dims)}

                    # Save the clean, matched data
                    if matched_inst is None:
                        matched_inst = inst.copy()
                        matched_inst.data = inst[im]
                    else:
                        matched_inst.concat_data(inst[im])

                    # Reset the clean flag
                    inst.clean_level = 'none'

        # Cycle the times
        if tinc.total_seconds() <= 86400.0:
            start += tinc
            if start + tinc > istart + dt.timedelta(days=1):
                istart += dt.timedelta(days=1)
        else:
            if start + tinc >= istart + dt.timedelta(days=1):
                istart += dt.timedelta(days=1)
            if istart >= start + tinc:
                start += tinc

    return matched_inst
