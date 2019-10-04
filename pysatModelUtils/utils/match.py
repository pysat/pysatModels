#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2019, AGB & pysat team
# Full license can be found in License.md
#-----------------------------------------------------------------------------
"""
Routines to match modelled and observational data

Routines
--------
collect_inst_model_pairs

"""

from __future__ import absolute_import
from __future__ import unicode_literals

import datetime as dt
import numpy as np
from os import path

import pysat

from . import extract

def load_model_xarray(ftime, model_inst=None, filename=None):
    """ Load and extract data from a model Instrument at the specified time

    Parameters
    ----------
    ftime : dt.datetime
        Desired time for model Instrument input
    model_inst : pysat.Instrument
        Model instrument object
    filename : str or NoneType
        Model filename, if the file is not include in the Instrument filelist.
        or a filename that requires time specification from ftime (default=None)

    Returns
    -------
    model_xarray : xarray.Dataset or NoneType
        Dataset from pysat Instrument object or None if there is no data

    """
    # Test the input
    if not hasattr(model_inst, 'load'):
        raise ValueError('must provide a pysat.Instrument object')

    # Format the filename, if needed
    if hasattr(ftime, 'strftime') and filename is not None:
        fname = ftime.strftime(filename)
    else:
        fname = filename

    # Load the model data, using the file if it exists
    if fname is not None and path.isfile(fname):
        model_inst.load(fname=fname)
    else:
        model_inst.load(date=ftime)

    # Extract the xarray Dataset, returning None if there is no data
    if model_inst.empty:
        model_xarray = None
    elif model_inst.pandas_format:
        model_xarray = model_inst.data.to_xarray()
    else:
       model_xarray = model_inst.data

    return model_xarray

    
def collect_inst_model_pairs(start, stop, tinc, inst, inst_download_kwargs={},
                             model_load_rout=load_model_xarray,
                             model_load_kwargs={"model_inst": None},
                             inst_clean_rout=None, inst_lon_name=None,
                             mod_lon_name=None, inst_name=[], mod_name=[],
                             mod_datetime_name=None, mod_time_name=None,
                             mod_units=[], sel_name=None, method='linear',
                             model_label='model', comp_clean='clean'):
    """Pair instrument and model data

    Parameters
    ----------
    start : dt.datetime
        Starting datetime
    stop : dt.datetime
        Ending datetime
    tinc : dt.timedelta
        Time incriment for model files
    inst : pysat.Instrument instance
        instrument object for which modelled data will be extracted
    inst_download_kwargs : dict
        optional keyword arguments for downloading instrument data (default={})
    model_load_rout : routine
        Routine to load model data into an xarray using datetime as arguement
        input input and other necessary data as keyword arguments.  If the
        routine requires a filename, ensure that the routine uses the datetime
        input to construct the correct filename, such as 'model_%Y%j.nc'
        (default=load_model_xarray)
    model_load_kwargs : dict
        string format that will construct the desired model filename from a
        datetime object.  The default will fail unless a model Instrument object
        is provided.  (default={"model_inst": None})
    inst_clean_rout : routine
        Routine to clean the instrument data
    inst_lon_name : string
        variable name for instrument longitude
    mod_lon_name : string
        variable name for model longitude
    inst_name : list of strings
        list of names of the data series to use for determing instrument
        location
    mod_name : list of strings
        list of names of the data series to use for determing model locations
        in the same order as inst_name.  These must make up a regular grid.
    mod_datetime_name : string
        Name of the data series in the model Dataset containing datetime info
    mod_time_name : string
        Name of the time coordinate in the model Dataset
    mod_units : list of strings
        units for each of the mod_name location attributes.  Currently
        supports: rad/radian(s), deg/degree(s), h/hr(s)/hour(s), m, km, and cm
    sel_name : list of strings or NoneType
        list of names of modelled data indices to append to instrument object,
        or None to append all modelled data (default=None)
    method : string
        Interpolation method.  Supported are 'linear', 'nearest', and
        'splinef2d'.  The last is only supported for 2D data and is not
        recommended here.  (default='linear')
    model_label : string
        name of model, used to identify interpolated data values in instrument
        (default="model")
    comp_clean : string
        Clean level for the comparison data ('clean', 'dusty', 'dirty', 'none')
        (default='clean')

    Returns
    -------
    matched_inst : pysat.Instrument instance
        instrument object and paired modelled data

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

    if len(inst_name) == 0:
        estr = 'Must provide instrument location attribute names as a list'
        raise ValueError(estr)

    if len(inst_name) != len(mod_name):
        estr = ''.join(['Must provide the same number of instrument and model ',
                       'location attribute names as a list'])
        raise ValueError(estr)

    if len(mod_name) != len(mod_units):
        raise ValueError('Must provide units for each model location attribute')

    if inst_clean_rout is None:
        raise ValueError('Need routine to clean the instrument data')

    # Download the instrument data, if needed
    # Could use some improvement, for not re-downloading times that you already
    # have
    if (stop - start).days != len(inst.files[start:stop]):
        inst.download(start=start, stop=stop, **inst_download_kwargs)

    # Cycle through the times, loading the model and instrument data as needed
    istart = start
    while start < stop:
        mdata = model_load_rout(start, **model_load_kwargs)

        if mdata is not None:
            # Get the range for model longitude
            if mod_lon_name in mdata.coords:
                lon_high = float(mdata.coords[mod_lon_name].max())
                lon_low = float(mdata.coords[mod_lon_name].min())
            elif mod_lon_name in mdata.data_vars:
                lon_high = float(np.nanmax(mdata.data_vars[mod_lon_name]))
                lon_low = float(np.nanmin(mdata.data_vars[mod_lon_name]))
            else:
                raise ValueError("".join(["unknown name for model longitude: ",
                                          mod_lon_name]))
            
            # Load the instrument data, if needed
            if inst.empty or inst.index[-1] < istart:
                inst.custom.add(pysat.utils.coords.update_longitude, 'modify',
                                low=lon_low, lon_name=inst_lon_name,
                                high=lon_high)
                inst.load(date=istart)

            if not inst.empty and inst.index[0] >= istart:
                added_names = extract.extract_modelled_observations(inst=inst,
                                        model=mdata, inst_name=inst_name,
                                        mod_name=mod_name,
                                        mod_datetime_name=mod_datetime_name,
                                        mod_time_name=mod_time_name,
                                        mod_units=mod_units, sel_name=sel_name,
                                        method=method, model_label=model_label)

                if len(added_names) > 0:
                    # Clean the instrument data
                    inst.clean_level = comp_clean
                    inst_clean_rout(inst)

                    im = list()
                    for aname in added_names:
                        # Determine the number of good points
                        if inst.pandas_format:
                            imnew = np.where(np.isfinite(inst[aname]))
                        else:
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
                        im = {kk: np.unique(im[i])
                              for i, kk in enumerate(inst.data.coords.keys())}

                    # Save the clean, matched data
                    if matched_inst is None:
                        matched_inst = pysat.Instrument
                        matched_inst.meta = inst.meta
                        matched_inst.data = inst[im]
                    else:
                        idata = inst[im]
                        matched_inst.data = \
                            inst.concat_data([matched_inst.data, idata])

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

    # Recast as xarray and add units
    if matched_inst is not None:
        if inst.pandas_format:
            matched_inst.data = matched_inst.data.to_xarray()
        for im in inst.meta.data.units.keys():
            if im in matched_inst.data.data_vars.keys():
                matched_inst.data.data_vars[im].attrs['units'] = \
                    inst.meta.data.units[im]

    return matched_inst

