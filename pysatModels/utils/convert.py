#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2019, pysat team
# Full license can be found in License.md
# -----------------------------------------------------------------------------
"""
Routines to support extracting pysat.Instrument data as xarray.Datasets

"""
from os import path


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
        or a filename that requires time specification from ftime
        (default=None)

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
    model_xarray = convert_pysat_to_xarray(model_inst)

    return model_xarray


def convert_pysat_to_xarray(inst):
    """ Extract data from a model Instrument as a Dataset with metadata

    Parameters
    ----------
    inst : pysat.Instrument
        Model instrument object

    Returns
    -------
    inst_data : xarray.Dataset or NoneType
        Dataset from pysat Instrument object or None if there is no data

    """

    # Extract the xarray Dataset, returning None if there is no data
    if inst.empty:
        inst_data = None
    elif inst.pandas_format:
        inst_data = inst.data.to_xarray()
    else:
        inst_data = inst.data

    # Add metadata to the xarray Dataset for all labels and data values
    if inst_data is not None:
        for mlabel in inst.meta.data.keys():
            for dkey in inst.meta.data[mlabel].keys():
                if dkey in inst_data.data_vars.keys():
                    inst_data.data_vars[dkey].attrs[mlabel] = inst.meta[dkey,
                                                                        mlabel]

    return inst_data
