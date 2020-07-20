# -*- coding: utf-8 -*-
"""
Supports loading data from files generated using TIEGCM
(Thermosphere Ionosphere Electrodynamics General Circulation Model) model.
TIEGCM file is a netCDF file with multiple dimensions for some variables.

Parameters
----------
platform : string
    'sami2py'
name : string
    'sami2'
tag : string
    None supported
sat_id : string
    None supported

Notes
-----
Loads into xarray format.

"""

# python 2/3 comptability
from __future__ import print_function
from __future__ import absolute_import

import datetime as dt
import functools
import os
import requests
import warnings

import xarray as xr
import pysat
from pysat.instruments.methods import general as mm_gen

platform = 'sami2py'
name = 'sami2'

# dictionary of data 'tags' and corresponding description
tags = {'': 'sami2py output file',
        'test': 'Standard output of sami2py for benchmarking'}
sat_ids = {'': ['']}
_test_dates = {'': {tag: dt.datetime(2019, 1, 1) for tag in tags.keys()}}
_test_download = {'': {'': False,
                       'test': False}}

# specify using xarray (not using pandas)
pandas_format = False

fname = 'sami2py_output_{year:04d}{month:02d}{day:02d}.nc'
supported_tags = {'': {'': fname,
                       'test': fname}}
list_files = functools.partial(mm_gen.list_files,
                               supported_tags=supported_tags)


def init(self):
    """Initializes the Instrument object with instrument specific values.

    Runs once upon instantiation.

    Parameters
    ----------
    self : pysat.Instrument
        This object

    """

    print("".join(["References and information about samip2y are available at ",
                   "https://sami2py.readthedocs.io/en/latest/introduction.html",
                   "#references"]))


def load(fnames, tag=None, sat_id=None, **kwargs):
    """Loads sami2py data using xarray.

    This routine is called as needed by pysat. It is not intended
    for direct user interaction.

    Parameters
    ----------
    fnames : array-like
        iterable of filename strings, full path, to data files to be loaded.
        This input is nominally provided by pysat itself.
    tag : string ('')
        tag name used to identify particular data set to be loaded.
        This input is nominally provided by pysat itself.
    sat_id : string ('')
        Satellite ID used to identify particular data set to be loaded.
        This input is nominally provided by pysat itself.
    **kwargs : extra keywords
        Passthrough for additional keyword arguments specified when
        instantiating an Instrument object. These additional keywords
        are passed through to this routine by pysat.

    Returns
    -------
    data : xarray.Dataset
        pysat formatted xarray Dataset
    metadata : pysat.Metadata
        Model run meta data

    Notes
    -----
    Any additional keyword arguments passed to pysat.Instrument
    upon instantiation are passed along to this routine.

    Examples
    --------
    ::


        inst = pysat.Instrument('sami2py', 'sami2')
        inst.load(2019, 1)

    """

    # load data
    data = xr.open_dataset(fnames[0])
    # move attributes to the Meta object
    # these attributes will be trasnferred to the Instrument object
    # automatically by pysat
    meta = pysat.Meta()
    for attr in data.attrs:
        setattr(meta, attr[0], attr[1])
    data.attrs = []

    # fill Meta object with variable information
    for key in data.variables.keys():
        attrs = data.variables[key].attrs
        meta[key] = attrs

    return data, meta


def download(date_array=None, tag=None, sat_id=None, data_path=None, user=None,
             password=None, **kwargs):
    """Placeholder for UCAR TIEGCM downloads. Doesn't do anything.

    Parameters
    ----------
    date_array : array-like
        list of datetimes to download data for. The sequence of dates need not
        be contiguous.
    tag : string ('')
        Tag identifier used for particular dataset. This input is provided by
        pysat.
    sat_id : string  ('')
        Satellite ID string identifier used for particular dataset. This input
        is provided by pysat.
    data_path : string (None)
        Path to directory to download data to.
    user : string (None)
        User string input used for download. Provided by user and passed via
        pysat. If an account
        is required for dowloads this routine here must error if user not
        supplied.
    password : string (None)
        Password for data download.
    **kwargs : dict
        Additional keywords supplied by user when invoking the download
        routine attached to a pysat.Instrument object are passed to this
        routine via kwargs.

    Notes
    -----
    This routine is invoked by pysat and is not intended for direct use by
    the end user.

    The test object generates the datetime requested by thte user, which may not
    match the date of the model run.

    """

    if tag == 'test':
        date = date_array[0]
        remote_url = 'https://github.com/sami2py/sami2py/'
        remote_path = 'blob/main/sami2py/tests/test_data/'
        # Need to tell github to show the raw data, not the webpage version
        fname = 'sami2py_output.nc?raw=true'
        # Use pysat-compatible name
        format_str = supported_tags[sat_id][tag]
        saved_local_fname = os.path.join(data_path,
                                         format_str.format(year=date.year,
                                                           month=date.month,
                                                           day=date.day))
        remote_path = '/'.join((remote_url.strip('/'), remote_path.strip('/'),
                                fname))
        req = requests.get(remote_path)
        if req.status_code != 404:
            open(saved_local_fname, 'wb').write(req.content)

    else:
        warnings.warn('Not implemented in this version.')
    return
