# -*- coding: utf-8 -*-
"""
Supports loading data from files generated using the sami2py model.
sami2py file is a netCDF file with multiple dimensions for some variables.
The sami2py project is at https://www.github.com/sami2py/sami2py

Properties
----------
platform
    'sami2py'
name
    'sami2'
tag
    None supported
inst_id
    None supported

"""

import datetime as dt
import functools
import logging
import os
import requests
import warnings

import xarray as xr

import pysat
from pysat.instruments.methods import general as mm_gen


logger = logging.getLogger(__name__)

platform = 'sami2py'
name = 'sami2'

# dictionary of data 'tags' and corresponding description
tags = {'': 'sami2py output file',
        'test': 'Standard output of sami2py for benchmarking'}
inst_ids = {'': ['', 'test']}
_test_dates = {'': {tag: dt.datetime(2019, 1, 1) for tag in tags.keys()}}
_test_download = {'': {'': False,
                       'test': True}}

# specify using xarray (not using pandas)
pandas_format = False

fname = 'sami2py_output_{year:04d}-{month:02d}-{day:02d}.nc'
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

    self.acknowledgements = " ".join(("This work uses the SAMI2 ionosphere",
                                     "model written and developed by the",
                                     "Naval Research Laboratory."))
    self.references = " ".join(("Huba, J.D., G. Joyce, and J.A. Fedder,",
                                     "Sami2 is Another Model of the Ionosphere",
                                     "(SAMI2): A new low‐latitude ionosphere",
                                     "model, J. Geophys. Res., 105, Pages",
                                     "23035-23053,",
                                     "https://doi.org/10.1029/2000JA000035,",
                                     "2000.\n",
                                     "Klenzing, J., Jonathon Smith, Michael",
                                     "Hirsch, & Angeline G. Burrell. (2020,",
                                     "July 17). sami2py/sami2py: Version 0.2.2",
                                     "(Version v0.2.2). Zenodo.",
                                     "http://doi.org/10.5281/zenodo.3950564"))
    logger.info(self.acknowledgements)

    return


# Required method
def clean(self):
    """Method to return SAMI data cleaned to the specified level

    Cleaning level is specified in inst.clean_level and pysat
    will accept user input for several strings. The clean_level is
    specified at instantiation of the Instrument object, though it may be
    updated to a more stringent level and re-applied after instantiation.
    The clean method is applied after default every time data is loaded.

    Note
    ----
    'clean' All parameters should be good, suitable for statistical and
            case studies
    'dusty' All paramers should generally be good though same may
            not be great
    'dirty' There are data areas that have issues, data should be used
            with caution
    'none'  No cleaning applied, routine not called in this case.

    """

    logger.info('Cleaning not supported for SAMI')
    return


def load(fnames, tag=None, inst_id=None, **kwargs):
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
    inst_id : string ('')
        Instrument ID used to identify particular data set to be loaded.
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

    Note
    ----
    Any additional keyword arguments passed to pysat.Instrument
    upon instantiation are passed along to this routine.

    Examples
    --------
    ::

        inst = pysat.Instrument('sami2py', 'sami2')
        inst.load(2019, 1)

    """

    # load data
    data, meta = pysat.utils.load_netcdf4(fnames, pandas_format=False)

    # add time variable for pysat compatibilty
    data['time'] = [dt.datetime(2019, 1, 1)
                    + dt.timedelta(seconds=int(val * 3600.0))
                    for val in data['ut'].values]

    return data, meta


def download(date_array=None, tag=None, inst_id=None, data_path=None, user=None,
             password=None, **kwargs):
    """Downloads sami2py data.  Currently only retrieves test data from github

    Parameters
    ----------
    date_array : array-like
        list of datetimes to download data for. The sequence of dates need not
        be contiguous.
    tag : string
        Tag identifier used for particular dataset. This input is provided by
        pysat. (default='')
    inst_id : string
        Instrument ID string identifier used for particular dataset. This input
        is provided by pysat. (default='')
    data_path : string
        Path to directory to download data to. (default=None)
    user : string
        User string input used for download. Provided by user and passed via
        pysat. If an account is required for dowloads this routine here must
        error if user not supplied. (default=None)
    password : string (None)
        Password for data download.
    **kwargs : dict
        Additional keywords supplied by user when invoking the download
        routine attached to a pysat.Instrument object are passed to this
        routine via kwargs.

    Note
    ----
    This routine is invoked by pysat and is not intended for direct use by
    the end user.

    The test object generates the datetime requested by the user, which may not
    match the date of the model run.

    """

    if tag == 'test':
        date = date_array[0]
        remote_url = 'https://github.com/sami2py/sami2py/'
        remote_path = 'blob/main/sami2py/tests/test_data/'
        # Need to tell github to show the raw data, not the webpage version
        fname = 'sami2py_output.nc?raw=true'
        # Use pysat-compatible name
        format_str = supported_tags[inst_id][tag]
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
        warnings.warn('Downloads currently only supported for test files.')
    return
