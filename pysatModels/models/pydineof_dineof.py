# -*- coding: utf-8 -*-
"""
Support exported model data from pyDINEOF.

Properties
----------
platform
    'pydineof'
name
    'dineof'
tag
    '', 'test'
inst_id
    ''

Note
----
pyDINEOFs is a Python package that interfaces with a version of
Data Interpolation Empirical Orthogonal Functions (DINEOFs).
This module couples into the systematic export pyDINEOF format and thus
should support all exports from the package.

Specific tags are not listed here as this method is intended to support
all pyDINEOF export models. Place the desired model (daily files)
at '{pysat_data_dir}/pydineof/dineof/{tag}'. Each model series is identified
using the `tag` keyword. It is presumed the default
naming scheme of 'dineof_{year:04d}-{month:02d}-{day:02d}.nc'
has been retained. Use the `file_format` option for custom filenames.

DINEOFs are a purely data based method that can analyze a data-set, with
data gaps, and extract a series of basis functions that optimally reproduce
the input data. The quality of the reconstruction is primarily determined
by the quantity and quality of the input data.

http://modb.oce.ulg.ac.be/mediawiki/index.php/DINEOF

References
----------
J.-M. Beckers and M. Rixen. EOF calculations and data filling from
incomplete oceanographic data sets. Journal of Atmospheric and
Oceanic Technology, 20(12):1839-­1856, 2003.

"""
import datetime as dt
import functools
import os
import requests
import warnings

import pysat

logger = pysat.logger

# ----------------------------------------------------------------------------
# Instrument attributes

platform = 'pydineof'
name = 'dineof'
tags = {'': 'pydineof output file',
        'test': 'Standard output of pydineof for benchmarking'}
inst_ids = {'': [tag for tag in tags.keys()]}

# Specify the use of xarray instead of pandas.
pandas_format = False

# ----------------------------------------------------------------------------
# Instrument test attributes

_test_dates = {'': {tag: dt.datetime(2009, 1, 1) for tag in tags.keys()}}
_test_download = {'': {'': False,
                       'test': True}}

# ----------------------------------------------------------------------------
# Instrument methods


def init(self):
    """Initialize the Instrument object with instrument specific values."""
    acks = ''.join(('The original DINEOF model code may be found at ',
                    'http://modb.oce.ulg.ac.be/mediawiki/index.php/DINEOF.',
                    'pyDINEOFs is stored online in a private repository at ',
                    'https://github.com/PhotonAudioLab/pyDINEOF'))
    self.acknowledgements = acks

    refs = ''.join(('J.-M. Beckers and M. Rixen. EOF calculations and data ',
                    'filling from incomplete oceanographic data sets. ',
                    'Journal of Atmospheric and Oceanic Technology, ',
                    '20(12):1839-­1856, 2003'))
    self.references = refs
    logger.info(self.acknowledgements)
    return


def clean(self):
    """Clean pydineof data to the specified level, unused."""

    logger.info('Cleaning not supported for DINEOFs')

    return


# ----------------------------------------------------------------------------
# Instrument functions
#
# Use local and default pysat methods

# Set the list_files routine
fname = 'dineof_{year:04d}-{month:02d}-{day:02d}.nc'
supported_tags = {'': {'': fname, 'test': fname}}
list_files = functools.partial(pysat.instruments.methods.general.list_files,
                               supported_tags=supported_tags)


def load(fnames, tag=None, inst_id=None, **kwargs):
    """Load pydineof data using xarray.

    Parameters
    ----------
    fnames : array-like
        iterable of filename strings, full path, to data files to be loaded.
        This input is nominally provided by pysat itself.
    tag : str or NoneType
        tag name used to identify particular data set to be loaded.
        This input is nominally provided by pysat itself. (default=None)
    inst_id : str or NoneType
        Instrument ID used to identify particular data set to be loaded.
        This input is nominally provided by pysat itself. (default=None)
    **kwargs : dict
        Pass-through for additional keyword arguments specified when
        instantiating an Instrument object. These additional keywords
        are passed through to this routine by pysat.

    Returns
    -------
    data : xarray.Dataset
        pysat formatted xarray Dataset
    meta : pysat.Metadata
        Model run meta data

    Note
    ----
    Any additional keyword arguments passed to pysat.Instrument
    upon instantiation are passed along to this routine.

    Examples
    --------
    ::

        inst = pysat.Instrument(inst_module=pysatModels.models.pydineof_dineof)
        inst.load(2019, 1)

    """

    # netCDF4 files were produced by xarray.
    # Returning an xarray.Dataset.
    return pysat.utils.load_netcdf4(fnames, epoch_name='time',
                                    pandas_format=False)


def download(date_array=None, tag=None, inst_id=None, data_path=None, **kwargs):
    """Download dineof data.  Currently only retrieves test data from github.

    Parameters
    ----------
    date_array : array-like or NoneType
        list of datetimes to download data for. The sequence of dates need not
        be contiguous. (default=None)
    tag : str or NoneType
        Tag identifier used for particular dataset. This input is provided by
        pysat. (default=None)
    inst_id : str or NoneType
        Instrument ID string identifier used for particular dataset. This input
        is provided by pysat. (default=None)
    data_path : str or NoneType
        Path to directory to download data to. (default=None)
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
        remote_url = 'https://github.com/pysat/pysatModels/'
        remote_path = 'blob/main/pysatModels/tests/test_data/'

        # Need to tell github to show the raw data, not the webpage version.
        fname = 'dineof-2009-01-01.nc?raw=true'

        # Use a pysat-compatible name.
        format_str = supported_tags[inst_id][tag]
        saved_local_fname = os.path.join(data_path,
                                         format_str.format(year=date.year,
                                                           month=date.month,
                                                           day=date.day))
        remote_path = '/'.join((remote_url.strip('/'), remote_path.strip('/'),
                                fname))
        with requests.get(remote_path) as req:
            if req.status_code != 404:
                if not os.path.isfile(saved_local_fname):
                    with open(saved_local_fname, 'wb') as open_f:
                        open_f.write(req.content)
            else:
                warnings.warn('Unable to find remote file: {:}'.format(
                    remote_path))
    else:

        warnings.warn('Downloads currently only supported for test files.')

    return
