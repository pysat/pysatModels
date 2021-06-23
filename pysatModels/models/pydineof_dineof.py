# -*- coding: utf-8 -*-
"""
Supports exported model data from pyDINEOF, a Python package
that interfaces with a version of Data Interpolation Empirical
Orthogonal Functions (DINEOFs). This module couples into
the systematic export pyDINEOF format and thus should support
all exports from the package.

Given the generality of the support, each model series is identified
using the `tag` keyword.

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

Properties
----------
platform : string
    pysat
name : string
    dineof
inst_id : string
    ['']
tag : string
    [*]

Note
----
Specific tags are not listed here as this method is intended to support
all pyDINEOF export models. Place the desired model (daily files)
at '{pysat_data_dir}/pydineof/dineof/{tag}'. It is presumed the default
naming scheme of 'dineof_{year:04d}-{month:02d}-{day:02d}.nc'
has been retained. Use the file_format option for custom filenames.
::
    imodule = pysatModels.instruments.pydineof_dineof
    model = pysat.Instrument(inst_module=imodule, tag=tag)

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
inst_ids = {'': ['', 'test']}

# specify using xarray (not using pandas)
pandas_format = False

# ----------------------------------------------------------------------------
# Instrument test attributes

_test_dates = {'': {tag: dt.datetime(2009, 1, 1) for tag in tags.keys()}}
_test_download = {'': {'': False,
                       'test': True}}

# ----------------------------------------------------------------------------
# Instrument methods


def init(self):
    """Initializes the Instrument object with instrument specific values.

    Runs once upon instantiation.

    Parameters
    ----------
    self : pysat.Instrument
        This object

    """
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


# Required method
def clean(self):
    """Method to return pydineof data cleaned to the specified level

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

    logger.info('Cleaning not supported for DINEOFs')

    return


# ----------------------------------------------------------------------------
# Instrument functions
#
# Use local and default pysat methods

# Set the list_files routine
# Set the list_files routine
fname = 'dineof_{year:04d}-{month:02d}-{day:02d}.nc'
supported_tags = {'': {'': fname, 'test': fname}}
list_files = functools.partial(pysat.instruments.methods.general.list_files,
                               supported_tags=supported_tags)


def load(fnames, tag=None, inst_id=None, **kwargs):
    """Loads pydineof data using xarray.

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
    meta : pysat.Metadata
        Model run meta data

    Note
    ----
    Any additional keyword arguments passed to pysat.Instrument
    upon instantiation are passed along to this routine.

    ::

        inst = pysat.Instrument(inst_module=pysatModels.models.pydineof_dineof)
        inst.load(2019, 1)

    """

    # netCDF4 files were produced by xarray
    # returning an xarray.Dataset
    return pysat.utils.load_netcdf4(fnames, epoch_name='time',
                                    pandas_format=False)


def download(date_array=None, tag=None, inst_id=None, data_path=None, **kwargs):
    """Downloads dineof data.  Currently only retrieves test data from github

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

        # Need to tell github to show the raw data, not the webpage version
        fname = 'dineof-2009-01-01.nc?raw=true'

        # Use a pysat-compatible name
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
            warnings.warn('Unable to find remote file: {:}'.format(remote_path))
    else:
        warnings.warn('Downloads currently only supported for test files.')

    return
