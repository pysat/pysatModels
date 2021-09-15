# -*- coding: utf-8 -*-
"""
Supports loading data from files generated using the main Aether model.
The aether file is a netCDF file with standard and possible custom outputs in
1D, 2D, and 3D (plus time). The Aether project is at
https://www.github.com/AetherModel/aether.

Properties
----------
platform
    'aether'
name
    'main'
tag
    '', 'run_name', 'test'
inst_id
    '3DALL', '3DNEU'

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

platform = 'aether'
name = 'main'
tags = {'': 'aether output file for unnamed run',
        'test': 'Standard output of Aether for benchmarking'}
inst_ids = {inst_id: list(tags.keys())
            for inst_id in ['3dall', '3dneu', '3dion']}

pandas_format = False
multi_file_day = True

# Local attributes
file_cadence = dt.timedelta(minutes=15)
fname = '_'.join(['{inst_id:s}', '{{year:04d}}{{month:02d}}{{day:02d}}',
                  '{{hour:02d}}{{minute:02d}}{{second:02d}}.nc'])
supported_tags = {inst_id: {tag: fname.format(inst_id=inst_id.upper())
                            for tag in tags.keys()}
                  for inst_id in inst_ids.keys()}

aether_epoch = dt.datetime(1965, 1, 1)

# ----------------------------------------------------------------------------
# Instrument test attributes

_test_dates = {inst_id: {tag: dt.datetime(2011, 3, 20) for tag in tags.keys()}
               for inst_id in inst_ids.keys()}
_test_download = {inst_id: {tag: True if tag == 'test' else False
                            for tag in tags.keys()}
                  for inst_id in inst_ids.keys()}

# ----------------------------------------------------------------------------
# Instrument methods


def init(self):
    """Initializes the Instrument object with instrument specific values.
    """

    self.acknowledgements = "".join(["Thanks to the Aether Dev Team"])
    self.references = "".join(["See Aether docs"])
    logger.info(self.acknowledgements)
    return


def clean(self):
    """Method to return Aether data cleaned to the specified level, unused
    """
    logger.info('Cleaning not supported for Aether')
    return


# ----------------------------------------------------------------------------
# Instrument functions
#
# Use local and default pysat methods

# Set the list_files routine
list_files = functools.partial(pysat.instruments.methods.general.list_files,
                               supported_tags=supported_tags,
                               file_cadence=file_cadence)


def load(fnames, tag=None, inst_id=None, **kwargs):
    """Loads Aether data using xarray.

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

    Examples
    --------
    ::

        inst = pysat.Instrument('aether', 'main')
        inst.load(2011, 79)

    """

    # Load data
    data, meta = pysat.utils.load_netcdf4(fnames, pandas_format=pandas_format)
    data['Time'] = [aether_epoch + dt.timedelta(seconds=etime * 1.0e-9)
                    for etime in data['Time'].values.astype(int)]
    data = data.rename({'Time': 'time'})
    
    return data, meta


def download(date_array=None, tag=None, inst_id=None, data_path=None, **kwargs):
    """Downloads aether data.  Currently only retrieves test data from github

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
        is provided by pysat.  (default=None)
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
        remote_url = 'https://github.com/AetherModel/TestData/'
        remote_path = 'blob/main/tests/test1/'

        # Need to tell github to show the raw data, not the webpage version
        local_name = '{:s}.nc?raw=true'.format(
            fname.format(inst_id=inst_id.upper()))

        # Use pysat-compatible name
        format_str = supported_tags[inst_id][tag]
        saved_local_fname = os.path.join(data_path,
                                         format_str.format(year=date.year,
                                                           month=date.month,
                                                           day=date.day,
                                                           hour=date.hour,
                                                           minute=date.minute,
                                                           second=date.second))
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
