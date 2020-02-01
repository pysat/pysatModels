# -*- coding: utf-8 -*-
"""
Produces fake instrument data for testing.
"""
from __future__ import print_function
from __future__ import absolute_import

import os
import numpy as np
import pandas as pds
import xarray as xr

import pysat

platform = 'pysat'
name = 'testmodel'

pandas_format = False
_test_dates = {'': {'': pysat.datetime(2009, 1, 1)}}


def init(self):
    self.new_thing = True


def load(fnames, tag=None, sat_id=None, malformed_index=False):
    """ Loads the test files

    Parameters
    ----------
    fnames : (list)
        List of filenames
    tag : (str or NoneType)
        Instrument tag (accepts '')
    sat_id : (str or NoneType)
        Instrument satellite ID (accepts '' or a number (i.e., '10'), which
        specifies the number of data points to include in the test instrument)
    malformed_index : bool (False)
        If True, the time index will be non-unique and non-monotonic.

    Returns
    -------
    data : (xr.Dataset)
        Testing data
    meta : (pysat.Meta)
        Metadata

    """

    # create an artifical model data set
    parts = os.path.split(fnames[0])[-1].split('-')
    yr = int(parts[0])
    month = int(parts[1])
    day = int(parts[2][0:2])
    date = pysat.datetime(yr, month, day)

    # Create one day of data at desired frequency
    index = pds.date_range(start=date, end=date+pds.DateOffset(seconds=86399),
                           freq='9000s')
    # Allow numeric string to select first set of data
    if sat_id.isnumeric() and (int(sat_id) < 86400):
        index = index[0:int(sat_id)]

    uts = index.hour*3600 + index.minute*60 + index.second

    latitude = np.linspace(-50, 50, 101)
    longitude = np.linspace(0, 360, 73)
    altitude = np.linspace(300, 500, 41)
    data = xr.Dataset({'uts': (('time'), index)},
                      coords={'time': index, 'latitude': latitude,
                              'longitude': longitude, 'altitude': altitude})

    slt = np.zeros([len(uts), len(longitude)])
    for i in range(len(uts)):
        for j in range(len(longitude)):
            slt[i, j] = np.mod(uts[i]/3600.0 + longitude[j]/15.0, 24.0)
    data['slt'] = (('time', 'longtiude'), slt)

    # Fake 3D data consisting of ones everywhere
    dummy1 = np.ones([len(uts), len(latitude), len(longitude)])
    data['dummy1'] = (('time', 'latitude', 'longitude'), dummy1)

    # Fake 4D data consisting of ones everywhere
    dummy2 = np.ones([len(uts), len(latitude), len(longitude), len(altitude)])
    data['dummy2'] = (('time', 'latitude', 'longitude', 'altitude'), dummy2)

    meta = pysat.Meta()

    return data, meta


def list_files(tag=None, sat_id=None, data_path=None, format_str=None,
               file_date_range=None):
    """Produce a fake list of files spanning a year

    Parameters
    ----------
    tag : (str)
        pysat instrument tag (default=None)
    sat_id : (str)
        pysat satellite ID tag (default=None)
    data_path : (str)
        pysat data path (default=None)
    format_str : (str)
        file format string (default=None)
    file_date_range : (pds.date_range)
        File date range (default=None)

    Returns
    -------
    Series of filenames indexed by file time

    """

    # Determine the appropriate date range for the fake files
    if file_date_range is None:
        start = _test_dates[''][''] - pds.DateOffset(years=1)
        stop = _test_dates[''][''] + pds.DateOffset(years=2) \
            - pds.DateOffset(days=1)
        file_date_range = pds.date_range(start, stop)

    index = file_date_range

    # Create the list of fake filenames
    names = [data_path + date.strftime('%Y-%m-%d') + '.nofile'
             for date in index]

    return pysat.Series(names, index=index)


def download(date_array, tag, sat_id, data_path=None,
             user=None, password=None):
    """ Download routine, not used since files are created locally"""
    pass
