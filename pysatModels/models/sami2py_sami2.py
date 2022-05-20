#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2022, pysat development team
# Full license can be found in License.md
# -----------------------------------------------------------------------------
"""Support loading data from files generated using the sami2py model.

sami2py file is a netCDF file with multiple dimensions for some variables.
The sami2py project is at https://github.com/sami2py/sami2py

Properties
----------
platform
    'sami2py'
name
    'sami2'
tag
    '', 'test'
inst_id
    ''

"""

import datetime as dt
import functools
from packaging import version as pack_version
import warnings
import xarray as xr

import pysat

from pysatModels.models.methods import general

logger = pysat.logger

# ----------------------------------------------------------------------------
# Instrument attributes

platform = 'sami2py'
name = 'sami2'
tags = {'': 'sami2py output file',
        'test': 'Standard output of sami2py for benchmarking'}
inst_ids = {'': [tag for tag in tags.keys()]}

# specify using xarray (not using pandas)
pandas_format = False

# ----------------------------------------------------------------------------
# Instrument test attributes

_test_dates = {'': {tag: dt.datetime(2019, 1, 1) for tag in tags.keys()}}
_test_download = {'': {'': False, 'test': True}}

# ----------------------------------------------------------------------------
# Instrument methods

clean = general.clean


def init(self):
    """Initialize the Instrument object with instrument specific values."""

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


# ----------------------------------------------------------------------------
# Instrument functions
#
# Use local and default pysat methods

# Set the list_files routine
sami_format = 'sami2py_output_{year:04d}-{month:02d}-{day:02d}.nc'
supported_tags = {'': {'': sami_format, 'test': sami_format}}
list_files = functools.partial(pysat.instruments.methods.general.list_files,
                               supported_tags=supported_tags)


def load(fnames, tag='', inst_id='', **kwargs):
    """Load sami2py data using xarray.

    This routine is called as needed by pysat. It is not intended
    for direct user interaction.

    Parameters
    ----------
    fnames : array-like
        Iterable of filename strings, full path, to data files to be loaded.
        This input is nominally provided by pysat itself.
    tag : str
        Tag name used to identify particular data set to be loaded.
        This input is nominally provided by pysat itself. (default='')
    inst_id : str
        Instrument ID used to identify particular data set to be loaded.
        This input is nominally provided by pysat itself. (default='')
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

        inst = pysat.Instrument('sami2py', 'sami2')
        inst.load(2019, 1)

    """

    # SAMI files only contain hours relative to the day in the time information.
    # Pull out dates from requested filenames and then use those dates
    # to construct a proper time.

    # SAMI template string defined as `sami_format` above.
    file_info = pysat.utils.files.parse_fixed_width_filenames(fnames,
                                                              sami_format)
    epochs = []
    for year, month, day in zip(file_info['year'], file_info['month'],
                                file_info['day']):
        epochs.append(dt.datetime(year, month, day))

    vstr = '3.0.2'  # TODO(#112) Remove support for backwards compatibility
    loaded_data = []
    loaded_meta = []
    for epoch, fname in zip(epochs, fnames):
        # Load data
        # TODO(#112) Remove backwards compatibility
        if pack_version.Version(pysat.__version__) < pack_version.Version(vstr):
            data, meta = pysat.utils.load_netcdf4([fname], pandas_format=False,
                                                  epoch_name='ut')
            data = data.rename({"ut": "time"})

            # Create datetimes from 'ut' variable
            data['time'] = [epoch
                            + dt.timedelta(seconds=int(val * 3600.0))
                            for val in data['time'].values]
        else:
            data, meta = pysat.utils.load_netcdf4([fname], pandas_format=False,
                                                  epoch_name='ut',
                                                  epoch_origin=epoch,
                                                  epoch_unit='h')

        # Store data/meta for each loop
        loaded_data.append(data)
        loaded_meta.append(meta)

        # Manually close link to file for peace of mind
        data.close()

    # Combine multiple datasets
    data = xr.combine_by_coords(loaded_data)

    # Combine metadata
    for loop_meta in loaded_meta[:-1]:
        meta.merge(loop_meta)

    return data, meta


def download(date_array, tag, inst_id, data_path):
    """Download sami2py data.

    Parameters
    ----------
    date_array : array-like
        List of datetimes to download data for. The sequence of dates need not
        be contiguous.
    tag : str
        Tag identifier used for particular dataset. This input is provided by
        pysat.
    inst_id : str
        Instrument ID string identifier used for particular dataset. This input
        is provided by pysat.
    data_path : str
        Path to directory to download data to.

    Note
    ----
    This routine is invoked by pysat and is not intended for direct use by
    the end user.

    The test object generates the datetime requested by the user, which may not
    match the date of the model run.

    Examples
    --------
    ::

        import datetime as dt
        import pysat

        inst = pysat.Instrument('sami2py', 'sami2', 'test')
        inst.download(start=dt.datetime(2020, 3, 8))

    """

    if tag == 'test':
        # Define the remote file data
        remote_url = ''.join(['https://github.com/sami2py/sami2py/blob',
                              '/main/sami2py/tests/test_data/'])
        fname = 'sami2py_output.nc?raw=true'  # Show raw data, not webpage ver.

        # Construct a format string
        format_str = supported_tags[inst_id][tag]

        # Download the remote test file
        general.download_test_data(remote_url, fname, data_path, date_array[0],
                                   format_str)
    else:
        warnings.warn('Downloads currently only supported for test files.')

    return
