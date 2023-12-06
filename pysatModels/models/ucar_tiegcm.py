# -*- coding: utf-8 -*-
"""Support loading data from files generated using TIEGCM model.

TIEGCM (Thermosphere Ionosphere Electrodynamics General Circulation Model)
file is a netCDF file with multiple dimensions for some variables. Support
is included for TIEGCM files from the University Corporation for Atmospheric
Research (UCAR) as well as TIEGCM files produced as part of the
NASA Ionospheric CONnections (ICON) Explorer Mission. TIEGCM files produced for
ICON are generated using Hough Mode Extensions generated from ICON measurements.
Download support is only available for ICON TIEGCM files.

Properties
----------
platform
    'ucar'
name
    'tiegcm'
tag
    ''
inst_id
    ''

"""

import datetime as dt
import functools
import os
import pandas as pds
import tempfile
import warnings
import zipfile

import pysat

from pysatModels.models.methods import general
from pysatNASA.instruments.methods import cdaweb as cdw

logger = pysat.logger

# ----------------------------------------------------------------------------
# Instrument attributes

platform = 'ucar'
name = 'tiegcm'
tags = {'': 'UCAR TIE-GCM file', 'icon': 'UCAR TIE-GCM file produced for ICON.'}
inst_ids = {'': list(tags.keys())}

# Specify using xarray (not using pandas)
pandas_format = False

# ----------------------------------------------------------------------------
# Instrument test attributes

_test_dates = {'': {'': dt.datetime(2019, 1, 1),
                    'icon': dt.datetime(2020, 1, 10)}}
_test_download = {'': {'': False, 'icon': True}}
_clean_warn = {'': {tag: general.clean_warn for tag in tags.keys()}}

# ----------------------------------------------------------------------------
# Instrument methods

clean = general.clean


def init(self):
    """Initialize the Instrument object with instrument specific values."""

    self.acknowledgements = "".join([
        "References and information about TIEGCM are available at ",
        "https://www.hao.ucar.edu/modeling/tgcm/index.php"])

    if self.tag == '':
        refs = [
            " ".join(("Dickinson, R. E., E. C. Ridley and R. G. Roble, A",
                      "three-dimensional general circulation model of the",
                      "thermosphere, J. Geophys. Res., 86, 1499-1512, 1981.")),
            " ".join(("Dickinson, R. E., E. C. Ridley and R. G. Roble,",
                      "Thermospheric general circulation with coupled",
                      "dynamics and composition, J. Atmos. Sci., 41, 205-219,",
                      "1984.")),
            " ".join(("Roble, R. G., and E. C. Ridley, An auroral model for",
                      "the NCAR thermospheric general circulation model",
                      "(TGCM), Annales Geophys., 5A, 369-382, 1987.")),
            " ".join(("Roble, R. G., E. C. Ridley and R. E. Dickinson, On the",
                      "global mean structure of the thermosphere, J. Geophys.",
                      "Res., 92, 8745-8758, 1987.")),
            " ".join(("Roble, R. G., E. C. Ridley, A. D. Richmond and R. E.",
                      "Dickinson, A coupled thermosphere/ionosphere general",
                      "circulation model, Geophys. Res. Lett., 15, 1325-1328,",
                      "1988.")),
            " ".join(("Richmond, A. D., E. C. Ridley and R. G. Roble, A",
                      "Thermosphere/Ionosphere General Circulation Model with",
                      "coupled electrodynamics, Geophys. Res. Lett., 19,",
                      "601-604, 1992.")),
            " ".join(("Roble, R. G., and E. C. Ridley, A",
                      "thermosphere-ionosphere-mesosphere-electrodynamics",
                      "general circulation model (TIME-GCM): equinox solar",
                      "cycle minimum simulations (30-500 km), Geophys. Res.",
                      "Lett., 21, 417-420, 1994.")),
            " ".join(("Roble, R. G., Energetics of the mesosphere and",
                      "thermosphere, AGU, Geophysical Monographs, eds. R. M.",
                      "Johnson and T. L. Killeen, 87, 1-22, 1995.")),
            " ".join(("Wang, W., M. Wiltberger, A. G. Burns, S. Solomon, T. L.",
                      "Killeen, N. Maruyama, and J. Lyon, Initial results",
                      "from the CISM coupled magnetosphere-ionosphere-",
                      "thermosphere (CMIT) model: thermosphere ionosphere ",
                      "responses, J. Atmos. Sol.-Terr. Phys., 66, 1425-1442,",
                      "doi:10.1016/j.jastp.2004.04.008, 2004.")),
            " ".join(("Solomon, S. C., and L. Y. Qian, Solar",
                      "extreme-ultraviolet irradiance for general circulation",
                      "models, J. Geophys. Res., 110, A10306,",
                      "doi:10.1029/2005JA011160, 2005.")),
            " ".join(("Qian, L., A. G. Burns, B. A. Emery, B. Foster, G. Lu,",
                      "A. Maute, A. D. Richmond, R. G. Roble, S. C. Solomon,",
                      "and W. Wangm, The NCAR TIE-GCM: A community model of",
                      "the coupled thermosphere/ionosphere system, in",
                      "Modeling the Ionosphere-Thermosphere System, AGU",
                      "Geophysical Monograph Series, 2014."))]

    elif self.tag == 'icon':
        refs = [" ".join(["Maute, A. Thermosphere-Ionosphere-Electrodynamics",
                          "General Circulation Model for the Ionospheric",
                          "Connection Explorer: TIEGCM-ICON. Space Sci Rev",
                          "212, 523–551 (2017).",
                          "https://doi.org/10.1007/s11214-017-0330-3"])]

    self.references = "\n".join((refs))
    logger.info(self.acknowledgements)
    return


# ----------------------------------------------------------------------------
# Instrument functions
#
# Use local and default pysat methods

# Set the list_files routine. TIEGCM files downloaded from CDAWeb have an
# intermediate day directory. Accounted for with the leading '*' for icon.
supported_tags = {
    '': {'': 'tiegcm_icon_merg2.0_totTgcm.s_{day:03d}_{year:4d}.nc',
         'icon': os.path.join('*', ''.join(['ICON_L4-3_TIEGCM_{year:04d}-',
                                            '{month:02d}-{day:02d}_',
                                            'v{version:02d}r{revision:03d}',
                                            '.NC']))}}
list_files = functools.partial(pysat.instruments.methods.general.list_files,
                               supported_tags=supported_tags)


def load(fnames, tag='', inst_id='', **kwargs):
    """Load TIE-GCM data using xarray.

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

        inst = pysat.Instrument('ucar', 'tiegcm')
        inst.load(2019, 1)

    """

    data, meta = pysat.utils.io.load_netcdf(fnames, pandas_format=False,
                                            epoch_name='time',
                                            decode_times=True)

    # Move misc parameters from xarray to the Instrument object via Meta.
    # Doing this after `meta` created ensures all metadata is still kept
    # even for moved variables. TODO(#117) Move to `meta.header`.
    meta.p0 = data['p0']
    meta.p0_model = data['p0_model']
    meta.grav = data['grav']
    meta.mag = data['mag']
    meta.timestep = data['timestep']

    # Remove these variables from xarray
    data = data.drop(['p0', 'p0_model', 'grav', 'mag', 'timestep'])

    return data, meta


# Remote filenames are different than final model filenames.
# `download_tags` used by both `download` and `list_remote_files`.
download_tags = {
    '': {'icon': {'remote_dir': '/pub/data/icon/l4/tiegcm/{year:04d}/',
                  'fname': ''.join(['icon_l4-3_tiegcm_{year:04d}-{month:02d}-',
                                    '{day:02d}_v{version:02d}r{revision:03d}',
                                    '.zip'])}}}


def download(date_array, tag, inst_id, data_path=None, **kwargs):
    """Download UCAR TIE-GCM from NASA CDAWeb.

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

    """

    if tag == '':
        warnings.warn('Not implemented, currently no support for Globus.')
    elif tag == 'icon':
        # Set up temporary directory for zip files
        temp_dir = tempfile.TemporaryDirectory()

        # Download using NASA CDAWeb methods in pysatNASA
        cdw.download(date_array, temp_dir.name, tag=tag, inst_id=inst_id,
                     supported_tags=download_tags)

        # Get a list of files in `temp_dir`
        dl_files = pysat.Files.from_os(
            temp_dir.name, format_str=download_tags[inst_id][tag]['fname'])

        # Decompress files
        for dl_fname in dl_files.values:
            dl_fname = os.path.split(dl_fname)[1]
            with zipfile.ZipFile(os.path.join(temp_dir.name, dl_fname),
                                 'r') as open_zip:
                open_zip.extractall(data_path)

        # Cleanup temporary directory
        temp_dir.cleanup()

    return


def list_remote_files(tag='', inst_id='', start=None, stop=None):
    """Return a Pandas Series of every file for chosen remote data.

    Parameters
    ----------
    tag : str
        Denotes type of file to load. (default='')
    inst_id : str
        Specifies the instrument ID. (default='')
    start : dt.datetime or NoneType
        Starting time for file list. A None value will start with the first
        file found. (default=None)
    stop : dt.datetime or NoneType
        Ending time for the file list.  A None value will stop with the last
        file found. (default=None)

    Returns
    -------
    files : pds.Series
        A series containing the verified available files.

    """

    if tag == '':
        warnings.warn('Not implemented, currently no support for Globus.')
        files = pds.Series([], dtype='a')
    elif tag == 'icon':
        files = cdw.list_remote_files(tag=tag, inst_id=inst_id, start=start,
                                      stop=stop, supported_tags=download_tags)

    return files
