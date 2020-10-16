# -*- coding: utf-8 -*-
"""
Supports loading data from files generated using TIEGCM
(Thermosphere Ionosphere Electrodynamics General Circulation Model) model.
TIEGCM file is a netCDF file with multiple dimensions for some variables.

Properties
----------
platform
    'ucar'
name
    'tiegcm'
tag
    None supported
sat_id
    None supported

"""

import datetime as dt
import functools
import logging
import warnings

import xarray as xr
import pysat
from pysat.instruments.methods import general as mm_gen

logger = logging.getLogger('pysat')


platform = 'ucar'
name = 'tiegcm'

# dictionary of data 'tags' and corresponding description
tags = {'': 'UCAR TIE-GCM file'}
# dictionary of satellite IDs, list of corresponding tags for each sat_ids
sat_ids = {'': ['']}
# good day to download test data for. Downloads aren't currently supported!
# format is outer dictionary has sat_id as the key
# each sat_id has a dictionary of test dates keyed by tag string
_test_dates = {'': {'': dt.datetime(2019, 1, 1)}}
_test_download = {'': {'': False}}

# specify using xarray (not using pandas)
pandas_format = False

fname = 'tiegcm_icon_merg2.0_totTgcm.s_{day:03d}_{year:4d}.nc'
supported_tags = {'': {'': fname}}
list_files = functools.partial(mm_gen.list_files,
                               supported_tags=supported_tags)

ack = " ".join(["References and information about TIEGCM are available at",
                "https://www.hao.ucar.edu/modeling/tgcm/index.php"])
refs = [" ".join(("Dickinson, R. E., E. C. Ridley and R. G. Roble, A",
                  "three-dimensional general circulation model of the",
                  "thermosphere, J. Geophys. Res., 86, 1499-1512, 1981.")),
        " ".join(("Dickinson, R. E., E. C. Ridley and R. G. Roble,",
                  "Thermospheric general circulation with coupled dynamics and",
                  "composition, J. Atmos. Sci., 41, 205-219, 1984.")),
        " ".join(("Roble, R. G., and E. C. Ridley, An auroral model for the",
                  "NCAR thermospheric general circulation model (TGCM),",
                  "Annales Geophys., 5A, 369-382, 1987.")),
        " ".join(("Roble, R. G., E. C. Ridley and R. E. Dickinson, On the",
                  "global mean structure of the thermosphere, J. Geophys.",
                  "Res., 92, 8745-8758, 1987.")),
        " ".join(("Roble, R. G., E. C. Ridley, A. D. Richmond and R. E.",
                  "Dickinson, A coupled thermosphere/ionosphere general",
                  "circulation model, Geophys. Res. Lett., 15, 1325-1328,",
                  "1988.")),
        " ".join(("Richmond, A. D., E. C. Ridley and R. G. Roble, A",
                  "Thermosphere/Ionosphere General Circulation Model with",
                  "coupled electrodynamics, Geophys. Res. Lett., 19, 601-604,",
                  "1992.")),
        " ".join(("Roble, R. G., and E. C. Ridley, A",
                  "thermosphere-ionosphere-mesosphere-electrodynamics general",
                  "circulation model (TIME-GCM): equinox solar cycle minimum",
                  "simulations (30-500 km), Geophys. Res. Lett., 21, 417-420,",
                  "1994.")),
        " ".join(("Roble, R. G., Energetics of the mesosphere and",
                  "thermosphere, AGU, Geophysical Monographs, eds. R. M.",
                  "Johnson and T. L. Killeen, 87, 1-22, 1995.")),
        " ".join(("Wang, W., M. Wiltberger, A. G. Burns, S. Solomon, T. L.",
                  "Killeen, N. Maruyama, and J. Lyon, Initial results from the",
                  "CISM coupled magnetosphere-ionosphere-thermosphere (CMIT)",
                  "model: thermosphere ionosphere responses, J. Atmos.",
                  "Sol.-Terr. Phys., 66, 1425-1442,",
                  "doi:10.1016/j.jastp.2004.04.008, 2004.")),
        " ".join(("Solomon, S. C., and L. Y. Qian, Solar extreme-ultraviolet",
                  "irradiance for general circulation models, J. Geophys.",
                  "Res., 110, A10306, doi:10.1029/2005JA011160, 2005.")),
        " ".join(("Qian, L., A. G. Burns, B. A. Emery, B. Foster, G. Lu, A.",
                  "Maute, A. D. Richmond, R. G. Roble, S. C. Solomon, and W.",
                  "Wangm, The NCAR TIE-GCM: A community model of the coupled",
                  "thermosphere/ionosphere system, in Modeling the",
                  "Ionosphere-Thermosphere System, AGU Geophysical Monograph",
                  "Series, 2014."))]


def init(self):
    """Initializes the Instrument object with instrument specific values.

    Runs once upon instantiation.

    Parameters
    ----------
    self : pysat.Instrument
        This object

    """

    self.meta.acknowledgements = ack
    self.meta.references = "\n".join((refs))
    logger.info(self.meta.acknowledgements)
    return


def load(fnames, tag=None, sat_id=None, **kwargs):
    """Loads TIEGCM data using xarray.

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


    data, meta = pysat.utils.load_netcdf4(fnames, pandas_format=False)

    # move misc parameters from xarray to the Instrument object via Meta
    # doing this after the meta ensures all metadata is still kept
    # even for moved variables
    meta.p0 = data['p0']
    meta.p0_model = data['p0_model']
    meta.grav = data['grav']
    meta.mag = data['mag']
    meta.timestep = data['timestep']
    # remove these variables from xarray
    data = data.drop(['p0', 'p0_model', 'grav', 'mag', 'timestep'])

    return data, meta


def download(date_array, tag, sat_id, data_path=None, user=None, password=None,
             **kwargs):
    """Placeholder for UCAR TIEGCM downloads. Doesn't do anything.

    Parameters
    ----------
    date_array : array-like
        list of datetimes to download data for. The sequence of dates need not
        be contiguous.
    tag : string
        Tag identifier used for particular dataset. This input is provided by
        pysat. (default='')
    sat_id : string
        Satellite ID string identifier used for particular dataset. This input
        is provided by pysat. (default='')
    data_path : string
        Path to directory to download data to. (default=None)
    user : string
        User string input used for download. Provided by user and passed via
        pysat. If an account is required for dowloads this routine here must
        error if user not supplied. (default=None)
    password : string
        Password for data download. (default=None)
    **kwargs : dict
        Additional keywords supplied by user when invoking the download
        routine attached to a pysat.Instrument object are passed to this
        routine via kwargs.

    Note
    ----
    This routine is invoked by pysat and is not intended for direct use by
    the end user.

    """

    warnings.warn('Not implemented in this version.')
    return
