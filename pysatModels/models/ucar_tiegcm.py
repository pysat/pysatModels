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
inst_id
    None supported

"""

import datetime as dt
import functools
import warnings

import pysat

logger = pysat.logger

# ----------------------------------------------------------------------------
# Instrument attributes

platform = 'ucar'
name = 'tiegcm'
tags = {'': 'UCAR TIE-GCM file'}
inst_ids = {'': ['']}

# specify using xarray (not using pandas)
pandas_format = False

# ----------------------------------------------------------------------------
# Instrument test attributes

_test_dates = {'': {'': dt.datetime(2019, 1, 1)}}
_test_download = {'': {'': False}}

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
    ack = " ".join(["References and information about TIEGCM are available at",
                    "https://www.hao.ucar.edu/modeling/tgcm/index.php"])
    refs = [" ".join(("Dickinson, R. E., E. C. Ridley and R. G. Roble, A",
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

    self.acknowledgements = ack
    self.references = "\n".join((refs))
    logger.info(self.acknowledgements)
    return


# Required method
def clean(self):
    """Method to return UCER/TIEGCM data cleaned to the specified level

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

    logger.info('Cleaning not supported or needed for TIEGCM')
    return

# ----------------------------------------------------------------------------
# Instrument functions
#
# Use local and default pysat methods


# Set the list_files routine
fname = 'tiegcm_icon_merg2.0_totTgcm.s_{day:03d}_{year:4d}.nc'
supported_tags = {'': {'': fname}}
list_files = functools.partial(pysat.instruments.methods.general.list_files,
                               supported_tags=supported_tags)


def load(fnames, tag=None, inst_id=None, **kwargs):
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

    Examples
    --------
    ::

        inst = pysat.Instrument('ucar', 'tiegcm')
        inst.load(2019, 1)

    """

    data, meta = pysat.utils.load_netcdf4(fnames, pandas_format=False)

    # Move misc parameters from xarray to the Instrument object via Meta
    # doing this after the meta ensures all metadata is still kept
    # even for moved variables
    meta.p0 = data['p0']
    meta.p0_model = data['p0_model']
    meta.grav = data['grav']
    meta.mag = data['mag']
    meta.timestep = data['timestep']

    # Remove these variables from xarray
    data = data.drop(['p0', 'p0_model', 'grav', 'mag', 'timestep'])

    return data, meta


def download(date_array, tag, inst_id, data_path=None, **kwargs):
    """Placeholder for UCAR TIEGCM downloads. Doesn't do anything.

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

    """

    warnings.warn('Not implemented in this version.')
    return
