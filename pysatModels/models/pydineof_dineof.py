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

import logging

import pysat

logger = logging.getLogger(__name__)

# the platform and name strings associated with this instrument
# need to be defined at the top level
platform = 'pydineof'
name = 'dineof'

# dictionary of data 'tags' and corresponding description
tags = {'*': 'Any pyDINEOF model export data set.'}
inst_ids = {'': ['*']}

# Define good days to download data for when pysat undergoes testing.
# format is outer dictionary has inst_id as the key
# each inst_id has a dictionary of test dates keyed by tag string
_test_dates = {'': {'': None}}

# Set to False to specify using xarray (not using pandas)
pandas_format = False


def init(self):
    """Initializes the Instrument object with instrument specific values.

    Runs once upon instantiation.

    """

    logger.info("DINEOF export models are produced by pyDINEOF.")
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
    return


def load(fnames, tag=None, inst_id=None):
    """Loads pyDINEOF data into xarray.

    This routine is called as needed by pysat. It is not intended
    for direct user interaction.

    Parameters
    ----------
    fnames : array-like
        iterable of filename strings, full path, to data files to be loaded.
        This input is nominally provided by pysat itself.
    tag : string ('')
        tag name used to identify particular data set to be loaded.
        This input is nominally provided by pysat itself. While
        tag defaults to None here, pysat provides '' as the default
        tag unless specified by user at Instrument instantiation.
    inst_id : string ('')
        Instrument ID used to identify particular data set to be loaded.
        This input is nominally provided by pysat itself.

    Returns
    -------
    data, metadata
        Data and Metadata are formatted for pysat. Data is an xarray
        DataSet while metadata is a pysat.Meta instance.

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

    # netCDF4 files were produced by xarray
    # returning an xarray.Dataset
    return pysat.utils.load_netcdf4(fnames, epoch_name='time',
                                    pandas_format=False)


def list_files(tag=None, inst_id=None, data_path=None, format_str=None):
    """Produce a list of files corresponding to pyDINEOF models.

    This routine is invoked by pysat and is not intended for direct
    use by the end user. Arguments are provided by pysat.

    Parameters
    ----------
    tag : string ('')
        tag name used to identify particular data set to be loaded.
        This input is nominally provided by pysat itself.
    inst_id : string ('')
        Instrument ID used to identify particular data set to be loaded.
        This input is nominally provided by pysat itself.
    data_path : string (None)
        Full path to directory containing files to be loaded. This
        is provided by pysat. The user may specify their own data path
        at Instrument instantiation and it will appear here.
    format_str : string (None)
        String template used to parse the datasets filenames. If a user
        supplies a template string at Instrument instantiation
        then it will appear here, otherwise defaults to None.

    Returns
    -------
    pandas.Series
        Series of filename strings, including the path, indexed by datetime.

    Examples
    --------
    ::
        If a filename is dineof_2009-01-01.nc then the template
        is 'dineof_{year:04d}-{month:02d}-{day:02d}.nc'

    Note
    ----
    The returned Series should not have any duplicate datetimes. If there are
    multiple versions of a file the most recent version should be kept and the
    rest discarded. This routine uses the pysat.Files.from_os constructor, thus
    the returned files are up to pysat specifications.

    Multiple data levels may be supported via the 'tag' input string.
    Multiple instruments via the inst_id string.


    """

    if format_str is None:
        # default string
        format_str = 'dineof_{year:04d}-{month:02d}-{day:02d}.nc'

    # use a pysat provided function to grab list of files from the
    # local file system that match the format defined above
    return pysat.Files.from_os(data_path=data_path, format_str=format_str)


def download(date_array, tag, inst_id, data_path=None, **kwargs):
    """Downloads are not supported.

    This routine is invoked by pysat and is not intended for direct use by the
    end user.

    Parameters
    ----------
    date_array : array-like
        list of datetimes to download data for. The sequence of dates need not
        be contiguous.
    tag : string ('')
        Tag identifier used for particular dataset. This input is provided by
        pysat.
    inst_id : string  ('')
        Instrument ID string identifier used for particular dataset. This input
        is provided by pysat.
    data_path : string (None)
        Path to directory to download data to.
    **kwargs : dict
        Additional keywords supplied by user when invoking the download
        routine attached to a pysat.Instrument object are passed to this
        routine via kwargs.

    """

    logger.warning('Downloads are not supported by pydineof_dineof.')

    return
