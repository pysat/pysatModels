# -*- coding: utf-8 -*-
"""General functions for model instruments."""

import os
import requests
import warnings

import pysatModels

clean_warn = {clean_level: [('logger', 'INFO', 'Cleaning not supported',
                             clean_level)]
              for clean_level in ['clean', 'dusty', 'dirty']}


def clean(inst):
    """Raise a low-level log message about lack of cleaning."""

    pysatModels.logger.info('Cleaning not supported for {:} {:}'.format(
        inst.platform, inst.name))

    return


def download_test_data(remote_url, remote_file, data_path, test_date=None,
                       format_str=None):
    """Download test data from an online repository.

    Parameters
    ----------
    remote_url : str
        URL of the target repository, including the path to the test file
    remote_file : str
        Remote file name
    data_path : str
        Path to directory where local file will be stored
    test_date : dt.datetime or NoneType
        Datetime for which the test file will be assigned, does not need to
        correspond to the test model run time. Only used if `format_str` is
        also provided. (default=None)
    format_str : str or NoneType
        Format string to construct a pysat-compatible filename or None to
        not change the filename (default=None)

    Note
    ----
    This routine is invoked by pysat and is not intended for direct use by
    the end user.

    The test object generates the datetime requested by the user, which may not
    match the date of the model run.

    """

    # Use pysat-compatible name
    if format_str is None or test_date is None:
        # Strip potential web info off the end of the remote filename
        local_fname = os.path.join(data_path, remote_file.split('?')[0])
    else:
        local_fname = os.path.join(
            data_path, format_str.format(year=test_date.year,
                                         month=test_date.month,
                                         day=test_date.day,
                                         hour=test_date.hour,
                                         minute=test_date.minute,
                                         second=test_date.second))

    # Ensure there are no double backslashes
    remote_path = '/'.join((remote_url.strip('/'), remote_file))

    # Get the remote file
    with requests.get(remote_path) as req:
        if req.status_code != 404:
            with open(local_fname, 'wb') as open_f:
                open_f.write(req.content)
        else:
            warnings.warn(''.join(('Unable to find remote file: ',
                                   remote_path)))

    return
