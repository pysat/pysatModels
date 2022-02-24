# Author: Angeline Burrell, NRL, 2019
"""Unit tests for `pysatModels.models.methods.general`."""

import datetime as dt
import logging
import os
import sys
import tempfile

import pysat

from pysatModels.models.methods import general


class TestMethodsGeneralLogging(object):
    """Unit tests for log messages raised by general methods."""

    def setup(self):
        """Set up the unit test environment."""
        self.ch = logging.StreamHandler()
        self.ch.setLevel(logging.INFO)
        self.model = pysat.Instrument("pysat", "testmodel")

        return

    def teardown(self):
        """Clean up the unit test environment."""

        del self.ch, self.model
        return

    def test_general_clean(self, caplog):
        """Test clean function provides expected logging information."""

        # Construct the expected log message
        msg = "Cleaning not supported for {:} {:}".format(self.model.platform,
                                                          self.model.name)

        # Run the clean function and retrieve logging messages
        with caplog.at_level(logging.INFO, logger='pysatModels'):
            general.clean(self.model)

        # Evaluate the logging output
        assert len(caplog.records) >= 1

        assert caplog.records[-1].message.find(msg) >= 0
        assert caplog.records[-1].levelname == "INFO"
        return


class TestMethodsGeneralDownload(object):
    """Unit tests for general methods handling downloads."""

    def setup(self):
        """Set up the unit test environment."""
        # TODO #100: remove if-statement when it is always triggered
        tkwargs = {}
        if sys.version_info.major >= 3 and sys.version_info.minor >= 10:
            tkwargs = {"ignore_cleanup_errors": True}
        self.tempdir = tempfile.TemporaryDirectory(**tkwargs)
        self.remote_url = ''.join(['https://github.com/pysat/pysatModels/',
                                   '/blob/main/pysatModels/tests/test_data/'])
        self.remote_file = 'dineof-2009-01-01.nc?raw=true'
        self.test_time = dt.datetime(2009, 1, 1)
        self.out_file = os.path.join(self.tempdir.name, 'dineof-2009-01-01.nc')

        return

    def teardown(self):
        """Clean up the unit test environment."""

        if os.path.isfile(self.out_file):
            os.remove(self.out_file)

        # Remove the temporary directory
        # TODO #100: Remove try/except when Python 3.10 is the lowest version
        try:
            self.tempdir.cleanup()
        except Exception:
            pass

        del self.tempdir, self.remote_url, self.remote_file, self.out_file
        del self.test_time
        return

    def test_download_test_data_reformat(self):
        """Test the download of remote test data with file renaming."""

        # Construct the missing input and expected output
        format_str = "dineof_{year:04d}{month:02d}{day:02d}.nc"
        self.out_file = os.path.join(self.tempdir.name, format_str.format(
            year=self.test_time.year, month=self.test_time.month,
            day=self.test_time.day))

        # Download the test file
        general.download_test_data(self.remote_url, self.remote_file,
                                   self.tempdir.name, self.test_time,
                                   format_str)

        # Ensure the expected file is present
        assert os.path.isfile(self.out_file), "did not download {:}".format(
            self.out_file)
        return

    def test_download_test_data_same(self):
        """Test the download of remote test data without file renaming."""

        # Download the test file
        general.download_test_data(self.remote_url, self.remote_file,
                                   self.tempdir.name)

        # Ensure the expected file is present
        assert os.path.isfile(self.out_file), "did not download {:}".format(
            self.out_file)
        return
