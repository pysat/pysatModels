"""Unit and Integration Tests for each instrument module.

Note
----
Imports test methods from pysat.tests.instrument_test_class

"""

import datetime as dt
import os
import shutil
import sys
import tempfile

# Import the test classes from pysat
import pysat
from pysat.tests.classes import cls_instrument_library as clslib

import pysatModels

# Retrieve the lists of Model instruments and testing methods
instruments = clslib.InstLibTests.initialize_test_package(
    clslib.InstLibTests, inst_loc=pysatModels.models)


class TestModels(clslib.InstLibTests):
    """Main class for instrument tests.

    Note
    ----
    All standard tests, setup, and teardown inherited from the core pysat
    instrument test class.

    """


class TestSAMILoadMultipleDays(object):
    """Test SAMI load code for multiple days."""

    def setup_class(self):
        """Initialize the testing setup once before all tests are run."""
        # Make sure to use a temporary directory so that the user setup is not
        # altered

        # TODO(#100): remove if-statement when it is always triggered
        tkwargs = {}
        if sys.version_info.major >= 3 and sys.version_info.minor >= 10:
            tkwargs = {"ignore_cleanup_errors": True}
        self.tempdir = tempfile.TemporaryDirectory(**tkwargs)
        self.saved_path = pysat.params['data_dirs']
        pysat.params.data['data_dirs'] = [self.tempdir.name]

        # Assign the location of the model Instrument sub-modules
        self.inst_loc = pysatModels.models
        return

    def teardown_class(self):
        """Clean up downloaded files and parameters from tests."""

        pysat.params.data['data_dirs'] = self.saved_path

        # Remove the temporary directory
        # TODO(#100): Remove try/except when Python 3.10 is the lowest version
        try:
            self.tempdir.cleanup()
        except Exception:
            pass

        del self.inst_loc, self.saved_path, self.tempdir
        return

    def test_multidayload(self):
        """Test loading multiple SAMI days."""

        inst = pysat.Instrument(inst_module=self.inst_loc.sami2py_sami2,
                                tag='test')
        test_date = self.inst_loc.sami2py_sami2._test_dates['']['test']
        inst.download(test_date, test_date)
        dl_file = os.path.join(inst.files.data_path, inst.files.files.values[0])

        # Create a second file, increment day.
        tmpl = 'sami2py_output_{year:04d}-{month:02d}-{day:02d}.nc'
        date = test_date + dt.timedelta(days=1)
        cp_file = os.path.join(inst.files.data_path,
                               tmpl.format(year=date.year, month=date.month,
                                           day=date.day))
        shutil.copy(dl_file, cp_file)

        # Update file list
        inst.files.refresh()

        # Load both files
        end_date = date + dt.timedelta(days=1)
        inst.load(date=test_date, end_date=end_date)
        data_two_days = inst.data
        meta_two_days = inst.meta

        # Load original download file
        inst.load(date=test_date)
        data_one_days = inst.data
        meta_one_days = inst.meta

        # Confirm longer date range is longer data
        assert len(data_two_days['slt']) == 2 * len(data_one_days['slt'])

        # Confirm metadata came out ok
        assert meta_two_days == meta_one_days

        # Confirm time index information
        assert data_two_days.indexes['time'].is_monotonic_increasing
        assert data_two_days.indexes['time'][0] >= test_date
        assert data_two_days.indexes['time'][0] < date
        assert data_two_days.indexes['time'][-1] > date
        assert data_two_days.indexes['time'][-1] < end_date

        return
