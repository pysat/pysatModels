"""Unit and Integration Tests for each instrument module.

Note
----
Imports test methods from pysat.tests.instrument_test_class

"""

import datetime as dt
import os
import pytest
import shutil

import pysat

import pysatModels

# Import the test classes from pysat
from pysat.tests.classes import cls_instrument_library as clslib


# Tell the standard tests which instruments to run each test on.
# Need to return instrument list for custom tests.
instruments = clslib.InstLibTests.initialize_test_package(
    clslib.InstLibTests, inst_loc=pysatModels.models)


class TestModels(clslib.InstLibTests):
    """Main class for instrument tests.

    Note
    ----
    All standard tests, setup, and teardown inherited from the core pysat
    instrument test class.

    """

    @pytest.mark.second
    def test_sami_multidayload(self):
        """Test loading multiple SAMI days."""

        inst = pysat.Instrument(inst_module=self.inst_loc.sami2py_sami2,
                                tag='test')
        test_date = self.inst_loc.sami2py_sami2._test_dates['']['test']
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
        assert data_two_days.indexes['time'].is_monotonic
        assert data_two_days.indexes['time'][0] >= test_date
        assert data_two_days.indexes['time'][0] < date
        assert data_two_days.indexes['time'][-1] > date
        assert data_two_days.indexes['time'][-1] < end_date

        return
