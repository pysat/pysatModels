"""Unit and Integration Tests for each instrument module.

Note
----
Imports test methods from pysat.tests.instrument_test_class

"""

from packaging import version as pack_version
import pytest
import sys
import tempfile

import pysat
from pysat.tests.instrument_test_class import InstTestClass

import pysatModels

# Retrieve the lists of Model instruments and testing methods
instruments = pysat.utils.generate_instrument_list(inst_loc=pysatModels.models)
method_list = [func for func in dir(InstTestClass)
               if callable(getattr(InstTestClass, func))]

# Search tests for iteration via pytestmark, update instrument list
for method in method_list:
    if hasattr(getattr(InstTestClass, method), 'pytestmark'):
        # Get list of names of pytestmarks
        mark_name = [mod_mark.name for mod_mark
                     in getattr(InstTestClass, method).pytestmark]

        # Add instruments from your library
        if 'all_inst' in mark_name:
            mark = pytest.mark.parametrize("inst_name", instruments['names'])
            getattr(InstTestClass, method).pytestmark.append(mark)
        elif 'download' in mark_name:
            mark = pytest.mark.parametrize("inst_dict",
                                           instruments['download'])
            getattr(InstTestClass, method).pytestmark.append(mark)
        elif 'no_download' in mark_name:
            mark = pytest.mark.parametrize("inst_dict",
                                           instruments['no_download'])
            getattr(InstTestClass, method).pytestmark.append(mark)


class TestModels(InstTestClass):
    """Main class for instrument tests.

    Note
    ----
    Uses class level setup and teardown so that all tests use the same
    temporary directory. We do not want to geneate a new tempdir for each test,
    as the load tests need to be the same as the download tests.

    """

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

        # TODO(#100): Remove try/except when Python 3.10 is the lowest version
        try:
            self.tempdir.cleanup()
        except Exception:
            pass

        del self.inst_loc, self.saved_path, self.tempdir
        return


class TestSAMIPysatVersion(object):
    """Test SAMI load code for pysat version differences across v3.0.2."""

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
        self.saved_ver = pysat.__version__
        pysat.params.data['data_dirs'] = [self.tempdir.name]

        # Assign the location of the model Instrument sub-modules
        self.inst_loc = pysatModels.models
        return

    def teardown_class(self):
        """Clean up downloaded files and parameters from tests."""

        pysat.params.data['data_dirs'] = self.saved_path
        pysat.__version__ = self.saved_ver

        # Remove the temporary directory
        # TODO(#100): Remove try/except when Python 3.10 is the lowest version
        try:
            self.tempdir.cleanup()
        except Exception:
            pass

        del self.inst_loc, self.saved_path, self.tempdir, self.saved_ver
        return

    def test_load_failure(self):
        """Test for SAMI load failure when faking a different pysat version."""

        if (pack_version.Version(pysat.__version__)
                >= pack_version.Version('3.0.2')):
            # Define target error variable label
            label = 'ut'

            # Replace reported version with one before 3.0.2
            vlabel = '3.0.1'

            # Expected error
            error = KeyError
        else:
            # Define target error variable label
            label = 'epoch_origin'

            # Expected error
            error = TypeError

            # Replace reported version with one including 3.0.2
            vlabel = '3.0.2'

        # Update reported pysat version
        pysat.__version__ = vlabel

        with pytest.raises(error) as verr:
            inst = pysat.Instrument(inst_module=self.inst_loc.sami2py_sami2,
                                    tag='test')
            inst.download(self.inst_loc.sami2py_sami2._test_dates['']['test'],
                          self.inst_loc.sami2py_sami2._test_dates['']['test'])
            inst.load(date=self.inst_loc.sami2py_sami2._test_dates['']['test'])

        assert str(verr).find(label) >= 0

        return
