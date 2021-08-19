import pytest
import tempfile

import pysat
from pysat.tests.instrument_test_class import InstTestClass

import pysatModels

# Test code - remove methods
del InstTestClass.test_load

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

    def setup_class(self):
        """Runs once before the tests to initialize the testing setu
        """
        # Make sure to use a temporary directory so that the user setup is not
        # altered
        self.tempdir = tempfile.TemporaryDirectory()
        self.saved_path = pysat.params['data_dirs']
        pysat.params.data['data_dirs'] = [self.tempdir.name]

        # Assign the location of the model Instrument sub-modules
        self.inst_loc = pysatModels.models

    def teardown_class(self):
        """Runs after every method to clean up previous testing
        """
        pysat.params.data['data_dirs'] = self.saved_path
        self.tempdir.cleanup()
        del self.inst_loc, self.saved_path, self.tempdir
