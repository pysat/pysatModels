import pytest

import pysatModels

import pysat.tests.test_instruments
from pysat.tests.test_instruments import generate_instrument_list
from pysat.tests.test_instruments import TestInstrumentsAll
from pysat.tests.test_instruments import TestInstrumentsDownload
from pysat.tests.test_instruments import TestInstrumentsNoDownload

instruments = generate_instrument_list(pysatModels.models.__all__,
                                       package='pysatModels.models')

InstClasses = [TestInstrumentsAll, TestInstrumentsDownload,
               TestInstrumentsNoDownload]

for InstClass in InstClasses:
    InstClass.__test__ = False
    method_list = [func for func in dir(InstClass)
                   if callable(getattr(InstClass, func))]
    # Search tests for iteration via pytestmark, update instrument list
    for method in method_list:
        if hasattr(getattr(InstClass, method), 'pytestmark'):
            Nargs =  len(getattr(InstClass, method).pytestmark)
            names = [getattr(InstClass, method).pytestmark[j].name
                     for j in range(0, Nargs)]
            getattr(InstClass, method).pytestmark.clear()
            if InstClass == TestInstrumentsAll:
                mark = pytest.mark.parametrize("name", instruments['names'])
                getattr(InstClass, method).pytestmark.append(mark)
            elif InstClass == TestInstrumentsDownload:
                mark = pytest.mark.parametrize("inst", instruments['download'])
                getattr(InstClass, method).pytestmark.append(mark)
            elif InstClass == TestInstrumentsNoDownload:
                mark = pytest.mark.parametrize("inst",
                                               instruments['no_download'])
                getattr(InstClass, method).pytestmark.append(mark)


class TestModels(TestInstrumentsAll):

    __test__ = True

    def setup(self):
        """Runs before every method to create a clean testing setup."""
        self.package = 'pysatModels.models'
        pass

    def teardown(self):
        """Runs after every method to clean up previous testing."""
        pass


class TestModelsDownload(TestInstrumentsDownload):

    __test__ = True

    def setup(self):
        """Runs before every method to create a clean testing setup."""
        self.package = 'pysatModels.models'
        pass

    def teardown(self):
        """Runs after every method to clean up previous testing."""
        pass


class TestModelsNoDownload(TestInstrumentsNoDownload):

    __test__ = True

    def setup(self):
        """Runs before every method to create a clean testing setup."""
        self.package = 'pysatModels.models'
        pass

    def teardown(self):
        """Runs after every method to clean up previous testing."""
        pass
