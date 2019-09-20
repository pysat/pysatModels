from nose.tools import raises

import pysat
import pysatModelUtils as mu


class TestBasics():
    def setup(self):
        """Runs before every method to create a clean testing setup."""
        self.testInst = pysat.Instrument('pysat', 'testing',
                                         clean_level='clean')
        self.start = pysat.datetime(2009, 1, 1)
        self.stop = pysat.datetime(2009, 1, 1)

    def teardown(self):
        """Runs after every method to clean up previous testing."""
        del self.testInst, self.start, self.stop

    @raises(ValueError)
    def test_collect_inst_model_pairs_wo_date(self):
        """Try to run without start or stop dates"""
        mu.collect_inst_model_pairs(inst=self.testInst)

    @raises(ValueError)
    def test_collect_inst_model_pairs_wo_inst(self):
        """Try to run without an instrument"""
        mu.collect_inst_model_pairs(start=self.start, stop=self.stop)

    @raises(ValueError)
    def test_collect_inst_model_pairs_wo_model(self):
        """Try to run without a model"""
        mu.collect_inst_model_pairs(start=self.start, stop=self.stop,
                                    inst=self.testInst)
