from __future__ import absolute_import, unicode_literals

from pytest import raises

import pysat

import pysatModelUtils.utils.match as match


class TestUtilsMatch():
    def setup(self):
        """Runs before every method to create a clean testing setup."""
        self.testInst = pysat.Instrument(platform=str('pysat'),
                                         name=str('testing'),
                                         clean_level='clean')
        self.start = pysat.datetime(2009, 1, 1)
        self.stop = pysat.datetime(2009, 1, 1)

    def teardown(self):
        """Runs after every method to clean up previous testing."""
        del self.testInst, self.start, self.stop

    def test_collect_inst_model_pairs_wo_date(self):
        """Try to run without start or stop dates"""
        with raises(ValueError) as verr:
            match.collect_inst_model_pairs(inst=self.testInst)

        assert verr.vaue.args[0].find('Must provide start and end time') >= 0

    def test_collect_inst_model_pairs_wo_inst(self):
        """Try to run without an instrument"""
        with raises(ValueError) as verr:
            match.collect_inst_model_pairs(start=self.start, stop=self.stop)

        assert verr.vaue.args[0].find('Must provide a pysat instrument') >= 0

    def test_collect_inst_model_pairs_wo_model(self):
        """Try to run without a model"""
        with raises(ValueError) as verr:
            match.collect_inst_model_pairs(start=self.start, stop=self.stop,
                                           inst=self.testInst)

        assert verr.vaue.args[0].find('Must provide model files') >= 0
