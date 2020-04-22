# Author: Angeline Burrell, NRL, 2019

from __future__ import absolute_import, unicode_literals

import numpy as np
import pytest

import pysat
from pysat.instruments import pysat_testmodel

import pysatModels.utils.extract as extract


@pytest.mark.skip("input requires a regular grid for the model")
class TestUtilsExtractObsViewModel:
    """ Unit tests for utils.extract.sat_view_through_model"""
    def setup(self):
        """ Runs before every method to create a clean testing setup."""
        self.args = [pysat.Instrument(platform=str('pysat'),
                                      name=str('testing'),
                                      sat_id='12', clean_level='clean'),
                     pysat.Instrument(platform=str('pysat'),
                                      name=str('testing'),
                                      sat_id='6', clean_level='clean'),
                     ['latitude', 'longitude'], ['dummy1', 'dummy2']]

    def teardown(self):
        del self.args

    def test_str_coords(self):
        """ Test string coordinate input"""
        self.args[2] = self.args[2][0]
        extract.instrument_view_through_model(*self.args)

        for label in self.args[3]:
            assert "model_{:s}".format(label) in self.args[1].data.columns

    def test_str_obs(self):
        """ Test string model observation input"""
        self.args[3] = self.args[3][0]
        extract.instrument_view_through_model(*self.args)

        assert "model_{:s}".format(self.args[3]) in self.args[1].data.columns


class TestUtilsExtractModObs:
    """ Unit tests for utils.extract.extract_modelled_observations """
    def setup(self):
        """Runs before every method to create a clean testing setup."""
        self.inst = pysat.Instrument(platform=str('pysat'),
                                     name=str('testing'), sat_id='10',
                                     clean_level='clean')
        self.model = pysat.Instrument(inst_module=pysat_testmodel)
        self.inst.load(yr=2009, doy=1)
        self.model.load(yr=2009, doy=1)
        self.input_args = [self.inst, self.model.data,
                           ["longitude", "latitude", "slt"],
                           ["longitude", "latitude", "slt"],
                           "time", "time", ["deg", "deg", "h"]]

    def teardown(self):
        """Runs after every method to clean up previous testing."""
        del self.inst, self.model, self.input_args

    @pytest.mark.parametrize("bad_index,bad_input,err_msg",
                             [(2, [], "Must provide instrument location"),
                              (2, ['glon', 'latitude', 'slt'],
                               "Unknown instrument location index"),
                              (3, [], "Must provide the same number"),
                              (6, [], "Must provide units for each "),
                              (4, "naname", "unknown model name for datetime"),
                              (5, "naname", "Unknown model time coordinate")])
    def test_bad_arg_input(self, bad_index, bad_input, err_msg):
        """ Test for expected failure with bad input arguments
        """
        self.input_args[bad_index] = bad_input

        with pytest.raises(ValueError) as verr:
            extract.extract_modelled_observations(*self.input_args)

        assert str(verr.value.args[0]).find(err_msg) >= 0

    @pytest.mark.parametrize("bad_key,bad_val,err_msg",
                             [("sel_name", ["slt"],
                               "No model data keys to interpolate"),
                              ("method", "not_a_method",
                               "interpn only understands the methods"),
                              ("model_label", 1, "Unknown format code "),
                              ("time_method", "fun", "unknown time method"),
                              ("pair_method", "fun",
                               "unknown pairing method")])
    def test_bad_kwarg_input(self, bad_key, bad_val, err_msg):
        """ Test for expected failure with bad kwarg input """
        kwargs = {bad_key: bad_val}

        with pytest.raises(Exception) as err:
            extract.extract_modelled_observations(*self.input_args, **kwargs)

        assert str(err.value.args[0]).find(err_msg) >= 0

    @pytest.mark.parametrize("sel_val", [["dummy1", "dummy2"], ["dummy1"]])
    def test_good_sel_name(self, sel_val):
        """ Test for success with different good selection name inputs"""
        out_keys = extract.extract_modelled_observations(*self.input_args,
                                                         sel_name=sel_val)
        for label in sel_val:
            assert "model_{:s}".format(label) in out_keys
        assert len(out_keys) == len(np.asarray(sel_val))

    # TODO: Add test for out-of-bounds data
    # TODO: Add tests for model data already in instrument

    def test_success(self):
        """ Test the extraction success"""
        out_keys = extract.extract_modelled_observations(*self.input_args)

        for label in self.model.data.data_vars.keys():
            if label not in self.input_args[3]:
                assert "model_{:s}".format(label) in out_keys
