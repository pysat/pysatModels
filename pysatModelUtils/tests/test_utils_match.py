# Author: Angeline Burrell, NRL, 2019

from __future__ import absolute_import, unicode_literals

import datetime as dt
import pytest

import pysat

import pysatModelUtils as pysat_mu
import pysatModelUtils.utils.match as match


class TestUtilsMatchCollectInstModPairs():
    """ Unit tests for utils.match.collect_inst_model_pairs """
    def setup(self):
        """Runs before every method to create a clean testing setup."""
        self.input_args = {"start": dt.datetime(2009, 1, 1),
                           "stop": dt.datetime(2009, 1, 1),
                           "tinc": dt.timedelta(days=1),
                           "inst": pysat.Instrument(platform=str('pysat'),
                                                    name=str('testing'),
                                                    clean_level='clean')}
        self.required_kwargs = {"model_load_kwargs":
                                {"model_inst":
                                 pysat.Instrument(platform=str('pysat'),
                                                  name=str('testing'),
                                                  clean_level='clean')},
                                "inst_lon_name": "longitude",
                                "mod_lon_name": "longitude",
                                "inst_name": ["longitude", "latitude", "slt"],
                                "mod_name": ["longitude", "latitude", "slt"],
                                "mod_datetime_name": "Epoch",
                                "mod_time_name": "uts",
                                "mod_units": ["deg", "deg", "h"]}

    def teardown(self):
        """Runs after every method to clean up previous testing."""
        del self.input_args, self.required_kwargs

    @pytest.mark.parametrize("del_key,err_msg",
                             [("model_load_kwargs",
                               "must provide a pysat.Instrument object"),
                              ("inst_lon_name", "Need longitude name for inst"),
                              ("mod_lon_name", "Need longitude name for model"),
                              ("inst_name", "Must provide instrument location"),
                              ("mod_name", "Must provide the same number"),
                              ("mod_units", "Must provide units for each "),
                              ("mod_datetime_name", "Need datetime coordinate"),
                              ("mod_time_name", "Need time coordinate")])
    def test_input_failure(self, del_key, err_msg):
        """ Test for expected failure when missing necessary input from kwargs
        """
        del self.required_kwargs[del_key]

        with pytest.raises(ValueError) as verr:
            match.collect_inst_model_pairs(**self.input_args,
                                           **self.required_kwargs)

        assert verr.value.args[0].find(err_msg) >= 0
    
