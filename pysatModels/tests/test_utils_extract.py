# Author: Angeline Burrell, NRL, 2019

from __future__ import absolute_import, unicode_literals

import logging
from io import StringIO
import numpy as np
import pytest

import pysat
from pysat.instruments import pysat_testmodel

import pysatModels as ps_mod
import pysatModels.utils.extract as extract


@pytest.mark.skip("input requires a regular grid for the model")
class TestUtilsExtractObsViewModel:
    """ Unit tests for utils.extract.sat_view_through_model"""
    def setup(self):
        """ Runs before every method to create a clean testing setup."""
        self.args = [pysat.Instrument(platform='pysat', name='testing',
                                      num_samples=12),
                     pysat.Instrument(platform='pysat', name='testing',
                                      num_samples=6),
                     ['latitude', 'longitude'], ['dummy1', 'dummy2']]

        return

    def teardown(self):
        del self.args

        return

    def test_str_coords(self):
        """ Test string coordinate input"""
        self.args[2] = self.args[2][0]
        extract.instrument_view_through_model(*self.args)

        for label in self.args[3]:
            assert "model_{:s}".format(label) in self.args[1].data.columns

        return

    def test_str_obs(self):
        """ Test string model observation input"""
        self.args[3] = self.args[3][0]
        extract.instrument_view_through_model(*self.args)

        assert "model_{:s}".format(self.args[3]) in self.args[1].data.columns

        return


class TestUtilsExtractModObs:
    """ Unit tests for utils.extract.extract_modelled_observations """

    def setup(self):
        """Runs before every method to create a clean testing setup."""
        self.inst = pysat.Instrument(platform='pysat', name='testing')
        self.model = pysat.Instrument(inst_module=pysat_testmodel)
        self.inst.load(date=pysat_testmodel._test_dates[''][''])
        self.model.load(date=pysat_testmodel._test_dates[''][''])
        self.model_label = 'tmodel'
        self.input_args = [self.inst, self.model.data,
                           ["longitude", "latitude", "altitude"],
                           ["longitude", "latitude", "altitude"],
                           "time", "time", ["deg", "deg", "km"]]
        self.input_kwargs = {"sel_name":
                             [kk for kk in self.model.data.data_vars
                              if len([dd for dd
                                      in self.model.data.data_vars[kk].dims
                                      if dd in self.input_args[3]])]}
        self.out = []
        self.log_capture = StringIO()
        ps_mod.logger.addHandler(logging.StreamHandler(self.log_capture))
        ps_mod.logger.setLevel(logging.INFO)

        return

    def teardown(self):
        """Runs after every method to clean up previous testing."""
        del self.inst, self.model, self.input_args, self.out, self.model_label
        del self.input_kwargs, self.log_capture

        return

    @pytest.mark.parametrize("bad_index,bad_input,err_msg",
                             [(2, [], "Must provide instrument location"),
                              (2, ['glon', 'latitude', 'altitude'],
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

        return

    @pytest.mark.parametrize("bad_key,bad_val,err_msg",
                             [("sel_name", ["altitude"],
                               "No model data keys to interpolate"),
                              ("method", "not_a_method",
                               "interpn only understands the methods"),
                              ("model_label", 1, "Unknown format code "),
                              ("time_method", "fun", "unknown time method"),
                              ("pair_method", "fun",
                               "unknown pairing method")])
    def test_bad_kwarg_input(self, bad_key, bad_val, err_msg):
        """ Test for expected failure with bad kwarg input """
        self.input_kwargs[bad_key] = bad_val

        with pytest.raises(Exception) as err:
            extract.extract_modelled_observations(*self.input_args,
                                                  **self.input_kwargs)

        assert str(err.value.args[0]).find(err_msg) >= 0

        return

    @pytest.mark.parametrize("sel_val", [["dummy1", "dummy2"], ["dummy1"]])
    def test_good_sel_name(self, sel_val):
        """ Test for success with different good selection name inputs"""
        self.input_kwargs = {"sel_name": sel_val,
                             "model_label": self.model_label}
        self.out = extract.extract_modelled_observations(*self.input_args,
                                                         **self.input_kwargs)
        for label in sel_val:
            assert "{:s}_{:s}".format(self.model_label, label) in self.out
        assert len(self.out) == len(np.asarray(sel_val))

        return

    def test_success_w_out_of_bounds(self):
        """ Test the extraction success for all variables without UT dependence
        """
        self.input_kwargs["model_label"] = self.model_label
        self.out = extract.extract_modelled_observations(*self.input_args,
                                                         **self.input_kwargs)
        lout = self.log_capture.getvalue()

        assert lout.find('One of the requested xi is out of bounds') >= 0

        for label in self.input_kwargs['sel_name']:
            if label not in self.input_args[3]:
                # Test each of the extracted model data columns
                tcol = "{:s}_{:s}".format(self.model_label, label)
                assert tcol in self.out
                assert tcol in self.inst.data.columns
                assert (self.inst.data[self.input_args[2][0]].shape
                        == self.inst.data[tcol].shape)
                assert len(self.inst.data[tcol][
                    ~np.isnan(self.inst.data[tcol])]) > 0

        return

    def test_failure_for_already_ran_data(self):
        """ Test the failure for all model variables already extracted """

        self.input_kwargs["model_label"] = self.model_label

        # Run everything successfully once
        extract.extract_modelled_observations(*self.input_args,
                                              **self.input_kwargs)

        # Run everything again, raising a value error
        with pytest.raises(ValueError) as err:
            extract.extract_modelled_observations(*self.input_args,
                                                  **self.input_kwargs)

            self.out = self.log_capture.getvalue()
            assert self.out.find('model data already interpolated') >= 0

        assert str(err.value.args[0]).find(
            'instrument object already contains all model data') >= 0

        return

    def test_success_for_some_already_ran_data(self):
        """ Test the success for some model variables already extracted """

        all_sel = list(self.input_kwargs['sel_name'])
        self.input_kwargs['model_label'] = self.model_label

        # Run through twice
        for i, selected in enumerate([all_sel[2:], all_sel]):
            self.input_kwargs['sel_name'] = selected
            self.out = extract.extract_modelled_observations(
                *self.input_args, **self.input_kwargs)
            lout = self.log_capture.getvalue()

        assert lout.find('model data already interpolated') >= 0

        for label in all_sel:
            if label not in self.input_args[3]:
                # Test each of the extracted model data columns
                tcol = "{:s}_{:s}".format(self.model_label, label)
                assert tcol in self.inst.data.columns
                assert (self.inst.data[self.input_args[2][0]].shape
                        == self.inst.data[tcol].shape)
                assert len(self.inst.data[tcol][
                    ~np.isnan(self.inst.data[tcol])]) > 0

        return


class TestUtilsExtractInstModView:
    """Unit tests for `utils.extract.instrument_view_through_model`."""

    def setup(self):
        """Runs before every method to create a clean testing setup."""
        self.inst = pysat.Instrument(platform='pysat', name='testing')
        self.model = pysat.Instrument(inst_module=pysat_testmodel)
        self.inst.load(date=pysat_testmodel._test_dates[''][''])
        self.model.load(date=pysat_testmodel._test_dates[''][''])
        self.model_label = 'tmodel'
        self.input_args = [self.inst, self.model.data,
                           ["latitude", "longitude", "altitude"],
                           ["latitude", "longitude", "altitude"],
                           "time", "time", ["deg", "deg", "km"]]
        self.input_kwargs = {"sel_name":
                             [kk for kk in self.model.data.data_vars
                              if len([dd for dd
                                      in self.model.data.data_vars[kk].dims
                                      if dd in self.input_args[3]])
                              and (len(self.model.data.data_vars[kk].dims)
                                   == len(self.input_args[3]) + 1)]}
        self.input_kwargs['methods'] = ['linear'] * len(
            self.input_kwargs['sel_name'])
        self.out = []
        self.log_capture = StringIO()
        ps_mod.logger.addHandler(logging.StreamHandler(self.log_capture))
        ps_mod.logger.setLevel(logging.INFO)

        return

    def teardown(self):
        """Runs after every method to clean up previous testing."""
        del self.inst, self.model, self.input_args, self.out, self.model_label
        del self.input_kwargs, self.log_capture

        return

    def test_standard_call(self):
        """Test for successful interpolation."""

        self.out = extract.instrument_view_through_model(*self.input_args,
                                                         **self.input_kwargs)
        for name in self.input_kwargs['sel_name']:
            assert ''.join(('model_', name)) in self.inst.data

        for name in self.out:
            assert name in self.inst.data

        return

    @pytest.mark.parametrize("bad_index,bad_input,err_msg",
                             [(2, [], 'Must provide inst_name as a list'),
                              (3, [], 'Must provide mod_name as a list'),
                              (2, ['glon', 'latitude', 'altitude'],
                               "Unknown instrument location index"),
                              (3, ['hi'], "Must provide the same number"),
                              (6, [], "Must provide units for each "),
                              (4, "naname", "unknown model name for datetime"),
                              (5, "naname", "Unknown model time coordinate")])
    def test_bad_arg_input(self, bad_index, bad_input, err_msg):
        """Test for expected failure with bad input arguments."""
        self.input_args[bad_index] = bad_input

        with pytest.raises(ValueError) as verr:
            extract.instrument_view_through_model(*self.input_args,
                                                  **self.input_kwargs)

        assert str(verr.value.args[0]).find(err_msg) >= 0

        return

    @pytest.mark.parametrize("bad_key,bad_val,err_msg",
                             [("sel_name", ["unknown_variable"],
                               "unknown_variable is not a valid model variabl"),
                              ("methods", ["not_a_method"],
                               'Methods only supports "linear" or "nearest".'),
                              ("methods", ['linear', 'linear', 'linear'],
                               "Must provide interpolation selection"),
                              ("model_label", 1, "expected str instance")])
    def test_bad_kwarg_input(self, bad_key, bad_val, err_msg):
        """Test for expected failure with bad kwarg input."""
        self.input_kwargs[bad_key] = bad_val

        with pytest.raises((ValueError, TypeError)) as err:
            extract.instrument_view_through_model(*self.input_args,
                                                  **self.input_kwargs)

        assert str(err.value.args[0]).find(err_msg) >= 0

        return

    def test_failure_for_already_ran_data(self):
        """Test the failure for all model variables already extracted."""

        self.input_kwargs["model_label"] = self.model_label

        # Run everything successfully once
        extract.instrument_view_through_model(*self.input_args,
                                              **self.input_kwargs)

        # Run everything again, raising a value error
        with pytest.raises(ValueError) as err:
            extract.instrument_view_through_model(*self.input_args,
                                                  **self.input_kwargs)

            self.out = self.log_capture.getvalue()
            assert self.out.find('model data already interpolated') >= 0

        assert str(err.value.args[0]).find(
            'instrument object already contains all model data') >= 0

        return

    def test_success_for_some_already_ran_data(self):
        """Test the success for some model variables already extracted."""

        all_sel = list(self.input_kwargs['sel_name'])
        all_sel.append('dummy72')
        self.input_kwargs['model_label'] = self.model_label

        self.model['dummy72'] = self.model['dummy2']

        # Run through twice
        for i, selected in enumerate([all_sel[1:], all_sel]):
            self.input_kwargs['sel_name'] = selected
            self.input_kwargs['methods'] = ['linear'] * len(selected)
            self.out = extract.instrument_view_through_model(
                *self.input_args, **self.input_kwargs)

            lout = self.log_capture.getvalue()

        assert lout.find('model data already interpolated') >= 0

        for label in all_sel:
            if label not in self.input_args[3]:
                # Test each of the extracted model data columns
                tcol = "{:s}_{:s}".format(self.model_label, label)
                assert tcol in self.inst.data.columns
                assert (self.inst.data[self.input_args[2][0]].shape
                        == self.inst.data[tcol].shape)
                assert len(self.inst.data[tcol][
                    ~np.isnan(self.inst.data[tcol])]) > 0

        return
