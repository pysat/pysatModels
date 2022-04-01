#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2022, pysat development team
# Full license can be found in License.md
# -----------------------------------------------------------------------------
"""Unit tests for `pysatModels.utils.extract`."""

import logging
import numpy as np
from packaging import version as pack_version
import pytest

import pysat
from pysat.instruments import pysat_testmodel

import pysatModels.utils.extract as extract


class TestUtilsExtractInstThroughMod(object):
    """Unit tests for `instrument_view_through_model`."""

    def setup(self):
        """Set up the unit test environment for each method."""

        self.inst = pysat.Instrument(platform='pysat', name='testing')
        self.model = pysat.Instrument(inst_module=pysat_testmodel, tag='')

        # Load the data in the instruments
        load_kwargs = {'date': pysat_testmodel._test_dates['']['']}
        if(pack_version.Version(pysat.__version__)
           > pack_version.Version('3.0.1')):
            load_kwargs['use_header'] = True

        self.inst.load(**load_kwargs)
        self.model.load(**load_kwargs)

        # Define the inputs
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
        self.input_kwargs['methods'] = [
            'linear' for i in self.input_kwargs['sel_name']]
        self.func = extract.instrument_view_through_model
        self.out = []
        return

    def teardown(self):
        """Clean up the unit test environment after each method."""

        del self.inst, self.model, self.input_args, self.out, self.input_kwargs
        return

    def eval_output(self):
        """Evaluate successful output."""

        if "model_label" not in self.input_kwargs:
            self.input_kwargs['model_label'] = 'model'

        for label in self.input_kwargs['sel_name']:
            if label not in self.input_args[3]:
                # Test each of the extracted model data columns
                tcol = "{:s}_{:s}".format(self.input_kwargs['model_label'],
                                          label)
                assert tcol in self.out, \
                    "missing expected output variable {:} in {:}".format(
                        tcol, self.out)
                assert tcol in self.inst.variables, \
                    "missing expected Instrument variable {:}".format(tcol)
                assert self.inst[label].shape == self.inst[tcol].shape, \
                    "unexpected output data ({:}) shape {:} != {:}".format(
                        tcol, self.inst[label].shape, self.inst[tcol].shape)
                assert len(self.inst[tcol].where(
                    ~np.isnan(self.inst.data[tcol])).values) > 0, \
                    "unexpected number of finite matching data values."
        return

    def test_standard_call(self):
        """Test for successful interpolation."""

        # Run the basic function call
        self.out = self.func(*self.input_args, **self.input_kwargs)

        # Evaluate the output
        self.eval_output()
        return

    @pytest.mark.parametrize("bad_index, bad_input, err_msg", [
        (2, [], 'Must provide inst_name as a list'),
        (3, [], 'Must provide mod_name as a list'),
        (2, ['glon', 'latitude', 'altitude'],
         "Unknown instrument location index"),
        (3, ['hi'], "Must provide the same number"),
        (6, [], "Must provide units for each "),
        (4, "naname", "unknown model name for datetime"),
        (5, "naname", "Unknown model time coordinate")])
    def test_bad_arg_input(self, bad_index, bad_input, err_msg):
        """Test for expected failure with bad input arguments.

        Parameters
        ----------
        bad_index : int
            Input argument index to replace
        bad_input : any type
            Bad input value
        err_msg : str
            Expected error message

        """

        self.input_args[bad_index] = bad_input

        with pytest.raises(ValueError) as verr:
            self.func(*self.input_args, **self.input_kwargs)

        assert str(verr).find(err_msg) >= 0
        return

    @pytest.mark.parametrize("bad_key, bad_val, err_msg", [
        ("sel_name", ["unknown_variable"],
         "unknown_variable is not a valid model variabl"),
        ("methods", ["not_a_method"],
         'Methods only supports "linear" or "nearest".'),
        ("methods", ['linear', 'linear', 'linear'],
         "Must provide interpolation selection"),
        ("model_label", 1, "expected str instance")])
    def test_bad_kwarg_input(self, bad_key, bad_val, err_msg):
        """Test for expected failure with bad kwarg input.

        Parameters
        ----------
        bad_key : str
            Input kwarg key to replace
        bad_val : any type
            Bad input value
        err_msg : str
            Expected error message

        """

        self.input_kwargs[bad_key] = bad_val

        with pytest.raises((ValueError, TypeError)) as err:
            self.func(*self.input_args, **self.input_kwargs)

        assert str(err).find(err_msg) >= 0
        return

    def test_failure_for_already_ran_data(self, caplog):
        """Test the failure for all model variables already extracted."""

        # Run everything successfully once
        self.func(*self.input_args, **self.input_kwargs)

        # Run everything again, raising a value error
        with caplog.at_level(logging.INFO, "pysatModels"):
            with pytest.raises(ValueError) as err:
                extract.extract_modelled_observations(*self.input_args,
                                                      **self.input_kwargs)

        # Evaluate the logging messages and error messages
        self.out = caplog.text
        assert self.out.find('model data already interpolated') >= 0

        assert str(err.value.args[0]).find(
            'instrument object already contains all model data') >= 0
        return

    def test_success_for_some_already_ran_data(self, caplog):
        """Test the success for some model variables already extracted."""

        all_sel = list(self.input_kwargs['sel_name'])
        all_sel.append('dummy72')

        self.model['dummy72'] = self.model['dummy2']

        # Run through twice
        with caplog.at_level(logging.INFO, logger='pysatModels'):
            for i, selected in enumerate([all_sel[1:], all_sel]):
                self.input_kwargs['sel_name'] = selected
                self.input_kwargs['methods'] = ['linear'] * len(selected)
                self.func(*self.input_args, **self.input_kwargs)

            self.out = caplog.text
            assert self.out.find('model data already interpolated') >= 0

        # Evaluate output
        self.input_kwargs['sel_name'] = pysat.utils.listify(all_sel[1])
        self.eval_output()
        return


@pytest.mark.skipif(pack_version.Version(pysat.__version__)
                    <= pack_version.Version('3.0.1'),
                    reason=''.join(('Requires test model in pysat ',
                                    ' v3.0.2 or later.')))
class TestUtilsExtractModObs(TestUtilsExtractInstThroughMod):
    """Unit tests for `utils.extract.extract_modelled_observations`."""

    def setup(self):
        """Set up the unit test environment for each method."""

        self.inst = pysat.Instrument(platform='pysat', name='testing')
        self.model = pysat.Instrument(inst_module=pysat_testmodel, tag='')

        # Load the data in the instruments
        load_kwargs = {'date': pysat_testmodel._test_dates['']['']}
        if(pack_version.Version(pysat.__version__)
           > pack_version.Version('3.0.1')):
            load_kwargs['use_header'] = True

        self.inst.load(**load_kwargs)
        self.model.load(**load_kwargs)

        # Define the inputs
        self.input_args = [self.inst, self.model.data,
                           ["longitude", "latitude", "altitude"],
                           ["longitude", "latitude", "altitude"],
                           "time", "time", ["deg", "deg", "km"]]
        self.input_kwargs = {"sel_name":
                             [kk for kk in self.model.data.data_vars
                              if len([dd for dd
                                      in self.model.data.data_vars[kk].dims
                                      if dd in self.input_args[3]])]}
        self.func = extract.extract_modelled_observations
        self.out = []
        return

    def teardown(self):
        """Clean up the unit test environment after each method."""

        del self.inst, self.model, self.input_args, self.out, self.input_kwargs
        return

    @pytest.mark.parametrize("bad_index, bad_input, err_msg", [
        (2, [], "Must provide instrument location"),
        (2, ['glon', 'latitude', 'altitude'],
         "Unknown instrument location index"),
        (3, [], "Must provide the same number"),
        (6, [], "Must provide units for each "),
        (4, "naname", "unknown model name for datetime"),
        (5, "naname", "Unknown model time coordinate")])
    def test_bad_arg_input(self, bad_index, bad_input, err_msg):
        """Test for expected failure with bad input arguments.

        Parameters
        ----------
        bad_index : int
            Input argument index to replace
        bad_input : any type
            Bad input value
        err_msg : str
            Expected error message

        """

        self.input_args[bad_index] = bad_input

        with pytest.raises(ValueError) as verr:
            self.func(*self.input_args, **self.input_kwargs)

        assert str(verr).find(err_msg) >= 0
        return

    @pytest.mark.parametrize("bad_key, bad_val, err_msg", [
        ("sel_name", ["altitude"], "No model data keys to interpolate"),
        ("method", "not_a_method", "interpn only understands the methods"),
        ("model_label", 1, "Unknown format code "),
        ("time_method", "fun", "unknown time method"),
        ("pair_method", "fun", "unknown pairing method")])
    def test_bad_kwarg_input(self, bad_key, bad_val, err_msg):
        """Test for expected failure with bad kwarg input.

        Parameters
        ----------
        bad_key : str
            Input kwarg key to replace
        bad_val : any type
            Bad input value
        err_msg : str
            Expected error message

        """

        self.input_kwargs[bad_key] = bad_val

        with pytest.raises((ValueError, TypeError)) as err:
            self.func(*self.input_args, **self.input_kwargs)

        assert str(err).find(err_msg) >= 0
        return

    @pytest.mark.parametrize("sel_val", [["dummy1", "dummy2"], ["dummy1"]])
    def test_good_sel_name(self, sel_val):
        """Test for success with different good selection name inputs.

        Parameters
        ----------
        sel_val : list
            List of selection variable names

        """

        self.input_kwargs = {"sel_name": sel_val, "model_label": "test_model"}
        self.out = self.func(*self.input_args, **self.input_kwargs)
        self.eval_output()
        return

    def test_success_w_out_of_bounds(self, caplog):
        """Test extraction success for all variables without UT dependence."""

        with caplog.at_level(logging.INFO, "pysatModels"):
            self.out = self.func(*self.input_args, **self.input_kwargs)

        assert caplog.text.find('One of the requested xi is out of bounds') >= 0
        self.eval_output()
        return

    def test_success_for_some_already_ran_data(self, caplog):
        """Test the success for some model variables already extracted."""

        all_sel = list(self.input_kwargs['sel_name'])

        # Run through twice
        with caplog.at_level(logging.INFO, 'pysatModels'):
            for i, selected in enumerate([all_sel[2:], all_sel]):
                self.input_kwargs['sel_name'] = selected
                self.out.extend(self.func(*self.input_args,
                                          **self.input_kwargs))

        # Test the logging message
        assert caplog.text.find('model data already interpolated') >= 0

        # Evaluate the output
        self.input_kwargs['sel_name'] = all_sel
        self.eval_output()
        return


@pytest.mark.skipif(pack_version.Version(pysat.__version__)
                    <= pack_version.Version('3.0.1'),
                    reason=''.join(('Requires test model in pysat ',
                                    ' v3.0.2 or later.')))
class TestUtilsExtractModObsXarray(TestUtilsExtractModObs):
    """Xarray unit tests for `utils.extract.extract_modelled_observations`."""

    def setup(self):
        """Set up the unit test environment for each method."""

        self.inst = pysat.Instrument(platform='pysat', name='testing_xarray')
        self.model = pysat.Instrument(inst_module=pysat_testmodel, tag='')

        # Load the data in the instruments
        load_kwargs = {'date': pysat_testmodel._test_dates['']['']}
        if(pack_version.Version(pysat.__version__)
           > pack_version.Version('3.0.1')):
            load_kwargs['use_header'] = True

        self.inst.load(**load_kwargs)
        self.model.load(**load_kwargs)

        # Define the inputs
        self.input_args = [self.inst, self.model.data,
                           ["longitude", "latitude", "altitude"],
                           ["longitude", "latitude", "altitude"],
                           "time", "time", ["deg", "deg", "km"]]
        self.input_kwargs = {"sel_name":
                             [kk for kk in self.model.data.data_vars
                              if len([dd for dd
                                      in self.model.data.data_vars[kk].dims
                                      if dd in self.input_args[3]])]}
        self.func = extract.extract_modelled_observations
        self.out = []
        return

    def teardown(self):
        """Clean up the unit test environment after each method."""

        del self.inst, self.model, self.input_args, self.out, self.input_kwargs
        return


@pytest.mark.skipif(pack_version.Version(pysat.__version__)
                    <= pack_version.Version('3.0.1'),
                    reason=''.join(('Requires test model in pysat ',
                                    ' v3.0.2 or later.')))
class TestUtilsExtractModObsXarray2D(TestUtilsExtractModObs):
    """Xarray unit tests for `utils.extract.extract_modelled_observations`."""

    def setup(self):
        """Set up the unit test environment for each method."""

        self.inst = pysat.Instrument(platform='pysat', name='testing2d_xarray')
        self.model = pysat.Instrument(inst_module=pysat_testmodel, tag='')

        # Load the data in the instruments
        load_kwargs = {'date': pysat_testmodel._test_dates['']['']}
        if(pack_version.Version(pysat.__version__)
           > pack_version.Version('3.0.1')):
            load_kwargs['use_header'] = True

        self.inst.load(**load_kwargs)
        self.model.load(**load_kwargs)

        # Define the inputs
        self.input_args = [self.inst, self.model.data,
                           ["longitude", "latitude", "altitude"],
                           ["longitude", "latitude", "altitude"],
                           "time", "time", ["deg", "deg", "km"]]
        self.input_kwargs = {"sel_name":
                             [kk for kk in self.model.data.data_vars
                              if len([dd for dd
                                      in self.model.data.data_vars[kk].dims
                                      if dd in self.input_args[3]])]}
        self.func = extract.extract_modelled_observations
        self.out = []
        return

    def teardown(self):
        """Clean up the unit test environment after each method."""

        del self.inst, self.model, self.input_args, self.out, self.input_kwargs
        return


class TestUtilsExtractInstModViewXarray(TestUtilsExtractInstThroughMod):
    """Xarray unit tests for `instrument_view_through_model`."""

    def setup(self):
        """Run before every method to create a clean testing setup."""

        self.inst = pysat.Instrument(platform='pysat', name='testing2d_xarray')
        self.model = pysat.Instrument(inst_module=pysat_testmodel)

        # Load the data in the instruments
        load_kwargs = {'date': pysat_testmodel._test_dates['']['']}
        if(pack_version.Version(pysat.__version__)
           > pack_version.Version('3.0.1')):
            load_kwargs['use_header'] = True

        self.inst.load(**load_kwargs)
        self.model.load(**load_kwargs)

        # Define the inputs
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
        self.func = extract.instrument_view_through_model
        self.out = []

        return

    def teardown(self):
        """Run after every method to clean up previous testing."""

        del self.inst, self.model, self.input_args, self.out, self.input_kwargs

        return


@pytest.mark.skipif(pack_version.Version(pysat.__version__)
                    <= pack_version.Version('3.0.1'),
                    reason=''.join(('Requires test model in pysat ',
                                    ' v3.0.2 or later.')))
class TestUtilsAltitudePressure(object):
    """Unit tests for `utils.extract.instrument_altitude_to_model_pressure`."""

    def setup(self):
        """Set up the unit test environment for each method."""

        self.inst = pysat.Instrument(platform='pysat', name='testing')
        self.model = pysat.Instrument(inst_module=pysat_testmodel,
                                      tag='pressure_levels')

        # Load the data in the instruments
        load_kwargs = {'date': pysat_testmodel._test_dates['']['']}
        if(pack_version.Version(pysat.__version__)
           > pack_version.Version('3.0.1')):
            load_kwargs['use_header'] = True

        self.inst.load(**load_kwargs)
        self.model.load(**load_kwargs)

        # Define the inputs
        self.input_args = [self.inst, self.model.data,
                           ["altitude", "latitude", "longitude"],
                           ["ilev", "latitude", "longitude"],
                           "time", "time", ['', "deg", "deg"],
                           'altitude', 'altitude', 'cm']
        self.input_kwargs = {}

        self.out = []
        return

    def teardown(self):
        """Clean up the unit test environment after each method."""

        del self.inst, self.model, self.input_args, self.out
        del self.input_kwargs
        return

    @pytest.mark.parametrize("bad_index,bad_input,err_msg",
                             [(2, [], "Must provide inst_name as a list"),
                              (2, [''], 'Must provide the same number of'),
                              (2, ['glon', 'latitude', 'altitude'],
                               "Unknown instrument location index"),
                              (3, [], "Must provide mod_name as a list"),
                              (3, [''], "Must provide the same number"),
                              (6, [], "Must provide units for each "),
                              (4, "naname", "unknown model name for datetime"),
                              (5, "naname", "Unknown model time coordinate"),
                              (7, 'navar', 'Unknown Instrument altitude key'),
                              (8, 'navar', 'Unknown Model altitude key')])
    def test_bad_arg_input(self, bad_index, bad_input, err_msg):
        """Test for expected failure with bad input arguments."""

        self.input_args[bad_index] = bad_input

        with pytest.raises(ValueError) as verr:
            extract.instrument_altitude_to_model_pressure(*self.input_args)

        assert str(verr).find(err_msg) >= 0
        return

    @pytest.mark.parametrize("tol_val", [10., 1., 0.1])
    @pytest.mark.parametrize("scale_val", [1000., 100., 50.])
    def test_good_translation_over_tolerance_and_scale(self, tol_val,
                                                       scale_val):
        """Test for success with different altitude tolerances and scales.

        Parameters
        ----------
        tol_val : float
            Tolerance for altitude matching
        scale_val : float
            Scale height used to normalize altitude differences

        """

        self.input_kwargs = {"tol": tol_val,
                             "scale": scale_val}
        self.out = extract.instrument_altitude_to_model_pressure(
            *self.input_args,
            **self.input_kwargs)

        # Calculate difference in altitude (Instrument and values extracted
        # from Model) and ensure it is less than specified tolerance.
        alt_diff = np.abs(self.inst[self.out[0]] - self.inst['altitude'])
        assert np.all(alt_diff <= tol_val)
        return

    def test_updated_metadata(self):
        """Test new pressure metadata is present in Instrument."""

        self.out = extract.instrument_altitude_to_model_pressure(
            *self.input_args)

        assert self.inst.meta['model_altitude', 'units'] == 'km'

        test_str = 'Interpolated Model altitude'
        assert self.inst.meta['model_altitude', 'notes'].find(test_str) >= 0

        test_str = 'Interpolated Model pressure'
        assert self.inst.meta['model_pressure', 'notes'].find(test_str) >= 0

        return

    def test_alternate_output_names(self):
        """Test alternate output labels work as expected."""
        self.input_kwargs = {'inst_out_alt': 'alter_altitude',
                             'inst_out_pres': 'alter_pressure'}

        self.out = extract.instrument_altitude_to_model_pressure(
            *self.input_args, **self.input_kwargs)

        assert 'alter_altitude' in self.inst.variables
        assert 'alter_altitude' == self.out[0]
        assert 'alter_pressure' in self.inst.variables
        assert 'alter_pressure' == self.out[1]

        alt_diff = np.abs(self.inst[self.out[0]] - self.inst['altitude'])

        # Test that the altitude difference is less than or equal to the
        # default tolerance value for the function.
        assert np.all(alt_diff <= 1.0)

        return


@pytest.mark.skipif(pack_version.Version(pysat.__version__)
                    <= pack_version.Version('3.0.1'),
                    reason=''.join(('Requires test model in pysat ',
                                    ' v3.0.2 or later.')))
class TestUtilsAltitudePressureXarray(TestUtilsAltitudePressure):
    """Xarray unit tests for `instrument_altitude_to_model_pressure`."""

    def setup(self):
        """Set up the unit test environment for each method."""

        self.inst = pysat.Instrument(platform='pysat', name='testing2d_xarray')
        self.model = pysat.Instrument(inst_module=pysat_testmodel,
                                      tag='pressure_levels')

        # Load the data in the instruments
        load_kwargs = {'date': pysat_testmodel._test_dates['']['']}
        if(pack_version.Version(pysat.__version__)
           > pack_version.Version('3.0.1')):
            load_kwargs['use_header'] = True

        self.inst.load(**load_kwargs)
        self.model.load(**load_kwargs)

        # Define the inputs
        self.input_args = [self.inst, self.model.data,
                           ["altitude", "latitude", "longitude"],
                           ["ilev", "latitude", "longitude"],
                           "time", "time", ['', "deg", "deg"],
                           'altitude', 'altitude', 'cm']
        self.input_kwargs = {}

        self.out = []
        return

    def teardown(self):
        """Clean up the unit test environment after each method."""

        del self.inst, self.model, self.input_args, self.out
        del self.input_kwargs
        return


@pytest.mark.skipif(pack_version.Version(pysat.__version__)
                    <= pack_version.Version('3.0.1'),
                    reason=''.join(('Requires `max_latitude` test Instrument ',
                                    'support in pysat v3.0.2 or later.')))
class TestUtilsExtractInstModIrregView(object):
    """Unit tests for `utils.extract.instrument_view_irregular_model`."""

    def setup(self):
        """Run before every method to create a clean testing setup."""

        self.inst = pysat.Instrument(platform='pysat', name='testing',
                                     num_samples=3, max_latitude=45.)
        self.model = pysat.Instrument(inst_module=pysat_testmodel,
                                      tag='pressure_levels',
                                      num_samples=96)

        # Load the data in the instruments
        load_kwargs = {'date': pysat_testmodel._test_dates['']['']}
        if(pack_version.Version(pysat.__version__)
           > pack_version.Version('3.0.1')):
            load_kwargs['use_header'] = True

        self.inst.load(**load_kwargs)
        self.model.load(**load_kwargs)

        # Define the inputs
        self.input_args = [self.inst, self.model.data,
                           ["altitude", "latitude", "longitude"],
                           ["ilev", "latitude", "longitude"],
                           "time", ["cm", "deg", "deg"], "ilev",
                           "altitude", [50., 10., 10.]]
        self.in_kwargs = {"sel_name": ["dummy_drifts", "altitude"]}
        self.out = []

        return

    def teardown(self):
        """Run after every method to clean up previous testing."""

        del self.inst, self.model, self.input_args, self.out, self.in_kwargs

        return

    def test_standard_call(self):
        """Test for successful interpolation."""

        self.out = extract.interp_inst_w_irregular_model_coord(*self.input_args,
                                                               **self.in_kwargs)
        for name in self.in_kwargs['sel_name']:
            assert ''.join(('model_', name)) in self.inst.data

        for name in self.out:
            assert name in self.inst.data

        # Ensure values are all finite
        for name in self.out:
            assert np.all(np.isfinite(self.inst[name]))

        return

    @pytest.mark.parametrize("bad_index,bad_input,err_msg",
                             [(2, [], 'Must provide inst_name as a list'),
                              (3, [], 'Must provide mod_name as a list'),
                              (2, ['glon', 'latitude', 'altitude'],
                               "Unknown instrument location index"),
                              (3, ['hi'], "Must provide the same number"),
                              (5, [], "Must provide units for each "),
                              (4, "naname", "unknown model name for datetime"),
                              (6, "lev", "mod_reg_dim must be a coordinate "),
                              (3, "lev", "mod_name must contain coordinate"),
                              (7, "not", "Unknown irregular model"),
                              (7, "lev", "Coordinate dimensions must"),
                              (8, [], "Must provide mod_var_delta "),
                              (8, ['hi'], 'Must provide the same number of')])
    def test_bad_arg_input(self, bad_index, bad_input, err_msg):
        """Test for expected failure with bad input arguments."""

        # Set the input arguments
        self.input_args[bad_index] = bad_input

        # Raise the desired error
        with pytest.raises(ValueError) as verr:
            extract.interp_inst_w_irregular_model_coord(*self.input_args,
                                                        **self.in_kwargs)

        # Evaluate the error message
        assert str(verr).find(err_msg) >= 0
        return

    @pytest.mark.parametrize("bad_key,bad_val,err_msg",
                             [("sel_name", ["unknown_variable"],
                               "Unknown model variable index unknown_variable"),
                              ("sel_name", [], 'Must provide sel_name as a li'),
                              ("model_label", 1, "expected str instance")])
    def test_bad_kwarg_input(self, bad_key, bad_val, err_msg):
        """Test for expected failure with bad kwarg input."""

        self.in_kwargs[bad_key] = bad_val

        with pytest.raises((ValueError, TypeError)) as err:
            extract.interp_inst_w_irregular_model_coord(*self.input_args,
                                                        **self.in_kwargs)

        assert str(err).find(err_msg) >= 0

        return


@pytest.mark.skipif(pack_version.Version(pysat.__version__)
                    <= pack_version.Version('3.0.1'),
                    reason=''.join(('Requires `max_latitude` test Instrument ',
                                    'support in pysat v3.0.2 or later.')))
class TestUtilsExtractInstModIrregViewXarray(TestUtilsExtractInstModIrregView):
    """Xarray unit tests for `instrument_view_irregular_model`."""

    def setup(self):
        """Run before every method to create a clean testing setup."""

        self.inst = pysat.Instrument(platform='pysat', name='testing2d_xarray',
                                     num_samples=3, max_latitude=45.)
        self.model = pysat.Instrument(inst_module=pysat_testmodel,
                                      tag='pressure_levels',
                                      num_samples=96)

        # Load the data in the instruments
        load_kwargs = {'date': pysat_testmodel._test_dates['']['']}
        if(pack_version.Version(pysat.__version__)
           > pack_version.Version('3.0.1')):
            load_kwargs['use_header'] = True

        self.inst.load(**load_kwargs)
        self.model.load(**load_kwargs)

        # Define the inputs
        self.input_args = [self.inst, self.model.data,
                           ["altitude", "latitude", "longitude"],
                           ["ilev", "latitude", "longitude"],
                           "time", ["cm", "deg", "deg"], "ilev",
                           "altitude", [50., 10., 10.]]
        self.in_kwargs = {"sel_name": ["dummy_drifts", "altitude"]}
        self.out = []

        return

    def teardown(self):
        """Run after every method to clean up previous testing."""

        del self.inst, self.model, self.input_args, self.out, self.in_kwargs

        return
