"""Unit tests for `pysatModels.utils.extract`."""

from io import StringIO
import logging
import numpy as np
from packaging import version as pack_version
import pytest

import pysat
from pysat.instruments import pysat_testmodel

import pysatModels as ps_mod
import pysatModels.utils.extract as extract


@pytest.mark.skip("input requires a regular grid for the model")

class TestUtilsExtractModObs(object):
    """Unit tests for `utils.extract.extract_modelled_observations`."""

    def setup(self):
        """Set up the unit test environment for each method."""

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
        """Clean up the unit test environment after each method."""

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
        """Test for expected failure with bad input arguments."""

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
        """Test for expected failure with bad kwarg input."""

        self.input_kwargs[bad_key] = bad_val

        with pytest.raises(Exception) as err:
            extract.extract_modelled_observations(*self.input_args,
                                                  **self.input_kwargs)

        assert str(err.value.args[0]).find(err_msg) >= 0
        return

    @pytest.mark.parametrize("sel_val", [["dummy1", "dummy2"], ["dummy1"]])
    def test_good_sel_name(self, sel_val):
        """Test for success with different good selection name inputs."""

        self.input_kwargs = {"sel_name": sel_val,
                             "model_label": self.model_label}
        self.out = extract.extract_modelled_observations(*self.input_args,
                                                         **self.input_kwargs)
        for label in sel_val:
            assert "{:s}_{:s}".format(self.model_label, label) in self.out
        assert len(self.out) == len(np.asarray(sel_val))
        return

    def test_success_w_out_of_bounds(self):
        """Test extraction success for all variables without UT dependence."""

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
        """Test the failure for all model variables already extracted."""

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
        """Test the success for some model variables already extracted."""

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


class TestUtilsExtractInstModView(object):
    """Unit tests for `utils.extract.instrument_view_through_model`."""

    def setup(self):
        """Run before every method to create a clean testing setup."""

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
        """Run after every method to clean up previous testing."""

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


@pytest.mark.skipif(pack_version.Version(pysat.__version__)
                    < pack_version.Version('3.0.2'),
                    reason=''.join(('Requires test model in pysat ',
                                    ' v3.1 or later.')))
class TestUtilsAltitudePressure(object):
    """Unit tests for `utils.extract.instrument_altitude_to_model_pressure`."""

    def setup(self):
        """Set up the unit test environment for each method."""

        self.inst = pysat.Instrument(platform='pysat', name='testing')
        self.model = pysat.Instrument(inst_module=pysat_testmodel,
                                      tag='pressure_levels')
        self.inst.load(date=pysat_testmodel._test_dates[''][''])
        self.model.load(date=pysat_testmodel._test_dates[''][''])

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

        assert str(verr.value.args[0]).find(err_msg) >= 0
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
                    < pack_version.Version('3.1.0'),
                    reason=''.join(('Requires `max_latitude` test Instrument ',
                                    'support in pysat v3.1 or later.')))
class TestUtilsExtractInstModIrregView(object):
    """Unit tests for `utils.extract.instrument_view_irregular_model`."""

    def setup(self):
        """Run before every method to create a clean testing setup."""

        self.inst = pysat.Instrument(platform='pysat', name='testing',
                                     num_samples=3, max_latitude=45.)
        self.model = pysat.Instrument(inst_module=pysat_testmodel,
                                      tag='pressure_levels',
                                      num_samples=96)
        self.inst.load(date=pysat_testmodel._test_dates[''][''])
        self.model.load(date=pysat_testmodel._test_dates[''][''])
        self.model_label = 'tmodel'
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

        del self.inst, self.model, self.input_args, self.out, self.model_label
        del self.in_kwargs

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

        self.input_args[bad_index] = bad_input

        with pytest.raises(ValueError) as verr:
            extract.instrument_view_irregular_model(*self.input_args,
                                                    **self.in_kwargs)

        assert str(verr.value.args[0]).find(err_msg) >= 0

        return

    @pytest.mark.parametrize("bad_key,bad_val,err_msg",
                             [("sel_name", ["unknown_variable"],
                               "Unknown model variable index unknown_variable"),
                              ("sel_name", [], 'Must provide sel_name as a li'),
                              ("model_label", 1, "expected str instance")])
    def test_bad_kwarg_input(self, bad_key, bad_val, err_msg, caplog):
        """Test for expected failure with bad kwarg input."""

        self.in_kwargs[bad_key] = bad_val

        with pytest.raises((ValueError, TypeError)) as err:
            extract.instrument_view_irregular_model(*self.input_args,
                                                    **self.in_kwargs)

        assert str(err.value.args[0]).find(err_msg) >= 0

        return
