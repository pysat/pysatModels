"""Unit tests for `pysatModels.utils.compare`."""

import logging
import numpy as np
import pytest
import xarray as xr

import pysatModels.utils.compare as compare


class TestUtilsCompare(object):
    """Unit tests for utils.compare."""

    def setup_method(self):
        """Set up the unit test environment for each method."""

        coords = {"lat": np.arange(-90, 90, 10), "alt": np.arange(200, 800, 10)}
        temp = np.zeros(shape=(len(coords['lat']), len(coords['alt'])))
        self.dat_key = 'sat_Ti'
        self.mod_key = 'mod_Ti'
        self.units = "C"
        self.pairs = xr.Dataset({self.dat_key: (["lat", "alt"], temp + 100.0,
                                                {"units": self.units}),
                                 self.mod_key: (["lat", "alt"], temp + 150.0,
                                                {"units": self.units})},
                                coords=coords)
        self.input_args = [self.pairs, [self.dat_key], [self.mod_key], ["all"]]
        self.out_stat = dict()
        self.out_units = dict()
        return

    def teardown_method(self):
        """Clean up the unit test environment after each method."""

        del self.dat_key, self.mod_key, self.units, self.pairs, self.input_args
        del self.out_stat, self.out_units
        return

    def eval_comp_output(self):
        """Evaluate the comparison output."""

        # Evaluate the output data dictionary structure
        assert self.out_stat.keys() == self.out_units.keys()
        assert self.dat_key in self.out_stat.keys()
        assert self.mod_key not in self.out_stat.keys()
        assert len(self.out_stat.keys()) == len(self.input_args[1])

        # Evaluate the unit output
        assert self.out_units[self.dat_key] == self.units

        # Evaluate the number of returned statistics (all requested)
        is_float = ['bias', 'meanPercentageError', 'medianLogAccuracy',
                    'symmetricSignedBias', 'meanSquaredError', 'RMSE',
                    'meanAbsError', 'medAbsError', 'nRMSE', 'MASE',
                    'medSymAccuracy', 'meanAPE', 'medAPE']
        not_float = ['logAccuracy', 'scaledError', 'forecastError', 'percError',
                     'absPercError']

        assert (len(self.out_stat[self.dat_key].keys())
                == len(is_float) + len(not_float)), \
            "unexpected number of statistics returned: {:}".format(
                [kk for kk in self.out_stat[self.dat_key].keys()])

        for method in is_float:
            assert isinstance(self.out_stat[self.dat_key][method], float) \
                or np.floating(self.out_stat[self.dat_key][method]), \
                'unexpected statistic value for {:s}: {:}'.format(
                    method, self.out_stat[self.dat_key][method])

        for method in not_float:
            assert len(self.out_stat[self.dat_key][method]) > 0
        return

    @pytest.mark.parametrize("input_ind,input_val,err_tst",
                             [(0, None, "must provide Dataset of paired"),
                              (1, [], "must provide equal number of instr"),
                              (1, ["Ti"], "unknown instrument data value"),
                              (2, ["Ti"], "unknown model data value"),
                              (3, ["nada"], "unknown statistical method(s)")])
    def test_compare_model_and_inst_input_failure(self, input_ind, input_val,
                                                  err_tst):
        """Test raises ValueError for badly shaped input.

        Parameters
        ----------
        input_ind : int
            Index of input to replace in `self.input_args`
        input_val :
            Data of any type to use as input, replacing the list value at
            `input_ind` in `self.input_args`
        err_tst : str
            Expected error message

        """
        self.input_args[input_ind] = input_val

        with pytest.raises(ValueError) as verr:
            compare.compare_model_and_inst(*self.input_args)

        assert str(verr).find(err_tst) >= 0

        return

    def test_compare_model_and_inst_warning_partial_inst_input(self, caplog):
        """Test raises logging warning for some unknown inputs varibles."""
        self.input_args[1].append('not_a_variable')
        self.input_args[2].append(self.mod_key)

        with caplog.at_level(logging.WARNING, logger='pysatModels'):
            self.out_stat, self.out_units = compare.compare_model_and_inst(
                *self.input_args)

        # Evaluate the logging output
        assert len(caplog.records) >= 1
        assert caplog.records[-1].message.find('removed 1 of the provided') >= 0
        assert caplog.records[-1].levelname == "WARNING"

        # Evalute the comparison output after removing the bad inputs
        self.input_args[1].pop()
        self.input_args[2].pop()
        self.eval_comp_output()

        return

    @pytest.mark.parametrize('input_ind', [1, 2, 3])
    def test_compare_model_and_inst_warning_cast_input(self, input_ind):
        """Test raises logging warning for some unknown inputs varibles.

        Parameters
        ----------
        input_ind : int
            Index of input to replace in `self.input_args`

        """
        # Save the original input and set the function input to be an element
        # instead of a list
        orig_arg = list(self.input_args[input_ind])
        self.input_args[input_ind] = orig_arg[0]

        # Run the code
        self.out_stat, self.out_units = compare.compare_model_and_inst(
            *self.input_args)

        # Evalute the comparison output after removing the bad inputs
        self.input_args[input_ind] = orig_arg
        self.eval_comp_output()

        return

    def test_compare_model_and_inst_warning_partial_mod_input(self, caplog):
        """Test raises logging warning for some unknown inputs varibles."""
        self.input_args[1].append(self.dat_key)
        self.input_args[2].append('not_a_variable')

        with caplog.at_level(logging.WARNING, logger='pysatModels'):
            self.out_stat, self.out_units = compare.compare_model_and_inst(
                *self.input_args)

        # Evaluate the logging output
        assert len(caplog.records) >= 1
        assert caplog.records[-1].message.find('removed 1 of the provided') >= 0
        assert caplog.records[-1].levelname == "WARNING"

        # Evalute the comparison output after removing the bad inputs
        self.input_args[1].pop()
        self.input_args[2].pop()
        self.eval_comp_output()

        return

    def test_compare_model_and_inst_data_success(self):
        """Test the successful comparison for all possible statistics."""

        # Get the output data
        self.out_stat, self.out_units = compare.compare_model_and_inst(
            *self.input_args)
        self.eval_comp_output()

        return
