# Author: Angeline Burrell, NRL, 2019
"""Unit tests for `pysatModels.utils.compare`."""

import numpy as np
import pytest
import xarray as xr

import pysatModels.utils.compare as compare


class TestUtilsCompare():
    """Unit tests for utils.compare."""

    def setup(self):
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
        return

    def teardown(self):
        """Clean up the unit test environment after each method."""

        del self.dat_key, self.mod_key, self.units, self.pairs, self.input_args
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

    def test_compare_model_and_inst_data_success(self):
        """Test the successful comparison for all possible statistics."""

        # Get the output data
        out_stat, out_units = compare.compare_model_and_inst(*self.input_args)

        # Evaluate the output data dictionary structure
        assert out_stat.keys() == out_units.keys()
        assert self.dat_key in out_stat.keys()
        assert self.mod_key not in out_stat.keys()
        assert len(out_stat.keys()) == 1

        # Evaluate the unit output
        assert out_units[self.dat_key] == self.units

        # Evaluate the number of returned statistics (all requested)
        is_float = ['bias', 'meanPercentageError', 'medianLogAccuracy',
                    'symmetricSignedBias', 'meanSquaredError', 'RMSE',
                    'meanAbsError', 'medAbsError', 'nRMSE', 'MASE',
                    'medSymAccuracy', 'meanAPE', 'medAPE']
        not_float = ['logAccuracy', 'scaledError', 'forecastError', 'percError',
                     'absPercError']

        assert (len(out_stat[self.dat_key].keys())
                == len(is_float) + len(not_float)), \
            "unexpected number of statistics returned: {:}".format(
                [kk for kk in out_stat[self.dat_key].keys()])

        for method in is_float:
            assert isinstance(out_stat[self.dat_key][method], float) \
                or np.floating(out_stat[self.dat_key][method]), \
                'unexpected statistic value for {:s}: {:}'.format(
                    method, out_stat[self.dat_key][method])

        for method in not_float:
            assert len(out_stat[self.dat_key][method]) > 0

        return
