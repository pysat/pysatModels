#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2022, pysat development team
# Full license can be found in License.md
# -----------------------------------------------------------------------------
"""Unit tests for `pysatModels.utils.testing`."""

import numpy as np
from packaging import version as pack_version
import pytest

import pysat
from pysat.instruments import pysat_testmodel

from pysatModels.utils.testing import compare_mod_name_coordinates


class TestUtilsCompareModName(object):
    """Unit tests for `compare_mod_name_coordinates`."""

    def setup(self):
        """Set up the unit test environment for each method."""

        self.model = pysat.Instrument(inst_module=pysat_testmodel, tag='')

        # Load the data in the instruments
        load_kwargs = {'date': pysat_testmodel._test_dates['']['']}
        if(pack_version.Version(pysat.__version__)
           > pack_version.Version('3.0.1')):
            load_kwargs['use_header'] = True

        self.model.load(**load_kwargs)

        return

    def teardown(self):
        """Clean up the unit test environment after each method."""

        del self.model
        return

    def test_failure_non_timelike_time(self):
        """Test the failure when first dimension not timelike."""

        # Replace time values with integers
        self.model['time'] = np.arange(0, len(self.model['time']))

        # Run, raising a value error
        with pytest.raises(ValueError) as err:
            compare_mod_name_coordinates(self.model['dummy1'],
                                         ['latitude', 'longitude'])

        assert str(err).find('does not appear to be a time dimension') >= 0

        return

    def test_failure_too_many(self):
        """Test the failure when `mod_name` too long."""

        # Run, raising a value error
        with pytest.raises(ValueError) as err:
            compare_mod_name_coordinates(self.model['dummy1'],
                                         ['latitude', 'longitude', 'altitude'])

        assert str(err).find('Provided too many ') >= 0

        return

    def test_compare_model_name_coordinates_success(self):
        """Ensure `compare_model_name_coordinates` works for proper inputs."""

        compare_mod_name_coordinates(self.model['dummy1'],
                                     ['latitude', 'longitude'])
        return

    @pytest.mark.parametrize("var,coords,msg,flag", [('dummy1',
                                                      ['longitude', 'latitude'],
                                                      'not in same order as',
                                                      False),
                                                     ('dummy1',
                                                      ['wrong', 'latitude'],
                                                      'not all within variable',
                                                      False),
                                                     ('dummy1',
                                                      ['latitude', 'longitude'],
                                                      'appear to be a time',
                                                      True)])
    def test_compare_model_name_coordinates_failure(self, var, coords, msg,
                                                    flag):
        """Ensure `compare_mod_name_coordinates` works for proper inputs.

        Parameters
        ----------
        var : str
            Variable to test.
        coords : list
            List of coordinate variable names.
        msg : str
            Failure message to test for.
        flag : bool
            If True, replace time index with integers.

        """

        if flag:
            self.model['time'] = np.arange(0, len(self.model['time']))

        with pytest.raises(ValueError) as err:
            compare_mod_name_coordinates(self.model[var], coords)

        assert str(err).find(msg) >= 0

        return


@pytest.mark.skipif(pack_version.Version(pysat.__version__)
                    <= pack_version.Version('3.0.1'),
                    reason=''.join(('Requires test model in pysat ',
                                    ' v3.0.2 or later.')))
class TestUtilsCompareModNamePressure(TestUtilsCompareModName):
    """Unit tests for `compare_mod_name_coordinates`."""

    def setup(self):
        """Set up the unit test environment for each method."""

        self.model = pysat.Instrument(inst_module=pysat_testmodel,
                                      tag='pressure_levels')

        # Load the data in the instruments
        load_kwargs = {'date': pysat_testmodel._test_dates['']['']}
        if(pack_version.Version(pysat.__version__)
           > pack_version.Version('3.0.1')):
            load_kwargs['use_header'] = True

        self.model.load(**load_kwargs)

        return
