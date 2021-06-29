#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2019, pysat team
# Full license can be found in License.md
# -----------------------------------------------------------------------------

import numpy as np
import os
import pytest
import xarray as xr

import pysat

from pysatModels.utils import convert


def eval_xarray_output(inst, xdata):
    """ Evaluate the convertion from pysat.Instrument to xarray.Dataset

    Parameters
    ----------
    inst : pysat.Instrument
        pysat Instrument object with or without loaded data
    xdata : xarray.Dataset or NoneType
        Data extracted from pysat.Instrument, with metadata if possible

    Raises
    ------
    AssertionError
        If data type is incorrect, if data variables don't match, or if
        metadata variables don't match

    """
    if inst.empty:
        # Test the output data type
        assert xdata is None
    else:
        # Test the output data type
        assert isinstance(xdata, xr.core.dataset.Dataset)

        # Test that the time index is a coordinate
        assert inst.index.name in xdata.coords

        # Test the data variables
        for kk in inst.variables:
            assert kk in xdata.data_vars or kk in xdata.coords

            # Test that the metadata was attached if this isn't a coordinate
            if kk not in xdata.coords:
                for mm in inst.meta.data.keys():
                    try:
                        assert xdata.data_vars[kk].attrs[mm] == inst.meta[kk,
                                                                          mm]
                    except AssertionError:
                        if np.isnan(inst.meta[kk, mm]):
                            assert np.isnan(xdata.data_vars[kk].attrs[mm])
                        elif np.isinfinite(inst.meta[kk, mm]):
                            assert np.isinfinite(xdata.data_vars[kk].attrs[mm])
                        else:
                            raise ValueError('Unknown metadata type')
    return


class TestUtilsConvertLoadModelXarray():
    """ Unit tests for utils.convert.load_model_xarray
    """
    def setup(self):
        """ Runs before every method to create a clean testing setup
        """
        self.ftime = pysat.instruments.pysat_testing_xarray._test_dates['']['']
        self.filename = "%Y-%m-%d.nofile"
        self.model_kwargs = {'platform': str('pysat'),
                             'name': str('testing_xarray'),
                             'num_samples': 12,
                             'clean_level': 'clean'}
        self.model_inst = None
        self.xout = None
        self.temp_file = 'None'

    def teardown(self):
        if os.path.isfile(self.temp_file):
            os.remove(self.temp_file)

        del self.ftime, self.model_kwargs, self.xout, self.filename
        del self.temp_file, self.model_inst

    def test_no_inst(self):
        """ Test failure when no instrument object is provided
        """
        with pytest.raises(ValueError) as verr:
            convert.load_model_xarray(self.ftime)

        assert verr.value.args[0].find("must provide a pysat.Instrument") >= 0
        return

    @pytest.mark.parametrize("fname", [(None), ('filename')])
    def test_load_filename(self, fname):
        """ Test success when loading through different filename options
        """
        if fname is not None:
            if hasattr(self, fname):
                fname = getattr(self, fname)
            self.temp_file = self.ftime.strftime(fname)

            # Create a temporary file
            with open(self.temp_file, 'w') as fout:
                fout.write('')

        # Load the test instrument data
        self.model_inst = pysat.Instrument(**self.model_kwargs)
        self.xout = convert.load_model_xarray(self.ftime, self.model_inst,
                                              filename=fname)

        # Evaluate the returned data
        eval_xarray_output(self.model_inst, self.xout)
        return

    @pytest.mark.parametrize("mkey, mval", [("name", "testing"),
                                            (None, None)])
    def test_load_inst(self, mkey, mval):
        """ Test success when loading different types of pysat Instruments
        """
        if mkey in self.model_kwargs.keys():
            self.model_kwargs[mkey] = mval
        self.model_inst = pysat.Instrument(**self.model_kwargs)
        self.xout = convert.load_model_xarray(self.ftime, self.model_inst)

        # Evaluate the returned data
        eval_xarray_output(self.model_inst, self.xout)
        return


class TestUtilsConvertPysatXarray():
    """ Unit tests for utils.convert.convert_pydat_to_xarray
    """
    def setup(self):
        """ Runs before every method to create a clean testing setup
        """
        self.ref_time = pysat.instruments.pysat_testing._test_dates['']['']

    def teardown(self):
        del self.ref_time

    @pytest.mark.parametrize("name, load", [('testing', False),
                                            ('testing', True),
                                            ('testmodel', True)])
    def test_convert_pysat_to_xarray(self, name, load):
        """ Test success when converting pysat data to an xarray.Dataset
        """
        # Initialize the pysat Instrument, loading if desired. Also ensure
        # load worked correctly
        inst = pysat.Instrument('pysat', name)

        if load:
            inst.load(date=self.ref_time)
            assert not inst.empty, "Instrument should not be empty"
        else:
            assert inst.empty, "Instrument should not have data"

        # Get the xarray data
        xdata = convert.convert_pysat_to_xarray(inst)

        # Evaluate the returned data
        eval_xarray_output(inst, xdata)

        return
