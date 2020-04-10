# Author: Angeline Burrell, NRL, 2019

from __future__ import absolute_import, unicode_literals

import datetime as dt
import logging
from io import StringIO
import os
import pytest
import xarray as xr

import pysat

import pysatModels as ps_mod
import pysatModels.utils.match as match


class TestUtilsMatchLoadModelXarray:
    """ Unit tests for utils.match.load_model_xarray
    """
    def setup(self):
        """ Runs before every method to create a clean testing setup
        """
        self.ftime = dt.datetime(2009, 1, 1)
        self.filename = "%Y-%m-%d.nofile"
        self.model_kwargs = {'platform': str('pysat'),
                             'name': str('testing_xarray'),
                             'sat_id': '12',
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
            match.load_model_xarray(self.ftime)

        assert verr.value.args[0].find("must provide a pysat.Instrument") >= 0

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
        self.xout = match.load_model_xarray(self.ftime, self.model_inst,
                                            filename=fname)

        # Test the output
        assert isinstance(self.xout, xr.core.dataset.Dataset)

        for kk in self.model_inst.data.data_vars:
            assert kk in self.xout.data_vars

        assert self.model_inst.index.name in self.xout.coords

    def test_load_pandas_inst(self):
        """ Test success when loading a panads pysat Instrument
        """
        self.model_kwargs["name"] = "testing"
        self.model_inst = pysat.Instrument(**self.model_kwargs)
        self.xout = match.load_model_xarray(self.ftime, self.model_inst)

        # Test the output
        assert isinstance(self.xout, xr.core.dataset.Dataset)

        for kk in self.model_inst.data.columns:
            assert kk in self.xout.data_vars

        assert self.model_inst.index.name in self.xout.coords

    def test_load_empty_inst(self):
        """ Test return value of None with empty instrument load
        """
        self.model_inst = pysat.Instrument(**self.model_kwargs)
        self.ftime = self.model_inst.files.files.index[0] - dt.timedelta(days=1)
        self.xout = match.load_model_xarray(self.ftime, self.model_inst)

        assert self.xout is None


class TestUtilsMatchCollectInstModPairs:
    """ Unit tests for utils.match.collect_inst_model_pairs
    """
    def setup(self):
        """Runs before every method to create a clean testing setup
        """
        self.inst = pysat.Instrument(platform=str('pysat'),
                                     name=str('testing'),
                                     clean_level='clean')
        self.inst.load(yr=2009, doy=1)
        self.input_args = [dt.datetime(2009, 1, 1), dt.datetime(2009, 1, 2),
                           dt.timedelta(days=1), self.inst]
        self.model = pysat.Instrument(platform=str('pysat'),
                                      name=str('testing2d_xarray'), sat_id='10',
                                      clean_level='clean')
        self.required_kwargs = {"model_load_kwargs":
                                {"model_inst": self.model},
                                "inst_clean_rout": lambda x: True,
                                "inst_lon_name": "longitude",
                                "mod_lon_name": "longitude",
                                "inst_name": ["longitude", "latitude", "slt"],
                                "mod_name": ["longitude", "latitude", "slt"],
                                "mod_datetime_name": "time",
                                "mod_time_name": "time",
                                "mod_units": ["deg", "deg", "h"]}
        self.log_capture = StringIO()
        ps_mod.logger.addHandler(logging.StreamHandler(self.log_capture))
        ps_mod.logger.setLevel(logging.INFO)
        self.out = None

    def teardown(self):
        """Runs after every method to clean up previous testing
        """
        del self.input_args, self.required_kwargs, self.inst, self.model
        del self.out, self.log_capture

    @pytest.mark.parametrize("mkey,mout",
                             [("verr", None), ("ierr", None),
                              ("other", "Unacceptable model load error")])
    def test_model_load_failure(self, mkey, mout):
        """ Test for expected failure when unable to load model data
        """
        def model_load_rout(stime, verr=False, ierr=False):
            if verr:
                raise ValueError('Acceptable model load error')
            elif ierr:
                raise IOError('Acceptable model load error')
            else:
                raise TypeError('Unacceptable model load error')

        self.required_kwargs['model_load_rout'] = model_load_rout

        if mkey in ['verr', 'ierr']:
            self.required_kwargs['model_load_kwargs'] = {mkey: True}
            self.out = match.collect_inst_model_pairs(*self.input_args,
                                                      **self.required_kwargs)
            lout = self.log_capture.getvalue()

            assert self.out is mout
            assert lout.find('unable to load model data at') >= 0
        else:
            self.required_kwargs['model_load_kwargs'] = {}
            with pytest.raises(TypeError, match=mout) as terr:
                match.collect_inst_model_pairs(*self.input_args,
                                               **self.required_kwargs)

    @pytest.mark.parametrize("del_key,err_msg",
                             [("inst_lon_name", "Need longitude name for inst"),
                              ("mod_lon_name", "Need longitude name for model"),
                              ("inst_name", "Must provide instrument location"),
                              ("mod_name", "Must provide the same number"),
                              ("mod_units", "Must provide units for each "),
                              ("mod_datetime_name", "Need datetime coordinate"),
                              ("mod_time_name", "Need time coordinate"),
                              ("inst_clean_rout", "Need routine to clean")])
    def test_input_failure(self, del_key, err_msg):
        """ Test for expected failure when missing necessary input from kwargs
        """
        del self.required_kwargs[del_key]

        with pytest.raises(ValueError, match=err_msg) as verr:
            match.collect_inst_model_pairs(*self.input_args,
                                           **self.required_kwargs)

    @pytest.mark.parametrize("cng_key,bad_val,err_msg",
                             [("mod_lon_name", "glon",
                               "unknown name for model longitude"),
                              ("mod_datetime_name", "dt",
                               "unknown model name for datetime"), ])
    def test_bad_input(self, cng_key, bad_val, err_msg):
        """ Test for expected failure with bad input
        """
        self.required_kwargs[cng_key] = bad_val

        with pytest.raises(ValueError, match=err_msg) as verr:
            match.collect_inst_model_pairs(*self.input_args,
                                           **self.required_kwargs)

    def test_bad_time(self):
        """ Test the match routine the times prevent any data from loading
        """
        self.input_args[1] = self.input_args[0]

        assert match.collect_inst_model_pairs(*self.input_args,
                                              **self.required_kwargs) is None

    @pytest.mark.parametrize("tinc_val", [dt.timedelta(days=1),
                                          dt.timedelta(days=2)])
    def test_tinc_success(self, tinc_val):
        """ Test the match success with different time increments
        """
        self.input_args[2] = tinc_val

        self.out = match.collect_inst_model_pairs(*self.input_args,
                                                  **self.required_kwargs)

        assert isinstance(self.out.data, xr.Dataset)

    def test_1D_model_success(self):
        """ Test the match success with one dimensional model data
        """
        self.model = pysat.Instrument(platform=str('pysat'),
                                      name=str('testing_xarray'), sat_id='10',
                                      clean_level='clean')
        self.required_kwargs["model_load_kwargs"] = {"model_inst": self.model}

        self.out = match.collect_inst_model_pairs(*self.input_args,
                                                  **self.required_kwargs)

        assert isinstance(self.out.data, xr.Dataset)

    def test_tinc_success_skip_download(self):
        """ Test the match success with skip_download key
        """
        self.required_kwargs['inst_download_kwargs'] = {'skip_download': True}

        self.out = match.collect_inst_model_pairs(*self.input_args,
                                                  **self.required_kwargs)

        assert isinstance(self.out.data, xr.Dataset)

    def test_tinc_inst_download_missing(self):
        """ Test the download data loop, which will fail to download anything
        """
        self.input_args[0] = self.inst.files.files.index[0]-dt.timedelta(days=1)
        self.input_args[1] = self.inst.files.files.index[0]

        self.out = match.collect_inst_model_pairs(*self.input_args,
                                                  **self.required_kwargs)

        assert self.out is None
