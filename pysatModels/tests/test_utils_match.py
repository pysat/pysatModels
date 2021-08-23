# Author: Angeline Burrell, NRL, 2019
"""Unit tests for `pysatModels.utils.match`."""

import datetime as dt
from io import StringIO
import logging
import numpy as np
import os
import pytest
import xarray as xr

import pysat

import pysatModels as ps_mod
import pysatModels.utils.match as match


class TestUtilsMatchLoadModelXarray(object):
    """Unit tests for utils.match.load_model_xarray."""

    def setup(self):
        """Set up the unit test environment for each method."""

        self.ftime = dt.datetime(2009, 1, 1)
        self.filename = "%Y-%m-%d.nofile"
        self.model_kwargs = {'platform': str('pysat'),
                             'name': str('testing_xarray'),
                             'num_samples': 12,
                             'clean_level': 'clean'}
        self.model_inst = None
        self.xout = None
        self.temp_file = 'None'
        return

    def teardown(self):
        """Clean up the unit test environment after each method."""

        if os.path.isfile(self.temp_file):
            os.remove(self.temp_file)

        del self.ftime, self.model_kwargs, self.xout, self.filename
        del self.temp_file, self.model_inst
        return

    def test_no_inst(self):
        """Test failure when no instrument object is provided."""

        with pytest.raises(ValueError) as verr:
            match.load_model_xarray(self.ftime)

        assert verr.value.args[0].find("must provide a pysat.Instrument") >= 0
        return

    @pytest.mark.parametrize("fname", [(None), ('filename')])
    def test_load_filename(self, fname):
        """Test success when loading through different filename options."""

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
        return

    def test_load_pandas_inst(self):
        """Test success when loading a panads pysat Instrument."""

        self.model_kwargs["name"] = "testing"
        self.model_inst = pysat.Instrument(**self.model_kwargs)
        self.xout = match.load_model_xarray(self.ftime, self.model_inst)

        # Test the output
        assert isinstance(self.xout, xr.core.dataset.Dataset)

        for kk in self.model_inst.data.columns:
            assert kk in self.xout.data_vars

        assert self.model_inst.index.name in self.xout.coords
        return

    def test_load_empty_inst(self):
        """Test return value of None with empty instrument load."""

        self.model_inst = pysat.Instrument(**self.model_kwargs)
        self.ftime = self.model_inst.files.files.index[0] - dt.timedelta(
            days=1)
        self.xout = match.load_model_xarray(self.ftime, self.model_inst)

        assert self.xout is None
        return


class TestUtilsMatchCollectInstModPairs(object):
    """Unit tests for utils.match.collect_inst_model_pairs."""

    def setup(self):
        """Set up the unit test environment for each method."""

        self.inst = pysat.Instrument(platform='pysat', name='testing')
        self.stime = pysat.instruments.pysat_testing._test_dates['']['']
        self.inst.load(date=self.stime)
        self.input_args = [self.stime, self.stime + dt.timedelta(days=1),
                           dt.timedelta(days=1), self.inst]
        self.ref_col = 'dummy1'
        self.model = pysat.Instrument(platform='pysat', name='testmodel',
                                      num_samples=10)
        self.required_kwargs = {"model_load_kwargs":
                                {"model_inst": self.model},
                                "inst_clean_rout": lambda x: True,
                                "inst_lon_name": "longitude",
                                "mod_lon_name": "longitude",
                                "inst_name": ["longitude", "latitude"],
                                "mod_name": ["longitude", "latitude"],
                                "mod_datetime_name": "time",
                                "mod_time_name": "time",
                                "mod_units": ["deg", "deg"]}
        self.log_capture = StringIO()
        ps_mod.logger.addHandler(logging.StreamHandler(self.log_capture))
        ps_mod.logger.setLevel(logging.INFO)
        self.out = None
        return

    def teardown(self):
        """Clean up the unit test environment after each method."""

        del self.input_args, self.required_kwargs, self.inst, self.model
        del self.out, self.log_capture, self.stime
        return

    @pytest.mark.parametrize("mkey,mout",
                             [("verr", None), ("ierr", None),
                              ("other", "Unacceptable model load error")])
    def test_model_load_failure(self, mkey, mout):
        """Test for expected failure when unable to load model data."""

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
            with pytest.raises(TypeError, match=mout):
                match.collect_inst_model_pairs(*self.input_args,
                                               **self.required_kwargs)
        return

    @pytest.mark.parametrize("del_key,err_msg",
                             [("inst_lon_name",
                               "Need longitude name for inst"),
                              ("mod_lon_name",
                               "Need longitude name for model"),
                              ("inst_name",
                               "Must provide instrument location"),
                              ("mod_name", "Must provide model location attr"),
                              ("mod_units", "Must provide model units as a "),
                              ("mod_datetime_name",
                               "Need datetime coordinate"),
                              ("mod_time_name", "Need time coordinate"),
                              ("inst_clean_rout", "Need routine to clean")])
    def test_input_failure(self, del_key, err_msg):
        """Test for expected failure when missing requried input from kwargs."""

        del self.required_kwargs[del_key]

        with pytest.raises(ValueError, match=err_msg):
            match.collect_inst_model_pairs(*self.input_args,
                                           **self.required_kwargs)

    @pytest.mark.parametrize("cng_key,bad_val,err_msg",
                             [("mod_lon_name", "glon",
                               "unknown name for model longitude"),
                              ("mod_datetime_name", "dt",
                               "unknown model name for datetime"), ])
    def test_bad_input(self, cng_key, bad_val, err_msg):
        """Test for expected failure with bad input."""

        self.required_kwargs[cng_key] = bad_val

        with pytest.raises(ValueError, match=err_msg):
            match.collect_inst_model_pairs(*self.input_args,
                                           **self.required_kwargs)
        return

    def test_bad_time(self):
        """Test the match routine the times prevent any data from loading."""

        self.input_args[1] = self.input_args[0]

        assert match.collect_inst_model_pairs(*self.input_args,
                                              **self.required_kwargs) is None
        return

    @pytest.mark.parametrize("tinc_val, einc, num",
                             [(dt.timedelta(days=1), dt.timedelta(days=1), 3),
                              (dt.timedelta(days=2), dt.timedelta(days=1), 3),
                              (dt.timedelta(days=1), dt.timedelta(days=2), 6)])
    def test_tinc_success(self, tinc_val, einc, num):
        """Test the match success with different time increments."""

        self.input_args[1] = self.stime + einc
        self.input_args[2] = tinc_val
        self.required_kwargs['model_label'] = 'tmodel'
        self.required_kwargs['sel_name'] = [self.ref_col]
        self.ref_col = '_'.join([self.required_kwargs['model_label'],
                                 self.ref_col])

        self.out = match.collect_inst_model_pairs(*self.input_args,
                                                  **self.required_kwargs)

        assert self.ref_col in [kk for kk in self.out.variables]
        assert len(self.out[self.ref_col]) == num
        assert np.all(np.isfinite(self.out[self.ref_col].values))
        assert self.out.index[0].date() == self.stime.date()
        assert self.out.index[-1] < self.stime + einc
        return

    @pytest.mark.parametrize("lin, lout, test_out",
                             [([-179.0, 179.0], [-180.0, 180.0], 3),
                              ([0.5, 359.0], [0.0, 360.0], 3),
                              ([-1.0, 210.0], None,
                               'unexpected longitude range')])
    def test_lon_output(self, lin, lout, test_out):
        """Test the match handling with different longitude range input."""

        def lon_model_load(ftime, model_inst=None, filename=None):
            mdata = match.load_model_xarray(ftime, model_inst, filename)

            if mdata is not None:
                if lout is None:
                    lin.append(mdata.dims['time'])
                    glon = {'glon': (('time'), np.linspace(*lin))}
                    mdata = mdata.assign(glon)
                else:
                    lin.append(mdata.dims['longitude'])
                    mdata.coords['longitude'] = np.linspace(*lin)
            return mdata

        self.required_kwargs['model_load_rout'] = lon_model_load
        self.required_kwargs['model_label'] = 'tmodel'
        self.required_kwargs['sel_name'] = [self.ref_col]
        self.ref_col = '{:s}_{:s}'.format(self.required_kwargs['model_label'],
                                          self.ref_col)

        if lout is None:
            hold_lon = self.required_kwargs['mod_lon_name']
            self.required_kwargs['mod_lon_name'] = 'glon'
            self.required_kwargs['mod_name'][0] = 'glon'

            with pytest.raises(ValueError, match=test_out):
                match.collect_inst_model_pairs(*self.input_args,
                                               **self.required_kwargs)

            self.required_kwargs['mod_lon_name'] = hold_lon
            self.required_kwargs['mod_name'][0] = hold_lon
        else:
            self.out = match.collect_inst_model_pairs(*self.input_args,
                                                      **self.required_kwargs)

            assert self.ref_col in [kk for kk in self.out.variables]
            assert len(self.out[self.ref_col]) == test_out
            assert np.all(np.isfinite(self.out[self.ref_col].values))
            assert self.out['longitude'].min() >= lout[0]
            assert self.out['longitude'].max() <= lout[1]
        return

    def test_success_skip_download(self):
        """Test the match success with skip_download key."""

        self.required_kwargs['inst_download_kwargs'] = {'skip_download': True}
        self.required_kwargs['model_label'] = 'tmodel'
        self.required_kwargs['sel_name'] = [self.ref_col]
        self.ref_col = '{:s}_{:s}'.format(self.required_kwargs['model_label'],
                                          self.ref_col)

        self.out = match.collect_inst_model_pairs(*self.input_args,
                                                  **self.required_kwargs)

        assert self.ref_col in [kk for kk in self.out.variables]
        assert len(self.out[self.ref_col]) > 0
        return

    def test__inst_download_missing(self):
        """Test the download data loop, which will fail to download anything."""

        self.input_args[0] = self.inst.files.files.index[0] - dt.timedelta(
            days=1)
        self.input_args[1] = self.inst.files.files.index[0]

        self.out = match.collect_inst_model_pairs(*self.input_args,
                                                  **self.required_kwargs)

        assert self.out is None
        return
