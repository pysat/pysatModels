Pair Modelled and Observed Data
===============================

One common analytical need is to obtain a combined data set of modelled and
observed observations at the same times and locations.  This can be done using
the ``collect_inst_model_pairs`` routine in
`utils.match <../utils.html#module-pysatModels.utils.match>`_.  This routine
takes a date range as input, and then extracts modelled observations at the
specified instrument location.  However, it does not interpolate in time.
Details about the interpolation of modelled data onto the instrument location
can be found in the
`utils.extract <../utils.html#module-pysatModels.utils.extract>`_ routine:
``extract_modelled_observations``.

In the example below, we load SAMI2 files that were created for testing and
pair them with DMSP IVM data from F18.  This example uses the external modules:

- pysat
- pysatMadrigal

::


   import datetime as dt
   import pysat
   import pysatMadrigal
   import pysatModels as ps_mod

   stime = dt.datetime(2019, 1, 1)
   tinc = dt.timedelta(days=1) # Required input is not used in this example
   input_kwargs = dict()

   # Initialize the model input information, including the keyword arguements
   # needed to load the model Instrument into an xarray Dataset using
   # the match.load_model_xarray routine
   sami2 = pysat.Instrument(inst_module=ps_mod.models.sami2py_sami2, tag='test')
   sami2.download(stime, stime) # Optional, only needs to be done once
   filename = sami2.files.files[stime]
   input_kwargs["model_load_kwargs"] = {'model_inst': sami2, 'filename': filename}

   # Initialize the DMSP IVM input, and ensure the best model interpolation
   # by extracting the clean model data after matching.
   dmsp_f18 = pysat.Instrument(inst_module=pysatMadrigal.instruments.dmsp_ivm, sat_id='f18', clean_level='none')
   dmsp_f18.download(stime, stime) # Skip this if you already have the data
   input_kwargs["inst_clean_rout"] = pysatMadrigal.instruments.dmsp_ivm.clean
   input_kwargs["inst_download_kwargs"] = {"skip_download": True}

   # Many keyword arguements are required, as they provide information on
   # the names of coordinates needed to extract data
   input_kwargs["inst_lon_name"] = "glon"
   input_kwargs["mod_lon_name"] = "glon"
   input_kwargs["inst_name"] = ["glon", "gdlat", "gdalt"]
   input_kwargs["mod_name"] = ["glon", "glat", "zalt"]
   input_kwargs["mod_units"] = ["deg", "deg", "km"]
   input_kwargs["mod_datetime_name"] = "time"
   input_kwargs["mod_time_name"] = "time"

   # Now we are ready to run using the defaults!
   matched_inst = ps_mod.utils.match.collect_inst_model_pairs(stime, stime+tinc, tinc, dmsp_f18, **input_kwargs)


This returns a pysat Instrument object with the DMSP IVM data and SAMI2
data at the same times and locations.  The DMSP IVM data has the same names as
the normal Instrument, and the SAMI2 data has the same name as the SAMI2
data, but with ``model_`` as a prefix to prevent confusion.  You can change
this prefix using the ``model_label`` keyword argument, allowing multiple
models to be matched to the same observational data set.
   
   
::

   # Using the results from the prior example
   print([cc for cc in matched_inst.data.columns if cc.find("model_") == 0])
