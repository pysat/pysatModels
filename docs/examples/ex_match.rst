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

In the example below, we load a series of TIEGCM files that were created every
15 minutes and are stored in a specified directory.  We pair them with CINDI
C/NOFS ion drift data that was previously downloaded using pysat.

::


   import datetime as dt
   import pysat
   import pysatModels as ps_mod

   stime = dt.datetime(2010, 1, 1)
   etime = dt.datetime(2010, 2, 1)
   tinc = dt.timedelta(minutes=15)
   input_kwargs = dict()

   # Initialize the model input information, including the keyword arguements
   # needed to load the model Instrument into an xarray Dataset using
   # the match.load_model_xarray routine
   filename = 'sample_experiment/tiegcm_%Y%j.%H%M.nc'
   tiegcm = pysat.Instrument(platform='ucar', name='tiegcm')
   input_kwargs["model_load_kwargs"] = {'model_inst': tiegcm,'filename': filename}

   # Initialize the CINDI C/NOFS input, and ensure the best model interpolation
   # by extracting the clean model data after matching
   cindi = pysat.Instrument(platform='cnofs', name='ivm', clean_level='none')
   cindi.download(stime, etime) # Skip this if you already have the data
   input_kwargs["inst_clean_rout"] = pysat.instruments.cnofs_ivm.clean
   input_kwargs["inst_download_kwargs"] = {"skip_download": True}

   # Many keyword arguements are required, as they provide information on
   # the names of coordinates needed to extract data
   input_kwargs["inst_lon_name"] = "glon"
   input_kwargs["mod_lon_name"] = "longitude"
   input_kwargs["inst_name"] = ["glon", "glat", "altitude"]
   input_kwargs["mod_name"] = ["longitude", "latitude", "altitude"]
   input_kwargs["mod_units"] = ["deg", "deg", "km"]
   input_kwargs["mod_datetime_name"] = "time"
   input_kwargs["mod_time_name"] = "time"

   # Now we are ready to run using the defaults!
   matched_inst = ps_mod.utils.match.collect_inst_model_pairs(stime, etime, tinc, cindi, **input_kwargs)


This returns a pysat Instrument object with the CINDI C/NOFS data and TIEGCM
data at the same times and locations.  The CINDI data has the same names as
the normal Instrument, and the TIEGCM data has the same name as the TIEGCM
data, but with ``model_`` as a prefix to prevent confusion.  You can change
this prefix using the ``model_label`` keyword argument, allowing multiple
models to be matched to the same observational data set.
   
   
::

   # Using the results from the prior example
   print([cc for cc in matched_inst.data.columns if cc.find("model_") == 0])
