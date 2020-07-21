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
pair them with DMSP IVM data from F18.  This is a regional example, rather than
a global example (because SAMI2 runs along a single magnetic meridian).  This
makes the input keyword arguments different and requires a custom function to
only load the satellite data along the desired magnetic meridian.

This example uses the external modules:

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

   # Define the custom function needed to limit the DMSP data
   lon_min = sami2['glon'].values.min()
   lon_max = sami2['glon'].values.max()
   input_kwargs["lon_pos"] = 0 # Adjust lon range before other custom functions
   def sami2_meridian(inst, lon_kwarg, min_long, max_long): 
       good_mask = (inst[lon_kwarg] >= min_long) & (inst[lon_kwarg] <= max_long)
       bad_index = inst.index[~good_mask]
       good_data = inst.data.drop(index=bad_index)
       inst.data = good_data
       return

   # Initialize the DMSP IVM input, and ensure the best model interpolation
   # by extracting the clean model data after matching.
   dmsp_f18 = pysat.Instrument(inst_module=pysatMadrigal.instruments.dmsp_ivm, sat_id='f18', clean_level='none')
   dmsp_f18.download(stime, stime) # Skip this if you already have the data
   dmsp_f18.custom.attach(sami2_meridian, kind='modify', at_pos='end', args=['glon', lon_min, lon_max])
   input_kwargs["inst_clean_rout"] = pysatMadrigal.instruments.dmsp_ivm.clean
   input_kwargs["inst_download_kwargs"] = {"skip_download": True}

   # Many keyword arguements are required, as they provide information on
   # the names of coordinates needed to extract data
   input_kwargs["inst_lon_name"] = "glon"
   input_kwargs["mod_lon_name"] = "glon"
   input_kwargs["inst_name"] = ["gdlat", "gdalt"] # Not considering longitude
   input_kwargs["mod_name"] = ["glat", "zalt"]    # in the data matching
   input_kwargs["mod_units"] = ["deg", "km"]
   input_kwargs["mod_datetime_name"] = "time"
   input_kwargs["mod_time_name"] = "time"

   # Now we are ready to run using the defaults!
   matched_inst = ps_mod.utils.match.collect_inst_model_pairs(stime, stime+tinc, tinc, dmsp_f18, **input_kwargs)


This returns a pysat Instrument object with the DMSP IVM data and SAMI2
data at the same times, latitudes, and altitudes along the SAMI2 meridian.
The DMSP IVM data has the same names as the normal Instrument, and the SAMI2
data has the same name as the SAMI2 data, but with ``model_`` as a prefix to
prevent confusion.  You can change this prefix using the ``model_label`` keyword
argument, allowing multiple models to be matched to the same observational data
set.
   
   
::

   # Using the results from the prior example
   print([cc for cc in matched_inst.data.variables.keys() if cc.find("model_") == 0])


This produces the output line: ``['model_deni', 'model_vsi', 'model_ti', 'model_te', 'model_slt']``


You can also match model and data results by location alone.  This is done by
setting the ``time_method`` keyword arguement to ``'max'``.
