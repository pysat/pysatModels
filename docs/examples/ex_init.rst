Loading Model Data
==================

Load Model Data into a pysat Instrument
---------------------------------------
pysatModels uses `pysat <https://github.com/pysat/pysat>`_ to load modelled data
sets.  As specified in the
`pysat tutorial <https://pysat.readthedocs.io/en/latest/tutorial.html>`_,
data may be loaded using the following commands.  TIE-GCM is used as an
example, and so to execute this code snippet the user will need to obtain a
TIE-GCM data file from `UCAR <https://www.hao.ucar.edu/modeling/tgcm/tie.php>`_.

::

   
   import pysat
   import pysatModels as ps_mod

   filename = 'tiegcm_filename.nc'
   tiegcm = pysat.Instrument(platform='ucar', name='tiegcm')
   tiegcm.load(fname=filename)


Load Model Data into an xarray Dataset
--------------------------------------

There are situations (such as when developing a new model) when it may be
inconvenient to create a pysat Instrument object for a modelled data set.  Many
of the pysatModels utilities allow xarray Datasets as input.  For these routines
or to retrieve an xarray Dataset for other purposes, you can use the
``load_model_xarray`` routine in
`utils.match <../utils.html#module-pysatModels.utils.match>`_.

In this example, the time is irrelevent because a full filename is provided:

::

   import datetime as dt
   import pysat
   import pysatModels as ps_mod

   ftime = dt.datetime(2010, 1, 1)
   filename = 'tiegcm_filename.nc'
   tiegcm = pysat.Instrument(platform='ucar', name='tiegcm')

   tg_dataset = ps_mod.utils.match.load_model_xarray(ftime, tiegcm, filename)


In this example, the filename includes temporal information, which is provided
within the loading function by the input time:

::

   import datetime as dt
   import pysat
   import pysatModels as ps_mod

   ftime = dt.datetime(2010, 1, 1)
   filename = 'tiegcm_%Y%j.nc'
   tiegcm = pysat.Instrument(platform='ucar', name='tiegcm')

   tg_dataset = ps_mod.utils.match.load_model_xarray(ftime, tiegcm, filename)


In this example, the routine takes advantage of the pysat file organization
system, and will return a ``NoneType`` object if no files are found for the
specified time:

   
::

   import datetime as dt
   import pysat
   import pysatModels as ps_mod

   ftime = dt.datetime(2010, 1, 1)
   tiegcm = pysat.Instrument(platform='ucar', name='tiegcm')

   tg_dataset = ps_mod.utils.match.load_model_xarray(ftime, tiegcm)
