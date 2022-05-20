.. _ex-init:

Loading Model Data
==================


.. _ex-init-loadinst:

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
   tiegcm = pysat.Instrument(inst_module=ps_mod.models.ucar_tiegcm)
   tiegcm.load(fname=filename)


.. _ex-init-loadxr:

Load Model Data into an xarray Dataset
--------------------------------------

There are situations (such as when developing a new model) when it may be
inconvenient to create a pysat Instrument object for a modelled data set.  Many
of the pysatModels utilities allow :py:class:`xarray.Dataset` objects as input.
For these routines or to retrieve an :py:class:`xarray.Dataset` for other
purposes, you can use the :py:func:`load_model_xarray` routine in
:ref:`utils-match`.

In this example, the time is irrelevent because a full filename is provided:

::

   import datetime as dt
   import pysat
   import pysatModels as ps_mod

   # Data directory definition is needed if you don't save the TIE-GCM file
   # to your pysat_data/ucar/tiegcm/ directory.  This definition assumes the
   # file specified by `filename` lives in your current working directory.
   data_dir = '.'

   # Define the file name, time, and initialize the instrument
   ftime = dt.datetime(2010, 1, 1)
   filename = 'tiegcm_filename.nc'
   tiegcm = pysat.Instrument(platform='ucar', name='tiegcm', data_dir=data_dir)

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
system, and will return a :py:class:`NoneType` object if no files are found for
the specified time:

   
::

   import datetime as dt
   import pysat
   import pysatModels as ps_mod

   ftime = dt.datetime(2010, 1, 1)
   tiegcm = pysat.Instrument(platform='ucar', name='tiegcm')

   tg_dataset = ps_mod.utils.match.load_model_xarray(ftime, tiegcm)
