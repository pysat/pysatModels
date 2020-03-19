Loading Model Data
==================

pysatModels uses `pysat <https://github.com/pysat/pysat>`_ to load modelled data
sets.  As specified in the
`pysat tutorial <https://pysat.readthedocs.io/en/latest/tutorial.html>`_,
data may be loaded using the following commands.  TIE-GCM is used as an
example, and so to execute this code snippet the user will need to obtain a
TIE-GCM data file from <UCAR <https://www.hao.ucar.edu/modeling/tgcm/tie.php>`_.

::
   import pysat
   import pysatModels as ps_mod

   filename = 'tiegcm_filename.nc'
   tiegcm = pysat.Instrument(platform='ucar', name='tiegcm')
   tiegcm.load(fname=filename)

   
