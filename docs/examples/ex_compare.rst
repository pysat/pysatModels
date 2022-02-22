.. _ex_compare:

Compare Paired Data Sets
========================

Model verification and validation is supported in pysatModels through
`pyForecastTools <https://github.com/drsteve/PyForecastTools>`_. Many of the
standard statistics used for these purposes can be run in a single go using
the utility :py:func:`pysatModels.utils.compare.compare_model_and_inst`. The
following example uses the paired data produced in example :ref:`ex_match-loc`.
The :py:data:`matched_inst` object at this point should display as shown below.

.. code:: python

	  print(matched_inst)
	  pysat Instrument object
	  -----------------------
	  Platform: 'cnofs'
	  Name: 'ivm'
	  Tag: ''
	  Instrument id: ''

	  Data Processing
	  ---------------
	  Cleaning Level: 'clean'
	  Data Padding: None
	  Keyword Arguments Passed to list_files: {}
	  Keyword Arguments Passed to load: {}
	  Keyword Arguments Passed to preprocess: {}
	  Keyword Arguments Passed to download: {}
	  Keyword Arguments Passed to list_remote_files: {}
	  Keyword Arguments Passed to clean: {}
	  Keyword Arguments Passed to init: {}
	  Custom Functions: 1 applied
	      0: <function update_longitude at 0x12fcb8310>
	       : Kwargs={'low': -180.0, 'lon_name': 'glon', 'high': 180.0}

          Local File Statistics
	  ---------------------
          Number of files: 930
	  Date Range: 01 January 2009 --- 01 April 2015

	  Loaded Data Statistics
	  ----------------------
	  Date: 01 January 2009
	  DOY: 001
	  Time range: 01 January 2009 00:00:00 --- 01 January 2009 23:59:59
	  Number of Times: 168596
	  Number of variables: 89

	  Variable Names:
	  RPAflag                         driftMeterflag                  ionVelocityX                    
                                              ...                                               
	  ECISC_index1                    LVLHSC_index1                   dineof_model_equator_model_data 

          pysat Meta object
	  -----------------
	  Tracking 15 metadata values
	  Metadata for 92 standard variables
	  Metadata for 0 ND variables


Now, we need to convert this :py:class:`pysat.Instrument` object to an
:py:class:`xarray.Dataset`.  This conversion is needed to simplify the
comparison analysis, since :py:attr:`pysat.Instrument.data` may be either
:py:class:`xarray.Dataset` or :py:class:`pandas.DataFrame` objects.  pysatModels
uses :py:class:`xarray.Dataset` as the base analysis class, because this class
is best suited for modelled output.  This can be easily done using
:py:func:`pysatModels.utils.convert.convert_pysat_to_xarray`. Before we do
that, though, we're going to update the units of the modelled data. This is
necessary for the comparison, since the
:py:func:`pysatModels.utils.compare.compare_model_and_inst` tests to make sure
the paired data have the same units.  It can handle converting between different
units of the same type, so we will specify that the modelled data is a velocity
in *cm/s*, while the observations are a velocity measured in *m/s*.

.. code:: python

	  from pysatModels.utils import convert

          inst_data_keys = ['ionVelmeridional']
          model_data_keys = ['dineof_model_equator_model_data']
	  matched_inst.meta[model_data_keys[0]] = {
              matched_inst.meta.labels.units: "cm/s"}
          paired_data = convert.convert_pysat_to_xarray(matched_inst)
          print(paired_data)

	  <xarray.Dataset>
	  Dimensions:                          (index: 168596)
	  Coordinates:
	    * index                            (index) datetime64[ns] 2009-01-01T00:00:...
	  Data variables: (12/89)
	    RPAflag                          (index) int16 4 4 4 4 4 4 4 ... 4 3 4 4 4 3
	    driftMeterflag                   (index) int16 0 0 0 0 0 0 0 ... 0 0 0 0 0 0
	    ionVelocityX                     (index) float32 nan nan nan ... nan nan nan
	    ionVelocityY                     (index) float32 633.8 589.8 ... 234.5 238.1
	    ionVelocityZ                     (index) float32 106.3 105.6 ... 119.6 118.5
	    vXvariance                       (index) float32 0.0 0.0 0.0 ... 0.0 0.0 0.0
	    ...                               ...
	    meridionalunitvectorX            (index) float32 -0.04526 ... 0.02975
	    meridionalunitvectorY            (index) float32 -0.2083 -0.208 ... -0.392
	    meridionalunitvectorZ            (index) float32 -0.977 -0.9771 ... -0.9195
	    ECISC_index1                     (index) float64 nan nan nan ... nan nan nan
	    LVLHSC_index1                    (index) float64 nan nan nan ... nan nan nan
	    dineof_model_equator_model_data  (index) float64 -0.06359 ... -1.612


Now we can compare the paired model and observed data.  The example below will
only calculate the bias-related statistics, but there are several options for
the :py:data:`method` keyword argument that allow single or groups of statistics
to be calculated in one call.

.. code:: python
	  from pysatModels.utils import compare

	  stat_dict, data_units = compare.compare_model_and_inst(
	      paired_data, inst_data_keys, model_data_keys,
	      methods=['all_bias'], unit_label=matched_inst.meta.labels.units)

	  print(data_units)
	  {'ionVelmeridional': 'm/s'}


Note the statistical output is in the units of the observed data set.  The
:py:data:`stat_dict` output is a *dict* with the observed data
variable name(s) as the first set of keys and the requested statistics for
each data type as a nested *dict*.

.. code:: python
	  print(stat_dict)

	  {'ionVelmeridional': {'symmetricSignedBias': masked,
	                        'meanPercentageError': -100.70728288065469,
				'bias': 18.772352902200083,
				'medianLogAccuracy': masked}}

Not all of the statistics were appropriate for the data set, as indicated by the
:py:exc:`RuntimeWarning` messages seen when running
:py:func:`~pysatModels.utils.compare.compare_model_and_inst`.  The values
show that, unsurprisingly, the random data from the test model file does not
agree well with the C/NOFS meridional **E** x **B** drifts.
