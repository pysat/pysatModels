.. _ex_extract:

Extract Observational-Style Data
================================

Comparison of model and Instrument data is supported in pysatModels, in part,
by enabling the extraction (or interpolation) of model output onto the same
locations as an observation-style data set. One common example is
'flying' a satellite through a model. The satellite locations are used
to extract relevant model data enabling direct comparison of observed and
modeled values.


Regular Grid Models
-------------------

:py:func:`pysatModels.utils.extract.extract_modeled_observations` supports
extracting values from models on a regular grid. The function can linearly
inteprolate model values onto instrument locations or use the nearest modeled
location. Uses :py:func:`scipy.interpolate.interpn` as the underlying
interpolation function.

Leaving description of @aburrell's function to @aburrell.


:py:func:`pysatModels.utils.extract.instrument_view_through_model` supports
interpolating values from regular grid models onto Instrument locations using
:py:func:`scipy.interpolate.RegularGridInterpolator`. Consider the following
example that interpolates model data onto a satellite data set using
pysat testing data sets.

.. code:: python

   import pysat
   import pysatModels

   # Load simulated satellite Instrument data set
   inst = pysat.Instrument('pysat', 'testing', max_latitude=45.)
   inst.load(2009, 1)

   # Load simulated regular-grid model Instrument
   model = pysat.Instrument('pysat', 'testmodel')
   model.load(2009, 1)

Looking at the loaded `model.data` we can see that the model is indeed regular.

.. code:: python

   In []: model.data
   Out[]:
   <xarray.Dataset>
   Dimensions:    (time: 96, latitude: 21, longitude: 73, altitude: 41)
   Coordinates:
     * time       (time) datetime64[ns] 2009-01-01 ... 2009-01-01T23:45:00
     * latitude   (latitude) float64 -50.0 -45.0 -40.0 -35.0 ... 40.0 45.0 50.0
     * longitude  (longitude) float64 0.0 5.0 10.0 15.0 ... 345.0 350.0 355.0 360.0
     * altitude   (altitude) float64 300.0 305.0 310.0 315.0 ... 490.0 495.0 500.0
   Data variables:
       uts        (time) float64 0.0 900.0 1.8e+03 ... 8.37e+04 8.46e+04 8.55e+04
       slt        (time, longitude) float64 0.0 0.3333 0.6667 ... 23.08 23.42 23.75
       mlt        (time, longitude) float64 0.2 0.5333 0.8667 ... 23.28 23.62 23.95
       dummy1     (time, latitude, longitude) float64 0.0 0.0 0.0 ... 0.0 3.0 6.0
       dummy2     (time, latitude, longitude, altitude) float64 0.0 0.0 ... 18.0

The `Coordinates` `time`, `latitude`, `longitude`, and `altitude` are all
one-dimensional and directly relevant to a physical satellite location. The
equivalent satellite variables are `latitude`, `longitude`, and `altitude`,
with `time` taken from the Instrument's associated datetimeindex
(`inst.data.index`).

.. code:: python

   In []: inst.variables
   Out[]:
   Index(['uts', 'mlt', 'slt', 'longitude', 'latitude', 'altitude', 'orbit_num',
          'dummy1', 'dummy2', 'dummy3', 'dummy4', 'string_dummy', 'unicode_dummy',
          'int8_dummy', 'int16_dummy', 'int32_dummy', 'int64_dummy',
          'model_dummy2'],
         dtype='object')

   In []: inst.data.index
   Out[]:
   DatetimeIndex(['2009-01-01 00:00:00', '2009-01-01 00:00:01',
                  '2009-01-01 00:00:02', '2009-01-01 00:00:03',
                  '2009-01-01 00:00:04', '2009-01-01 00:00:05',
                  '2009-01-01 00:00:06', '2009-01-01 00:00:07',
                  '2009-01-01 00:00:08', '2009-01-01 00:00:09',
                  ...
                  '2009-01-01 23:59:50', '2009-01-01 23:59:51',
                  '2009-01-01 23:59:52', '2009-01-01 23:59:53',
                  '2009-01-01 23:59:54', '2009-01-01 23:59:55',
                  '2009-01-01 23:59:56', '2009-01-01 23:59:57',
                  '2009-01-01 23:59:58', '2009-01-01 23:59:59'],
                 dtype='datetime64[ns]', name='Epoch', length=86400, freq=None)


Interpolating `model` data onto `inst` is accomplished via

.. code:: python

   new_data_keys = pysatModels.utils.extract.instrument_view_through_model(inst,
                              model.data, ['latitude', 'longitude', 'altitude'],
                              ['latitude', 'longitude', 'altitude'], 'time',
                              'time', ['deg', 'deg', 'km'], ['dummy2'])

where `inst` and `model.data` provide the required :py:class:`pysat.Instrument`
object and :py:class:`xarray.Dataset`. The ::

   ['latitude', 'longitude', 'altitude']

term provides the content and ordering of the coordinates for model variables
to be interpolated. The subsequent ::

   ['latitude', 'longitude', 'altitude']

term provides the equivalent content from the satellite's data set, in the same
order as the model coordinates. In this case, the same lables are used for
both the satellite and modeled data sets. The ::

   'time', 'time'

terms cover the model labels used for time variable and coordinate. The ::

   ['deg', 'deg', 'km']

term covers the units for the model dimensions (latitude/longitude/altitude).
Units for the corresponding information from `inst` are taken directly from the
:py:class:`pysat.Instrument` object. The final presented input::

    ['dummy2']

is a list of model variables that will be interpolated onto `inst`.

The results of ::

    inst[new_data_keys].plot(title='Interpolation Example')

are shown below.



Irregular Grid Models
---------------------

Some models aren't on a regular grid, or may not be a regular grid across
the coordinates of interest. Consider an alternative model data set,

.. code:: python

    In []: model = pysat.Instrument('pysat', 'testmodel', tag='pressure_levels')

    In []: model.load(2009, 1)

    In []: model.data
    Out[]:
    <xarray.Dataset>
    Dimensions:       (time: 24, latitude: 72, longitude: 144, lev: 57, ilev: 57)
    Coordinates:
      * time          (time) datetime64[ns] 2009-01-01 ... 2009-01-01T23:00:00
      * latitude      (latitude) float64 -88.75 -86.25 -83.75 ... 83.75 86.25 88.75
      * longitude     (longitude) float64 -180.0 -177.5 -175.0 ... 172.5 175.0 177.5
      * lev           (lev) float64 -7.0 -6.75 -6.5 -6.25 -6.0 ... 6.25 6.5 6.75 7.0
      * ilev          (ilev) float64 -6.875 -6.625 -6.375 ... 6.625 6.875 7.125
    Data variables:
        uts           (time) float64 0.0 3.6e+03 7.2e+03 ... 7.92e+04 8.28e+04
        altitude      (time, ilev, latitude, longitude) float64 0.0 0.0 ... 5.84e+07
        dummy_drifts  (time, ilev, latitude, longitude) float64 0.0 0.0 ... 83.01
        slt           (time, longitude) float64 12.0 12.17 12.33 ... 10.67 10.83
        mlt           (time, longitude) float64 12.2 12.37 12.53 ... 10.87 11.03
        dummy1        (time, latitude, longitude) float64 0.0 0.0 0.0 ... 0.0 9.0

Model variables, such as `dummy_drifts`, are regular over
`(time, ilev, latitude, longitude)`, where `ilev` is a constant pressure level.
Unfortunately, the observational data in `inst` doesn't contain pressure level
as a simulated/measured parameter. However, `altitude` is present in the model
data but varies over all four coordinates. Interpolating `dummy_drifts`
onto `inst` requires either adding an appropriate value for `ilev` into `inst`,
or iterpolating model variables using the irregular variable `altitude` instead
of `ilev`.

Altitude to Pressue
^^^^^^^^^^^^^^^^^^^

:py:func:`pysatModels.utils.extract.instrument_altitude_to_model_pressure`
will use information in a model to generate approrpiate pressure levels for a
supplied altitude in an observational-like data set.

.. code:: python

    from pysatModels.utils.extract import instrument_altitude_to_model_pressure as iamp
    keys  = iamp(inst, model.data, ["altitude", "latitude", "longitude"],
                 ["ilev", "latitude", "longitude"],
                 "time", "time", ['', "deg", "deg"],
                 'altitude', 'altitude', 'cm')

The function will guess a pressure level for all locations in `inst` and then
use the regular mapping from pressure to altitude to obtain the equivalent altitude
from the model. The pressure is adjusted up/down an increment based upon the
comparison and the process is repeated until the target tolerance (default
is 1 km) is achieved. The keys for the model derived pressure and altitude
values added to `inst` are returned from the function.

.. code:: python

    In []: keys
    Out[]: ['model_altitude', 'model_pressure']

    In []: inst['model_pressure']
    Out[]:
    Epoch
    2009-01-01 00:00:00    3.104662
    2009-01-01 00:00:01    3.104652
    2009-01-01 00:00:02    3.104642
    2009-01-01 00:00:03    3.104632
    2009-01-01 00:00:04    3.104623
                             ...
    2009-01-01 23:59:55    2.494845
    2009-01-01 23:59:56    2.494828
    2009-01-01 23:59:57    2.494811
    2009-01-01 23:59:58    2.494794
    2009-01-01 23:59:59    2.494776
    Name: model_pressure, Length: 86400, dtype: float64

    In []: inst['model_altitude'] - inst['altitude']
    Out[]:
    Epoch
    2009-01-01 00:00:00   -0.744426
    2009-01-01 00:00:01   -0.744426
    2009-01-01 00:00:02   -0.744425
    2009-01-01 00:00:03   -0.744424
    2009-01-01 00:00:04   -0.744424
                             ...
    2009-01-01 23:59:55   -0.610759
    2009-01-01 23:59:56   -0.610757
    2009-01-01 23:59:57   -0.610754
    2009-01-01 23:59:58   -0.610751
    2009-01-01 23:59:59   -0.610749
    Length: 86400, dtype: float64

Using the added `model_pressure` information model values may not be interpolated
onto `inst` using regular grid methods.

.. code:: python

    new_keys = pysatModels.utils.extract.instrument_view_through_model(inst,
               model.data, ['model_pressure', 'latitude', 'longitude'],
               ['ilev', 'latitude', 'longitude'], 'time', 'time',
               ['', 'deg', 'deg'], ['dummy_drifts'])

.. code:: python

    In []: new_data_keys
    Out[]: ['model_dummy_drifts']

    In []: inst['model_dummy_drifts']
    Out[]:
    Epoch
    2009-01-01 00:00:00    30.289891
    2009-01-01 00:00:01    30.305303
    2009-01-01 00:00:02    30.320704
    2009-01-01 00:00:03    30.336092
    2009-01-01 00:00:04    30.351469
                             ...
    2009-01-01 23:59:55    63.832658
    2009-01-01 23:59:56    63.868358
    2009-01-01 23:59:57    63.904047
    2009-01-01 23:59:58    63.939724
    2009-01-01 23:59:59    63.975389
    Name: model_dummy_drifts, Length: 86400, dtype: float64

Attached image.

The time to translate altitude to model pressure is ~3 s, and the regular
interpolation takes an additional ~300 ms.

Irregular Variable
^^^^^^^^^^^^^^^^^^

More generally,
:py:func:`pysatModels.utils.extract.interp_inst_w_irregular_model_coord` can
deal with irregular coordinates when interpolating onto an observational-like
data set using :py:func:`scipy.interpolate.griddata`. The `model` loaded above
is regular against pressure level, latitude, and longitude, however it is
irregular with

.. code:: python

    keys = pysatModels.utils.extract.interp_inst_w_irregular_model_coord(inst,
                model.data, ["altitude", "latitude", "longitude"],
                ["ilev", "latitude", "longitude"],
                "time", ["cm", "deg", "deg"], "ilev",
                "altitude", [50., 10., 10.],
                sel_name=["dummy_drifts", "altitude"])

where `inst` and `model.data` provide the required :py:class:`pysat.Instrument`
object and :py:class:`xarray.Dataset`. The ::

   ["altitude", "latitude", "longitude"]

term provides the content and ordering of the spatial locations for `inst`.
The subsequent ::

   ["ilev", "latitude", "longitude"]

term provides the equivalent regular dimension labels from `model.data`,
in the same order as the underlying model dimensions. While this function
does operate on irregular data it also needs information on the underlying
regular memory structure of the variables. The ::

   "time"

terms cover the model label used for the datetime coordinate. The ::

   ["cm", "deg", "deg"]

term covers the units for the model information (altitude/latitude/longitude)
that maps to the `inst` information in `["altitude", "latitude", "longitude"]`.
Note that the "cm" covers units for 'altitude' in `model.data`, the variable
that will replace 'ilev', while "deg" and "deg" covers the units for
the "latitude" and "longitude" dimensions. Units for the corresponding information from `inst` are taken directly from the
:py:class:`pysat.Instrument` object. The ::

    "ilev"

identifies the regular model dimension that will be replaced with irregular
data for interpolation. The ::

    "altitude"

identifies the irregular model variable that will replace the regular coordinate.
The ::

    [50., 10., 10.]

term is used to define a half-window for each of the `inst` locations, in units
from `inst`, used to downselect data from `model.data` to reduce computational
requirements. In this case a window of +/-50 km in altitude,
+/-10 degrees in latitude,
and +/-10 degrees in longitude is used. The keyword argument ::

    sel_name=["dummy_drifts", "altitude"]

identifies the `model.data` variables that will be interpolated onto `inst`.

The results of ::

    inst[keys].plot(title='Interpolation Example')

are shown below.





Text below saved as a development reference while writing.

Routines to extract observational-style data from model output.

Model verification and validation is supported in pysatModels through
`pyForecastTools <https://github.com/drsteve/PyForecastTools>`_. Many of the
standard statistics used for these purposes can be run in a single go using
the utility :py:func:`pysatModels.utils.compare.compare_model_and_inst`. The
following example uses the paired data produced in example :ref:`ex_match-loc`.
The :py:data:`matched_inst` object at this point should display as shown below.

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

Note the statistical output is in the units of the observed data set.  The
:py:data:`stat_dict` output is a *dict* with the observed data
variable name(s) as the first set of keys and the requested statistics for
each data type as a nested *dict*.


Not all of the statistics were appropriate for the data set, as indicated by the
:py:exc:`RuntimeWarning` messages seen when running
:py:func:`~pysatModels.utils.compare.compare_model_and_inst`.  The values
show that, unsurprisingly, the random data from the test model file does not
agree well with the C/NOFS meridional **E** x **B** drifts.
