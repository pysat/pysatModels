.. _install:

Installation
============

The following instructions will allow you to install pysatModels.


.. _install-prereq:

Prerequisites
-------------

.. image:: images/poweredbypysat.png
    :width: 150px
    :align: right
    :alt: powered by pysat Logo, blue planet with orbiting python


pysatModels uses common Python modules, as well as modules developed by and for
the Space Physics community.  This module officially supports Python 3.6+.

 ============== =================
 Common modules Community modules
 ============== =================
  numpy         pysat>=3.0.2
  pandas        pyForecastTools
  requests
  scipy
  xarray
 ============== =================


.. _install-opt:

Installation Options
--------------------

You may now obtain pysatModels from PyPi, Zenodo, or the GitHub repository.  If
you use PyPi, simply call:

::

   pip install pysatModels


If you use GitHub or Zenodo, you need to first obtain the package and then
install it.

1. Clone the git repository or download the repository from Zenodo and unzip
   the compressed files.  To clone the git repository, use the command below.
::


   git clone https://github.com/pysat/pysatModels.git


2. Install pysatModels from the repository folder, once it is in the desired
   location. Change directories into the repository folder and run the setup.py
   file. There are a few ways you can do this:

   A. Install on the system (root privileges required)::


        sudo python3 setup.py install
   B. Install at the user level::


        python3 setup.py install --user
   C. Install with the intent to develop locally::


        python3 setup.py develop --user

.. _post-install:
Post Installation
-----------------

After installation, you may register the :py:mod:`pysatModel` model
:py:class:`Instrument` sub-modules with pysat.  If this is your first time using
pysat, check out the `quickstart guide
<https://pysat.readthedocs.io/en/latest/quickstart.html>`_ for pysat. Once pysat
is set up, you may choose to register the the :py:mod:`pysatModel` model
:py:class:`Instruments` sub-modules by:

.. code:: python


   import pysat
   import pysatModels as pymod

   pysat.utils.registry.register_by_module(pymod.models)

You may then use the pysat :py:attr:`platform` and :py:attr:`name` keywords to
initialize the model :py:class:`Instrument` instead of the
:py:attr:`inst_module` keyword argument.
