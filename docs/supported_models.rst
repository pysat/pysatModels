.. _supmodels:

Supported Models
================


.. _supmodels-dineof:

DINEOF
------

Supports the pyDINEOF model output.  pyDINEOF is a Python implementation of the
Data INterpolating Emperical Orthogonal Function method. Information about this
package can obtained by contacting
`Russell Stoneback <https://github.com/rstoneback>`_.

.. automodule:: pysatModels.models.pydineof_dineof
   :members:


.. _supmodels-sami2:

SAMI2
-----

Supports the SAMI2 (Sami2 is Another Model of the Ionosphere 2) model through the
sami2py interface. Sami2py is a python module that runs the SAMI2 model, as well
as archives, loads and plots the resulting modeled values. SAMI2 is a model
developed by the Naval Research Laboratory to simulate the motions of plasma
in a 2D ionospheric environment along a dipole magnetic field
[Huba et al, 2000].  Information about this model can be found at the
`sami2py github page <https://github.com/sami2py/sami2py>`_,
along with a list of the
`SAMI2 principle papers <https://sami2py.readthedocs.io/en/latest/introduction.html#references>`_.

.. automodule:: pysatModels.models.sami2py_sami2
   :members:


.. _supmodels-tiegcm:

TIE-GCM
-------

Supports the UCAR (University Corporation for Atmospheric Research) model,
Thermosphere-Ionosphere-Electrodynamics General Circulation Model (TIE-GCM).
Information about this model can be found at the
`UCAR TIE-GCM website <https://www.hao.ucar.edu/modeling/tgcm/tie.php>`_,
along with a list of the principle papers and references.

.. automodule:: pysatModels.models.ucar_tiegcm
   :members:
