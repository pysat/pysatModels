#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2022, pysat development team
# Full license can be found in License.md
# -----------------------------------------------------------------------------
"""Routines for loading model data into a pysat Instrument object."""

# Import key modules and skip F401 testing in flake8
from pysatModels.models import methods  # noqa: F401
from pysatModels.models import pydineof_dineof  # noqa: F401
from pysatModels.models import sami2py_sami2  # noqa: F401
from pysatModels.models import ucar_tiegcm  # noqa: F401

__all__ = ['pydineof_dineof', 'sami2py_sami2', 'ucar_tiegcm']
