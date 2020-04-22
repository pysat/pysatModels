#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2019, AGB & pysat team
# Full license can be found in License.md
# -----------------------------------------------------------------------------
"""
pysatModels.models
======================

Routines for loading model data into a pysat Instrument object

"""

from __future__ import absolute_import
from __future__ import unicode_literals

# Import key modules and skip F401 testing in flake8
from pysatModels.models import ucar_tiegcm  # noqa: F401

__all__ = ['ucar_tiegcm']
