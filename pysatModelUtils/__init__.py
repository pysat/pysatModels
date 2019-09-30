#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2019, AGB & pysat team
# Full license can be found in License.md
#-----------------------------------------------------------------------------
"""
pysatModelUtils
===============

Model utilities designed to facilitate studies that integrate observational
and modelled data sets.

"""

from __future__ import absolute_import
from __future__ import unicode_literals

import logging
import os

from . import (utils)
from . import (models)

# set the version
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'version.txt')) as version_file:
    __version__ = version_file.read().strip()

# Define a logger object to allow easier log handling
logging.raiseExceptions = False
logger = logging.getLogger('pysatModelUtils_logger')
