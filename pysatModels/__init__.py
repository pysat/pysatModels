#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2019, AGB & pysat team
# Full license can be found in License.md
# -----------------------------------------------------------------------------
"""Core library for pysatModels.

Model utilities designed to facilitate studies that integrate observational
and modelled data sets.

"""

import logging
import os

# Import key modules and skip F401 testing in flake8
from pysatModels import models  # noqa: F401
from pysatModels import utils  # noqa: F401

# Set the version
local_dir = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(local_dir, 'version.txt')) as version_file:
    __version__ = version_file.read().strip()

del local_dir, version_file

# Define a logger object to allow easier log handling
logging.raiseExceptions = False
logger = logging.getLogger('pysatModels_logger')
handler = logging.StreamHandler()
formatter = logging.Formatter('%(name)s %(levelname)s: %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.WARNING)
