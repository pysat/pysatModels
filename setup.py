#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2019, Authors
# Full license can be found in License.md and AUTHORS.md
# -----------------------------------------------------------------------------

import os
import codecs
from setuptools import setup, find_packages


here = os.path.abspath(os.path.dirname(__file__))
with codecs.open(os.path.join(here, 'description.txt'), encoding='utf-8') as f:
    long_description = f.read()
version_filename = os.path.join('pysatModels', 'version.txt')
with codecs.open(os.path.join(here, version_filename)) as version_file:
    version = version_file.read().strip()

# change setup.py for readthedocs - commented for now
# on_rtd = os.environ.get('READTHEDOCS') == 'True'
install_requires = ['pysat', 'scipy', 'pandas', 'xarray', 'numpy']

# Run setup

setup(name='pysatModels',
      version=version,
      url='github.com/pysat/pysatModels',
      author='Angeline G. Burrell, Jeff Klenzing, Russell Stoneback',
      author_email='angeline.burrell@nrl.navy.mil',
      description='Model-data comparisons for the pysat ecosystem',
      long_description=long_description,
      packages=find_packages(),
      classifiers=[
          "Development Status :: 4 - Beta",
          "Topic :: Scientific/Engineering :: Physics",
          "Intended Audience :: Science/Research",
          "License :: BSD",
          "Natural Language :: English",
          "Programming Language :: Python :: 2.7",
          "Programming Language :: Python :: 3.5",
          "Programming Language :: Python :: 3.6",
          "Programming Language :: Python :: 3.7",
          "Programming Language :: Python :: 3.8",
          "Operating System :: MacOS :: MacOS X",
      ],
      include_package_data=True,
      zip_safe=False,
      install_requires=install_requires,
      )
