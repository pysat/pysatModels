<div align="left">
        <img height="0" width="0px">
        <img width="20%" src="docs/images/pysatModels.png" alt="pysatModels" title="pysatModels"</img>
</div>

# pysatModels
[![Documentation Status](https://readthedocs.org/projects/pysatmodels/badge/?version=latest)](https://pysatmodels.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://github.com/github/docs/actions/workflows/main.yml/badge.svg)](https://github.com/github/docs/actions/workflows/main.yml/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/pysat/pysatModels/badge.svg?branch=main)](https://coveralls.io/github/pysat/pysatModels?branch=main)

pysatModels handles model-centric data loading through pysat and contains a
variety of tools to perform model-data analysis, including model validation.

Come join us on Slack! An invitation to the pysat workspace is available 
in the 'About' section of the
[pysat GitHub Repository.](https://github.com/pysat/pysat)

# Installation

The following instructions will allow you to install pysatModels and give some
examples on how to use the routines

## Prerequisites

pysatModels uses common Python modules, as well as modules developed by and for
the Space Physics community.  This module officially supports Python 3.6+.

| Common modules | Community modules |
| -------------- | ----------------- |
| numpy          | pysat             |
| pandas         | pyForecastTools   |
| requests       |                   |
| scipy          |                   |
| xarray         |                   |

## Installation through PyPi

```
pip install pysatmodels
```

## Installation from GitHub

1. Clone the git repository:
   `git clone https://github.com/pysat/pysatModels.git`

2. Install pysatModels:
   Change directories into the repository folder and run the setup.py file.
   There are a few ways you can do this:
   A. Install on the system (root privileges required):
      `sudo python3 setup.py install`
   B. Install at the user level:
      `python3 setup.py install --user`
   C. Install with the intent to develop locally:
      `python3 setup.py develop --user`

# Running the code

To allow pysat to manage your model data, you will need to specify a directory
to store data.  This only needs to be done the first time pysat is run.  In
python, run:

```
pysat.params['data_dirs'] = 'path/to/directory/that/exists'
```

Examples for using the different pysatModels utilities are included in the
[documentation](https://pysatmodels.readthedocs.io/en/latest/).
