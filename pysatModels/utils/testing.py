#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2022, pysat development team
# Full license can be found in License.md
# -----------------------------------------------------------------------------
"""Supporting functions for testing user inputs."""


def compare_mod_name_coordinates(data, mod_name):
    """Compare ordering of mod_name against data coordinates, ignoring time.

    First dimension for `data` needs to be the main DateTimeIndex.

    Parameters
    ----------
    data : xarray.Dataseries
        Variable from xarray to compare against `mod_name`.
    mod_name : list
        List of coordinate names.

    Raises
    ------
    ValueError if `mod_name` and `data` aren't consistent.

    """

    # Get a list of all spatial coordinates, time is presumed to be first.
    # Coordinates aren't in actual data ordering though :(
    # Get the ordering of coordinates by using the ordering of data.dims.
    # Collect dims for coordinates and then find the association with data.dims
    # to construct correct order for data coordinates.
    raw_coords = list(data.coords)
    coord_dims = []
    for coord in raw_coords:
        coord_dims.append(data[coord].dims[0])

    # Get list of dimension names
    dims = list(data.dims)

    # Check first dimension is full of datetime data
    if str(data[data.dims[0]].data.dtype).find('datetime') < 0:
        estr = ''.join(['First dimension of ', data.dims[0], ' does not ',
                        'appear to be a time dimension.'])
        raise ValueError(estr)
    else:
        # Drop time-like dimension
        dims = dims[1:]

    # Construct dim order for coordinates
    coords = []
    for dim in dims:
        for i, coord in enumerate(coord_dims):
            if dim == coord:
                coords.append(raw_coords[i])

    # Check actual coordinates all present in `mod_name`
    missing_coords = []
    for coord in coords:
        if coord not in mod_name:
            missing_coords.append(coord)

    if len(missing_coords) > 0:
        estr = ''.join(['Provided coordinates ', repr(mod_name), ' are not ',
                        'all within variable coordinates ',
                        repr(missing_coords), '.'])
        raise ValueError(estr)

    # Check number of provided coords against known data coordinates.
    if len(mod_name) > len(coords):
        estr = 'Provided too many coordinates in `mod_name`.'
        raise ValueError(estr)

    # Check order for both sets of coordinates.
    for i, name in enumerate(mod_name):
        if name != coords[i]:
            estr = ''.join(['Provided coordinates ', repr(mod_name),
                            ' not in same order as variable coordinates ',
                            repr(coords)])
            raise ValueError(estr)

    return
