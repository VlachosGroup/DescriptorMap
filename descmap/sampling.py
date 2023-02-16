"""Functionality related to sampling"""
from warnings import warn
from collections import namedtuple

import numpy as np
from scipy.stats import variation
import pandas as pd
from pyDOE import lhs

DescPoint = namedtuple('DescPoint',
                       ('name', 'i', 'val', 'field_len', 'sampling'))

def get_linspace_map(descriptors_data):
    """Returns linear spacing mapping
    
    Parameters
    ----------
        descriptors_data : list of dict
            Each element of the list has information related to a descriptor.
            The dictionary is expected to have the fields: ``name``,
            ``low_value``, ``high_value``, ``n``.
    Returns
    -------
        linspace_map : list of list of DescPoint namedtuples
            Linear space mapping of all the descriptors. Each element of the
            top list has info about each descriptor. The next list contains the
            runs. The DescPoint namedtuples has the attribute ``name``, ``i``,
            ``val``, and ``field_len``.

            e.g. [
                    [('A', 0, 0.1, 1), ('A', 1, 0.2, 1), ('A', 2, 0.3, 1)],
                    [('B', 0, 0.5, 1), ('B', 1, 1., 1)]
                 ]
    """
    linspace_map = []
    for descriptor_data in descriptors_data:
        # Individual descriptor mapping
        desc_map = [] 
        name = descriptor_data['name']
        field_len = len(str(descriptor_data['n']))
        for i, val in enumerate(np.linspace(descriptor_data['low_value'],
                                            descriptor_data['high_value'],
                                            descriptor_data['n'])):
            desc_point = DescPoint(name=name, i=i, val=val, field_len=field_len,
                                   sampling='linear')
            desc_map.append(desc_point)
        linspace_map.append(desc_map)
    return linspace_map

def get_lhs_map(descriptors_data):
    """Returns Latin hypercube sampling mapping
    
    Parameters
    ----------
        descriptors_data : list of dict
            Each element of the list has information related to a descriptor.
            The dictionary is expected to have the fields: ``name``,
            ``low_value``, ``high_value``, ``n``. Note that all ``n`` are
            expected to be the same. If there is a discrepancy, assumes the
            max value.
    Returns
    -------
        lhs_map : list of list of DescPoint namedtuples
            Latin hypercube sampling mapping of all the descriptors. Each
            element of the top list has info about each descriptor.
            The next list contains the runs. The DescPoint namedtuple has the
            attributes ``name``, ``i``, ``val``, and ``field_len``.

            e.g. [
                    [('A', 0, 0.37, 1), ('A', 1, 0.98, 1), ('A', 2, 0.24, 1)],
                    [('B', 0, 0.60, 1), ('B', 1, 0.03, 1), ('B', 1, 0.95, 1)]
                 ]
    """
    '''Determine number of samples'''
    descriptors_n = np.array([row['n'] for row in descriptors_data])
    samples = np.max(descriptors_n)
    # Warn user if # of data points are inconsistent
    if not np.isclose(variation(descriptors_n), 0.):
        warn_msg = ('Number of samples for each descriptor in Latin hypercube '
                    'sampling must remain constant. Using maximum value, {}.'
                    ''.format(samples))
        warn(warn_msg)

    '''Determine number of descriptors'''
    n = len(descriptors_data)

    '''Get field length'''
    field_len = len(str(samples))

    '''Create the LHS array'''
    lhs_array = lhs(n=n, samples=samples)

    '''Format the LHS array to desired mapping format'''
    lhs_map = []
    for i, descriptor_data in enumerate(descriptors_data):
        # Individual descriptor mapping
        desc_map = []
        name = descriptor_data['name']
        low_val = descriptor_data['low_value']
        high_val = descriptor_data['high_value']
        val_range = high_val - low_val
        for j, val in enumerate(lhs_array[:, i]):
            out_val = val*val_range + low_val
            desc_point = DescPoint(name=name, i=j, val=out_val,
                                   field_len=field_len, sampling='lhs')
            desc_map.append(desc_point)
        lhs_map.append(desc_map)
    return lhs_map

sampling_map = {
    'linear': get_linspace_map,
    'lhs': get_lhs_map
}
"""dict: Keys represent sampling type. Values represent function handles."""

