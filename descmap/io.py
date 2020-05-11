"""Input/Output functionality"""

import os
import json
from pathlib import Path

import numpy as np
from pmutt.io.excel import read_excel

from descmap.job import get_descriptors_shape, get_tot_job_num

def write_out_file(in_path, out_path, fields):
    """Helper method to write submission files
    
    Parameters
    ----------
        in_path : str
            Path to read the input file to modify. Fields to be replace should
            be proceeded by '__'. e.g. __omkm_path.
        out_path : str
            Path to write the output file with modified fields filled in.
        fields : dict
            Keys of dictionary are the fields to replace (without '__') and the
            values are the strings to replace the text.        
    """

    '''Read template lines'''
    with open(in_path, 'r') as in_ptr:
        out_lines = ''.join(in_ptr.readlines())

    '''Process output lines'''
    for field, val in fields.items():
        out_lines = out_lines.replace('__{}'.format(field), str(val))
    
    '''Write output lines'''
    with open(out_path, 'w', newline='\n') as out_ptr:
        out_ptr.write(out_lines)

def organize_excel_inputs(in_path='./inputs.xlsx'):
    """Reads Excel sheet and creates accessible JSON
    
    Parameters
    ----------
        in_path : str
            Path to Excel sheet
    Returns
    -------
        out_dict : dict
            Output dictionary whose keys are important sheets. The values types
            vary based on information found in the spreadsheet.
    """
    in_path = Path(in_path).as_posix()
    '''Organize inputs into a dictionary'''
    out_dict = {
        'reactor': read_excel(in_path, sheet_name='reactor')[0],
        'units': read_excel(in_path, sheet_name='units')[0],
        'descriptors': read_excel(in_path, sheet_name='descriptors'),
        'fields': read_excel(in_path, sheet_name='fields')[0],
        # 'job': read_excel(in_path, sheet_name='job')[0],
        'analysis': read_excel(in_path, sheet_name='analysis')[0],
        'phases': read_excel(in_path, sheet_name='phases'),
        'options': read_excel(in_path, sheet_name='options')[0],
        'id': {}
    }

    '''Create folders'''
    Path('./setup').mkdir(exist_ok=True, parents=True)
    Path('./omkm').mkdir(exist_ok=True, parents=True)
    Path('./analysis').mkdir(exist_ok=True, parents=True)
    Path('./log').mkdir(exist_ok=True, parents=True)

    '''Add extra parameters for fields'''
    n_col = [row['n'] for row in out_dict['descriptors']]
    sampling = out_dict['fields']['sampling']
    out_dict['fields']['n_jobs'] = get_tot_job_num(n_col=n_col, sampling=sampling)
    out_dict['fields']['desc_shape'] = get_descriptors_shape(n_col=n_col,
                                                             sampling=sampling)
    out_dict['fields']['n_concurrent'] = int(np.min([out_dict['fields']['n_jobs'],
                                                  out_dict['fields']['n_concurrent']]))
    return out_dict

# def get_rel_paths(paths, start, suffix=''):
#     """Helper method to return relative paths

#     Parameters
#     ----------
#         paths : dict
#             Paths to convert to relative paths. The key is a str identifying
#             what the path represents and the value is the path
#         start : str
#             Path to start from
#         suffix : str, optional
#             Suffix to add onto end. Default is ''
#     Returns
#     -------
#         rel_paths : dict
#             Relative paths.
#     """
#     rel_paths = {}
#     for name, path in paths.items():
#         rel_path = os.path.relpath(path=path, start=start)
#         # Standardize Path
#         rel_path = Path(rel_path).as_posix()
#         rel_paths['{}{}'.format(name, suffix)] = rel_path
#     return rel_paths
