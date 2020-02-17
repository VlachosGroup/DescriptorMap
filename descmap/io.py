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

def organize_excel_inputs(in_path='./input.xlsx'):
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
        'job': read_excel(in_path, sheet_name='job')[0],
        'analysis': read_excel(in_path, sheet_name='analysis')[0],
        'phases': read_excel(in_path, sheet_name='phases'),
        'id': {}
    }

    '''Add path data'''
    paths = read_excel(in_path, sheet_name='paths')[0]
    # For Excel sheet
    paths['excel_in'] = in_path
    paths['json_in'] = Path(in_path).with_suffix('.json').as_posix()

    # For setup
    Path(paths['setup_path']).mkdir(exist_ok=True, parents=True)
    paths['setup_python_in'] = Path(paths['template_path'], 'python',
                                    'setup.py').as_posix()
    paths['setup_python_out'] = Path(paths['setup_path'], 'setup.py').as_posix()
    paths['setup_sub_in'] = Path(paths['template_path'], 'sub', 'setup.qs').as_posix()
    paths['setup_sub_out'] = Path(paths['setup_path'], 'setup.qs').as_posix()

    # For OMKM
    Path(paths['omkm_path']).mkdir(exist_ok=True, parents=True)
    paths['omkm_sub_in'] = Path(paths['template_path'], 'sub', 'JA_omkm.qs').as_posix()
    paths['omkm_sub_out'] = Path(paths['omkm_path'], 'JA_omkm.qs').as_posix()

    # For analysis
    Path(paths['analysis_path']).mkdir(exist_ok=True, parents=True)
    paths['analysis_python_in'] = Path(paths['template_path'], 'python',
                                       'analysis.py').as_posix()
    paths['analysis_python_out'] = Path(paths['analysis_path'],
                                        'analysis.py').as_posix()
    paths['analysis_sub_in'] = Path(paths['template_path'], 'sub', 
                                    'analysis.qs').as_posix()
    paths['analysis_sub_out'] = Path(paths['analysis_path'],
                                     'analysis.qs').as_posix()
    out_dict['paths'] = paths

    # Log folder
    Path(paths['log_path']).mkdir(exist_ok=True, parents=True)

    '''Add extra parameters for job'''
    n_col = [row['n'] for row in out_dict['descriptors']]
    sampling = out_dict['job']['sampling']
    out_dict['job']['n_jobs'] = get_tot_job_num(n_col=n_col, sampling=sampling)
    out_dict['job']['desc_shape'] = get_descriptors_shape(n_col=n_col,
                                                          sampling=sampling)
    out_dict['job']['n_concurrent'] = int(np.min([out_dict['job']['n_jobs'],
                                                  out_dict['job']['n_concurrent']]))
    return out_dict

def get_rel_paths(paths, start, suffix=''):
    """Helper method to return relative paths

    Parameters
    ----------
        paths : dict
            Paths to convert to relative paths. The key is a str identifying
            what the path represents and the value is the path
        start : str
            Path to start from
        suffix : str, optional
            Suffix to add onto end. Default is ''
    Returns
    -------
        rel_paths : dict
            Relative paths.
    """
    rel_paths = {}
    for name, path in paths.items():
        rel_path = os.path.relpath(path=path, start=start)
        # Standardize Path
        rel_path = Path(rel_path).as_posix()
        rel_paths['{}{}'.format(name, suffix)] = rel_path
    return rel_paths
