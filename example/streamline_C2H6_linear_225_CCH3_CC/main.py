import os
from os.path import relpath
import json
import shutil
import platform
import subprocess
from pathlib import Path

import numpy as np
from pmutt.io.excel import read_excel

from descmap.job import submit_job
from descmap.io import organize_excel_inputs, write_out_file

def process_setup(inputs):
    """Writes and submits job related to setup
    
    Parameters
    ----------
        inputs : dict
            Dictionary containing all the inputs from input.xlsx
    """
    Path('./setup').mkdir(exist_ok=True, parents=True)
    
    '''Writes Python script for setup'''
    write_out_file(in_path='./templates/setup.py',
                   out_path='./setup/setup.py',
                   fields=inputs['fields'])

    '''Writes submission script for setup'''
    write_out_file(in_path='./templates/setup.qs',
                   out_path='./setup/setup.qs',
                   fields=inputs['fields'])

    '''Submits setup job'''
    setup_cmd = ['sbatch', '--parsable', './setup/setup.qs']
    inputs['id']['setup'] = submit_job(setup_cmd)

def process_omkm(inputs):
    """Writes and submits job related to OpenMKM runs
    
    Parameters
    ----------
        inputs : dict
            Dictionary containing all the inputs from input.xlsx
    """
    Path('./omkm').mkdir(exist_ok=True, parents=True)
    '''Writes submission script for setup'''
    write_out_file(in_path='./templates/JA_omkm.qs',
                   out_path='./omkm/JA_omkm.qs',
                   fields=inputs['fields'])

    '''Submits OpenMKM job'''
    omkm_cmd = ['sbatch',
                '--parsable',
                '--dependency=afterok:{}'.format(inputs['id']['setup']),
                './omkm/JA_omkm.qs']
    inputs['id']['omkm'] = submit_job(omkm_cmd)

def process_analysis(inputs):
    """Writes and submits job related to analysis
    
    Parameters
    ----------
        inputs : dict
            Dictionary containing all the inputs from input.xlsx
    """
    Path('./analysis').mkdir(exist_ok=True, parents=True)
    '''Writes Python script for setup'''
    write_out_file(in_path='./templates/analysis.py',
                   out_path='./analysis/analysis.py',
                   fields=inputs['fields'])

    '''Writes submission script for analysis'''
    write_out_file(in_path='./templates/analysis.qs',
                   out_path='./analysis/analysis.qs',
                   fields=inputs['fields'])

    '''Submits analysis job'''
    analysis_cmd = ['sbatch',
                    '--parsable',
                    '--dependency=afterany:{}'.format(inputs['id']['omkm']),
                    './analysis/analysis.qs']
    inputs['id']['analysis'] = submit_job(analysis_cmd)

if __name__ == "__main__":
    '''Change working directory to script directory'''
    try:
        os.chdir(os.path.dirname(__file__))
    except FileNotFoundError:
        pass

    '''Read Excel sheet and organize the inputs'''
    inputs = organize_excel_inputs('inputs.xlsx')

    '''Setup'''
    process_setup(inputs=inputs)
    '''OMKM'''
    process_omkm(inputs=inputs)
    '''Analysis'''
    process_analysis(inputs=inputs)

    '''Write JSON file with inputs'''
    with open('./inputs.json', 'w', newline='\n') as f_ptr:
        json.dump(inputs, f_ptr, indent=2)
