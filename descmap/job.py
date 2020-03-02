"""Functionality related to job submissions"""

import os
import platform
import subprocess

import numpy as np
from pmutt.io.excel import read_excel

from descmap.errors import raise_invalid_sampling_method

def get_descriptors_shape(n_col, sampling):
    """Calculate the shape of the descriptor space
    
    Parameters
    ----------
        n_col : list of int
            Each row corresponds to number of runs for each dimension.
        sampling : str
            Type of sampling. Accepted values are 'linear' and 'lhs'
    Returns
    -------
        descriptors_shape : tuple
            Number of points in each dimension
    """
    if sampling == 'lhs':
        descriptors_shape = [int(np.max(n_col))]*len(n_col)
    elif sampling == 'linear':
        descriptors_shape = n_col
    else:
        raise_invalid_sampling_method(sampling)
    return descriptors_shape
    
def get_tot_job_num(n_col, sampling):
    """Calculates the total number of jobs
    
    Parameters
    ----------
        n_col : list of int
            Each row corresponds to number of runs for each dimension.
        sampling : str
            Type of sampling. Accepted values are 'linear' and 'lhs'
    Returns
    -------
        n_jobs : int
            Total number of jobs
    Raises
    ------
        ValueError:
            Raised if ``sampling`` is not supported.
    """
    if sampling.lower() == 'lhs':
        n_jobs = int(np.amax(n_col))
    elif sampling.lower() == 'linear':
        n_jobs = int(np.prod(n_col))
    else:
        raise_invalid_sampling_method(sampling)
    return n_jobs

def submit_job(cmd, test_run=False, verbose=True):
    """Submits the job using subprocess.Popen
    
    Parameters
    ----------
        cmd : list of str
            Command to run. Normally begins with 'sbatch' and ends with the
            path to the 'qs' script
    Returns
    -------
        id : str
            ID of submitted job
    """
    if test_run:
        proc_id = '0000'
        verbose_msg = '\tTest run. cmd: {}'.format(' '.join(cmd))
    elif platform.system() == 'Windows':
        proc_id = '0000'
        verbose_msg = '\tWindows detected: cmd {}'.format(' '.join(cmd))
    else:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        proc_id = proc.stdout.read().decode('utf-8').strip()
        verbose_msg = '\tSubmitted job. Job ID: {}'.format(proc_id)

    # Print message if requested
    if verbose:
        print(verbose_msg)
    return proc_id
