# Scripts used as templates/analysis.py
import os
import json
import yaml
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotly.subplots import make_subplots
from sklearn.linear_model import LinearRegression
import plotly.graph_objects as go

from pmutt.io.excel import read_excel
from pmutt import get_molecular_weight
from vunits import constants as c
from vunits.quantity import Quantity
from descmap.analysis.plot import *
from descmap.analysis.chemkin import *

#%%
if __name__ == "__main__":
    # Change folder to script folder
    try:
        os.chdir(os.path.dirname(__file__))
    except FileNotFoundError:
        pass
       
    '''Paths'''
    main_path = os.path.dirname(os.getcwd())
    in_path = Path(r'../inputs.json')
    output_path = os.path.join(main_path, 'outputs.xlsx')
    analysis_path = os.path.join(main_path, "analysis", "")
    
    coverage_path = os.path.join(analysis_path, "cov", "")
    cov_path = Path(os.path.join(analysis_path, "cov", ""))
    cov_path.mkdir(parents = True, exist_ok = True)
    
    '''Read inputs'''
    with open(in_path, 'r') as f_ptr:
        inputs = json.load(f_ptr)

    '''Read relevant species'''
    reactant_name = inputs['analysis']['Reactant']
    product_name = inputs['analysis']['Product']
    desc_name = inputs['descriptors'][0]['name']
    n_desc = len(desc_name)

    '''Process simulations'''
    jobs_data = pd.read_excel('../job_summary.xlsx')
    jobs_data.set_index('Path', inplace = True)
    jobs_labels = ['Path: {}'.format(job_name) for job_name in jobs_data.index]
    jobs_name = jobs_data['GCN'].tolist()
    n_jobs = len(jobs_name)    

    '''Read MKM results'''
    tof_list = get_tof(main_path, jobs_name)
    conv_list = get_conversion(main_path, jobs_name)
    cov_dict_list = get_coverage(main_path, jobs_name, inputs)
    
    '''Get species list'''
    species_list = get_species(main_path)
    
    '''Plot Quantities of Interest'''
    conv_title = 'Conversion of {}'.format(reactant_name)
    tof_title = 'TOF of {}'.format(reactant_name)
    
    os.chdir(analysis_path) # Change current path to the analysis folder
    
    plot_1d_simple(path_out = './conv.html',
                   x = jobs_name,
                   y = conv_list,
                   title = conv_title,
                   x_label = str(desc_name),
                   y_label = 'Conversion (%)',
                   tickformat = '0.3f')
    
    plot_1d_simple(path_out = './tof.html',
                   x = jobs_name,
                   y = tof_list,
                   title = tof_title,
                   x_label = str(desc_name),
                   y_label = 'TOF (s^-1)',
                   tickformat = '0.4f')

    plot_coverage(jobs_name, cov_dict_list, species_list, coverage_path)
    
    write_model_output(main_path, desc_name, jobs_name, output_path, inputs)
    
