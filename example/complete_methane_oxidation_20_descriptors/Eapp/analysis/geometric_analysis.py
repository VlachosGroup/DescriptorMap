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
import plotly.express as px
import plotly.io as pio

from pmutt.io.excel import read_excel
from pmutt import get_molecular_weight
from vunits import constants as c
from vunits.quantity import Quantity
from descmap.analysis import plot_1d_volcano
from descmap.analysis_geometric import *

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
    
    Ea_path = os.path.join(analysis_path, "Eapp", "")
    Eapp_path = Path(os.path.join(analysis_path, "Eapp", ""))
    Eapp_path.mkdir(parents = True, exist_ok = True)
    
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
    temp_list = get_temperature(main_path, jobs_name)
    Eapp_list = get_Eapp(tof_list, temp_list, jobs_name)
    rounded_Eapp = [round(x,2) for x in Eapp_list] # round the reaction order for plotting
    
    '''Plot Quantities of Interest'''
    Eapp_title = 'Apparent Activation Energy (kcal/mol)'
    
    os.chdir(analysis_path) # Change current path to the analysis folder
    
    # Eapp vs GCN
    plot_1d(path_out = './Eapp.html',
            x = jobs_name,
            y = rounded_Eapp,
            title = Eapp_title,
            x_label = str(desc_name),
            y_label = 'Eapp (kcal/mol)',
            tickformat = '0.2f')
  
    # For each GCN: plot the regression line and save the corresponding figure in one folder
    for i in range(n_jobs):
        gcn = jobs_name[i]
        temp = temp_list[i]
        tof = tof_list[i]
        
        x_data = np.array([1/x for x in temp])
        y_data = np.log(tof)
        
        # Prepare plot arguments
        gcn_string = f'{gcn:.3f}'
        gcn_path = Path(str(Path(gcn_string))+".html") # Keep decimals to save figs for each GCN
        gcn_path = Eapp_path.joinpath(gcn_path)
        title = 'GCN={}'.format(gcn)
        
        plot_ols(path_out=gcn_path.as_posix(),
                 x = x_data,
                 y = y_data,
                 title = title,
                 x_label = '1/T (K-1)',
                 y_label = 'ln(TOF)',
                 tickformat = '0.4f')

    # Write outputs into excel
    write_Eapp(main_path, jobs_name, output_path)
    
