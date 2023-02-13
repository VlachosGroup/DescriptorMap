"""Functionality related to MKM kinetics analysis using geometric descriptors"""
import os
import re
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path
import xlsxwriter
import pmutt.constants as cons

# Import packages for plotting
from plotly.subplots import make_subplots
import chart_studio.plotly as py
import plotly.graph_objects as go
import plotly.express as px
import matplotlib
import matplotlib.pyplot as plt

'''Get QoI from Chemkin outputs'''
def get_tof(main_path, jobs_name):
    """Extract TOF from Chemkin outputs
    
    Parameters
    ----------
        main_path : str
            Path to descmap main folder (MKM parent folder)
        jobs_name : list
            List of descriptor names 

    Returns
    -------
        tof_list : list of TOF extracted from each MKM model
    """
    n_jobs = len(jobs_name)
    tof_list = []
    
    # Read TOF from MKM tube_conv.out file
    for i in range(n_jobs):
        output_path = os.path.join(main_path, "omkm", "{}".format(i), "OUT.d")
        os.chdir(output_path)
        
        with open('tube_conv.out', 'r') as fid:
            tof_file = fid.read()
        tof_line = tof_file.splitlines()[1:]
        tof = []
        if len(tof_line) == 1:
            for line in tof_line:
                data = re.split(' +', line)
            tof_list.append(abs(float(data[4])))
        else:
            for line in tof_line:
                data = re.split(' +', line)
                tof.append(abs(float(data[4])))
            tof_list.append(tof)
    return tof_list

def get_temperature(main_path, jobs_name):
    """Extract temperatures from Chemkin outputs
    
    Parameters
    ----------
        main_path : str
            Path to descmap main folder (MKM parent folder)
        jobs_name : list
            List of descriptor names 

    Returns
    -------
        temp_list : list of temperatures extracted from each MKM model
    """
    n_jobs = len(jobs_name)
    temp_list = []
    
    # Read reaction temperature from tube_conv.out file
    for i in range(n_jobs):
        output_path = os.path.join(main_path, "omkm", "{}".format(i), "OUT.d")
        os.chdir(output_path)
        
        with open('tube_conv.out', 'r') as fid:
            temp_file = fid.read()
        temp_line = temp_file.splitlines()[1:]
        temp = []
        if len(temp_line) == 1:
            for line in temp_line:
                data = re.split(' +', line)
            temp_list.append(abs(float(data[2])))
        else:
            for line in temp_line:
                data = re.split(' +', line)
                temp.append(abs(float(data[2])))
            temp_list.append(temp)        
    return temp_list

def get_conversion(main_path, jobs_name):
    """Extract reaction conversion from Chemkin outputs
    
    Parameters
    ----------
        main_path : str
            Path to descmap main folder (MKM parent folder)
        jobs_name : list
            List of descriptor names 

    Returns
    -------
        conv_list : list of conversions extracted from each MKM model
    """
    n_jobs = len(jobs_name)
    conv_list = []
    
    # Read reactant conversion from tube_conv.out file
    for i in range(n_jobs):
        output_path = os.path.join(main_path, "omkm", "{}".format(i), "OUT.d")
        os.chdir(output_path)
        
        with open('tube_conv.out', 'r') as fid:
            conv_file = fid.read()
        conv_line = conv_file.splitlines()[1:]
        conv = []
        if len(conv_line) == 1:
            for line in conv_line:
                data = re.split(' +', line)
            conv_list.append(abs(float(data[3])))
        else:
            for line in conv_line:
                data = re.split(' +', line)
                conv.append(abs(float(data[3])))
            conv_list.append(conv)          
    return conv_list
    
def get_coverage(main_path, jobs_name, inputs):
    """Extract coverages of surface species from Chemkin outputs
    
    Parameters
    ----------
        main_path : str
            Path to descmap main folder (MKM parent folder)
        jobs_name : list
            List of descriptor names 
        inputs: dict
            Inputs read from excel sheet. Keys correspond to sheet names.
    Returns
    -------
        cov_dict_list : list of dict
            Each dictionary corresponds to one MKM model. Keys represent species
            and values represent surface coverages.
    """
    n_jobs = len(jobs_name)
    
    # Read flow rate from inputs
    flow_rate = inputs['reactor']['flow_rate']    
    cov_dict_list = []
    
    # Read species coverages from tube_cov.out file
    for i in range(n_jobs):
        output_path = os.path.join(main_path, "omkm", "{}".format(i), "OUT.d")
        os.chdir(output_path)
        
        with open('tube_cov_ss.out', 'r') as fid:
            cov_file = fid.read()
        spe_line = cov_file.splitlines()[1]
        spe_name = re.split(' +', spe_line)
        # Remove the first two and last two elements to get the species list
        species = spe_name[2:-2]
        
        cov_lines = cov_file.splitlines()[2:]
        cov_dict = {}
        
        for j in range(len(species)):
            coverage = []
            volume = []
            res_time = []
            for line in cov_lines:
                vol = re.split(' +', line)[1]
                time = float(vol) / flow_rate
                cov_data = re.split(' +', line)[2:]
                
                coverage.append(float(cov_data[j]))
                volume.append(float(vol))
                res_time.append(time)
            
            spe_key = species[j]
            cov_dict['Time (s)'] = res_time
            cov_dict['Reactor length (cm^3)'] = volume
            cov_dict[spe_key] = coverage
        
        cov_dict_list.append(cov_dict)
    return cov_dict_list

def get_species(main_path):
    """Get names of all surfae species
    
    Parameters
    ----------
        main_path : str
            Path to descmap main folder (MKM parent folder)

    Returns
    -------
        species : list of surface species' names
    """
    output_path = os.path.join(main_path, "omkm", "0", "OUT.d")
    os.chdir(output_path)
    # Read surface species from tube_cov_ss.out file
    with open('tube_cov_ss.out', 'r') as fid:
        cov_file = fid.read()
    spe_line = cov_file.splitlines()[1]
    spe_name = re.split(' +', spe_line)
    # Remove the first two and last two elements to get the species list
    species = spe_name[2:-2]
    return species

def get_mole_fraction(main_path, jobs_name, reactant_name):
    """Extract reactant mole fraction flow in the reactor
    
    Parameters
    ----------
        main_path : str
            Path to descmap main folder (MKM parent folder)
        jobs_name : list
            List of descriptor names 
        reactant_name: str
            Name of reactant to be investigated
    Returns
    -------
        mole_list : list of mole fractions extracted from each MKM model
    """
    mole_list = []
    n_jobs = len(jobs_name)
    
    for i in range(n_jobs):
        in_path = os.path.join(main_path, "omkm", "{}".format(i), "INP.d")
        os.chdir(in_path)
        
        with open('tube_mole.inp', 'r') as in_fid:
            mole_file = in_fid.read()
        
        if reactant_name == 'CH4':
            ch4_line = mole_file.splitlines()[10]
            ch4_mole = np.array(re.split(' +', ch4_line)[1:]).astype(float)
            mole_list.append(ch4_mole)
        
        elif reactant_name == 'O2':
            o2_line = mole_file.splitlines()[11]
            o2_mole = np.array(re.split(' +', o2_line)[1:]).astype(float)
            mole_list.append(o2_mole)
        
        else:
            err_msg = 'Reactant {} does not exist'.format(reactant_name)
            raise ValueError(err_msg)
    return mole_list

def get_reaction_order(mole_list, tof_list, jobs_name):
    """Calculate reaction orders 
    
    Parameters
    ----------
        mole_list : list
            Extracted input mole fractions
        tof_list : list
            Extracted TOF from MKM outputs
        jobs_name: list
            List of descriptor names
    Returns
    -------
        reaction_order : list of reaction orders
    """
    n_jobs = len(jobs_name)
    reaction_order = []
    
    for i in range(n_jobs):
        mole = mole_list[i]
        tof = tof_list[i]
        
        # Calculate reaction orders
        x_data = np.log(mole)
        y_data = np.log(tof)
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_data, y_data)
        reaction_order.append(slope)
    return reaction_order
        
def get_Eapp(tof_list, temp_list, jobs_name):
    """Calculate apparent activation energies (kcal/mol)
    
    Parameters
    ----------
        tof_list : list
            Extracted TOF from MKM outputs
        temp_list : list
            Reaction temperatures
        jobs_name: list
            List of descriptor names
    Returns
    -------
        ea_list : list of calculated apparent activation energies
    """
    n_jobs = len(jobs_name)
    ea_list = []
    
    for i in range(n_jobs):
        temp = temp_list[i]
        tof = tof_list[i]
        
        # Calculate apparent activation energies
        x_data = [1/x for x in temp]
        y_data = np.log(tof)
        
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_data, y_data)
        app_ea = (-slope)*cons.R(units='kcal/mol/K')
        ea_list.append(app_ea)
    return ea_list

'''Plot QoI'''
def plot_1d(path_out, x, y, title, x_label, y_label, tickformat,
            ymin_cutoff = None, ymax_cutoff = None):
    """Creates a volcano plot using Plotly
    
    Parameters
    ----------
        path_out : str
            Name of file (ending with .html extension)
        x : list
            x data to plot
        y : list
            y data to plot
        title : str
            Title of plot
        x_label : str
            x axis label
        y_label : str
            y axis label
        ymin_cutoff : float, optional
            Minimum cutoff y value. The minimum of ``y`` will be used if it is
            higher than ``ymin_cutoff`` or if ``ymin_cutoff`` is not specified.
        ymax_cutoff : float, optional
            Maximum cutoff y value. The maximum of ``y`` will be used if it is
            higher than ``ymax_cutoff`` of ir ``ymax_cutoff`` is not specified.
        tickformat: str
            Specify number of decimals to keep in ticks
    """
    layout={'title': {'text': title},
            'xaxis': {'title': x_label,
                      'tickformat': tickformat,
                      'ticks': 'outside',
                      'mirror': True,
                      'showline': True},
            'yaxis': {'title': y_label,
                      'tickformat': tickformat,
                      'ticks': 'outside',
                      'mirror': True,
                      'showline': True,
                      'linewidth': 2.},
            'legend': {'x': 0., 'y': 1}}
    fig = go.Figure(go.Scatter(x=x, y=y, mode='markers+text',
                               marker=dict(color='red', size=12),
                               text=y, textposition='bottom center',
                               textfont=dict(size=12)),
                               layout=layout)
    fig.update_yaxes(range = [ymin_cutoff, ymax_cutoff])
    fig.write_html(path_out)
    fig.write_image(path_out.replace('html', 'png'),
                    scale=10, width=6, height=8)
    fig.write_image(path_out.replace('html', 'svg'),
                    scale=10, width=6, height=8)

def plot_ols(path_out, x, y, title, x_label, y_label, tickformat,
            ymin_cutoff = None, ymax_cutoff = None):
    """Creates ordinary least squares (OLS) regression plot using Plotly
    
    Parameters
    ----------
        path_out : str
            Name of file (ending with .html extension)
        x : list
            x data to plot
        y : list
            y data to plot
        title : str
            Title of plot
        x_label : str
            x axis label
        y_label : str
            y axis label
        ymin_cutoff : float, optional
            Minimum cutoff y value. The minimum of ``y`` will be used if it is
            higher than ``ymin_cutoff`` or if ``ymin_cutoff`` is not specified.
        ymax_cutoff : float, optional
            Maximum cutoff y value. The maximum of ``y`` will be used if it is
            higher than ``ymax_cutoff`` of ir ``ymax_cutoff`` is not specified.
        tickformat: str
            Specify number of decimals to keep in ticks
    """
    layout={'title': {'text': title},
            'xaxis': {'title': x_label,
                      'tickformat': tickformat,
                      'ticks': 'outside',
                      'mirror': True,
                      'showline': True},
            'yaxis': {'title': y_label,
                      'tickformat': tickformat,
                      'ticks': 'outside',
                      'mirror': True,
                      'showline': True,
                      'linewidth': 2.},
            'legend': {'x': 0., 'y': 1}}
    
    fig = go.Figure(px.scatter(x=x, y=y, trendline='ols', trendline_color_override='red'))
    fig.update_yaxes(range = [ymin_cutoff, ymax_cutoff])
    fig.update_layout(layout)
    fig.write_html(path_out)
    fig.write_image(path_out.replace('html', 'png'),
                    scale=10, width=6, height=8)
    fig.write_image(path_out.replace('html', 'svg'),
                    scale=10, width=6, height=8)

def plot_coverage(jobs_name, cov_dict_list, species_list, path):
    """Creates surface coverage plots using Matplotlib
    
    Parameters
    ----------
        jobs_name : list
            List of descriptor names
        cov_dict_list : list
            List of dictionaries containing surface species
        species_list : list
            List of surface species
        path : str
            Path to save output figures
    """
    n_jobs = len(jobs_name)
    
    for i in range(n_jobs):
        gcn = jobs_name[i]
        cov_dict = cov_dict_list[i]
        
        plt.figure()
        
        for j in species_list:
            spe_cov = cov_dict[j]
            r_time = cov_dict['Time (s)']
            
            if max(spe_cov) >= 1e-5:    # only plot species with coverage larger than 1e-5
                plt.plot(r_time, spe_cov, label = j)
        
        gcn_string = f'{gcn:.3f}'
        plt.legend(loc = 'best')
        plt.xlabel('Reaction Time (s)')
        plt.ylabel('Surface Coverage (ML)')
        plt.savefig(path + 'GCN{}_coverage.png'.format(gcn_string), transparent = False)

'''Write results to output excel'''
def write_model_output(main_path, desc_name, jobs_name, output_path, inputs):
    """Write Chemkin outputs to Excel
    
    Parameters
    ----------
        main_path : str
            Path to descmap main folder (MKM parent folder)
        desc_name : list
            List of descriptor names 
        jobs_name : list
            List of job_array names 
        output_path : str
            Path to save output excel file 
        inputs: dict
            Inputs read from excel sheet. Keys correspond to sheet names.
    """
    tof_list = get_tof(main_path, jobs_name)
    temp_list = get_temperature(main_path, jobs_name)
    conv_list = get_conversion(main_path, jobs_name)
    
    tof_data_dict = {'TOF': tof_list,
                     'Conversion': conv_list,
                     'Temperature': temp_list}
    tof_data_df = pd.DataFrame(data = tof_data_dict, index = jobs_name)
    tof_data_df.index.name = desc_name
    
    # Create a Pandas Excel writer using XlsxWriter as the engine
    writer = pd.ExcelWriter(output_path, engine = 'xlsxwriter')
    tof_data_df.to_excel(writer, sheet_name = 'Summary')
    
    n_jobs = len(jobs_name)
    cov_dict_list = get_coverage(main_path, jobs_name, inputs)
    for i in range(n_jobs):
        gcn = jobs_name[i]
        gcn_string = f'{gcn:.3f}'
        
        cov_dict = cov_dict_list[i]
        cov_df = pd.DataFrame(data = cov_dict)
        cov_df.to_excel(writer, sheet_name = 'GCN={}'.format(gcn_string))
    writer.save()
    
def write_rxn_order(main_path, jobs_name, reactant_name, output_path):
    """Write calculated reaction orders to Excel
    
    Parameters
    ----------
        main_path : str
            Path to descmap main folder (MKM parent folder)
        jobs_name : list
            List of job_array names 
        reactant_name : str
            Name of reactant to calculate the reaction order 
        output_path : str
            Path to save output excel file 
    """
    mole_list = get_mole_fraction(main_path, jobs_name, reactant_name)
    tof_list = get_tof(main_path, jobs_name)
    rxn_order_list = get_reaction_order(mole_list, tof_list, jobs_name)
    rounded_order = [round(x,2) for x in rxn_order_list]
    
    rxn_order_dict = {'Reaction Order': rounded_order}
    rxn_order_df = pd.DataFrame(data = rxn_order_dict, index = jobs_name)
    rxn_order_df.index.name = 'GCN'
    
    n_jobs = len(jobs_name)
    writer = pd.ExcelWriter(output_path, engine = 'xlsxwriter')
    rxn_order_df.to_excel(writer, sheet_name = 'Summary')
    
    for i in range(n_jobs):
        gcn = jobs_name[i]
        gcn_string = f'{gcn:.3f}'
        
        mole = mole_list[i]
        tof = tof_list[i]
        tof_dict = {'Mole Fraction': mole,
                    'TOF': tof}
        tof_df = pd.DataFrame(data = tof_dict)
        tof_df.to_excel(writer, sheet_name = 'GCN={}'.format(gcn_string))
    writer.save()

def write_Eapp(main_path, jobs_name, output_path):
    """Write calculated apparent activation energies to Excel
    
    Parameters
    ----------
        main_path : str
            Path to descmap main folder (MKM parent folder)
        jobs_name : list
            List of job_array names 
        output_path : str
            Path to save output excel file 
    """
    tof_list = get_tof(main_path, jobs_name)
    temp_list = get_temperature(main_path, jobs_name)
    Eapp_list = get_Eapp(tof_list, temp_list, jobs_name)
    
    Eapp_dict = {'Eapp (kcal/mol)': Eapp_list}
    Eapp_df = pd.DataFrame(data = Eapp_dict, index = jobs_name)
    Eapp_df.index.name = 'GCN'
    
    n_jobs = len(jobs_name)
    writer = pd.ExcelWriter(output_path, engine = 'xlsxwriter')
    Eapp_df.to_excel(writer, sheet_name = 'Summary')
    
    for i in range(n_jobs):
        gcn = jobs_name[i]
        gcn_string = f'{gcn:.3f}'
        
        tof = tof_list[i]
        temp = temp_list[i]
        Eapp_sub_dict = {'Temperature': temp,
                         'TOF': tof}
        Eapp_sub_df = pd.DataFrame(data = Eapp_sub_dict)
        Eapp_sub_df.to_excel(writer, sheet_name = 'GCN={}'.format(gcn_string))
    writer.save()
    