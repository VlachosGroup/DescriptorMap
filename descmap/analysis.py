import os

import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go

from vunits import constants as c
from vunits.quantity import Quantity
from pmutt.cantera.units import Units

def get_reactor_inputs(reactor_inputs, phases, units):
    """Extracts the reactor inputs
    
    Parameters
    ----------
        reactor_inputs : dict
            Conditions of the reactor. Keys represent reactor conditions and
            values represent quantity
        phases : list of dict
            Phases to consider for site density calculation. At least one
            dictionary is expected to have 'site_density' as the key.
        units : dict
            Units corresponding to interested quantities. If an entry is not
            provided, uses `pMuTT's default units for OpenMKM`_.
    Returns
    -------
        reactor : dict of `Quantity`_ objects
            Reactor conditions containing the keys:
            
            - 'Q'
            - 'P'
            - 'V'
            - 'T'
            - 'cat_abyv'
            - 'site_den'
            - 'mol_catalyst'

    .. _`pMuTT's default units for OpenMKM`: https://vlachosgroup.github.io/pMuTT/api/kinetic_models/cantera_units/pmutt.cantera.units.Units.html#pmutt.cantera.units.Units
    .. _`Quantity`: https://vlachosgroup.github.io/vunits/api/quantity/quantity.html#quantity-class
    """
    # Create unit quantities
    default_units = Units().to_CTI_dict()
    unit_qty = {}
    for unit_type in default_units.keys():
        unit_type = unit_type.replace('_unit', '')
        try:
            unit_qty[unit_type] = Quantity.from_units(units=units[unit_type])
        except KeyError:
            unit_qty[unit_type] = \
                    Quantity.from_units(units=default_units[unit_type])
    # Create output dictionary
    reactor = {
        'Q': reactor_inputs['flow_rate']*unit_qty['length']**3/unit_qty['time'],
        'P': reactor_inputs['P']*unit_qty['pressure'],
        'V': reactor_inputs['V']*unit_qty['length']**3,
        'T': Quantity(mag=reactor_inputs['T'], K=1.),
        'cat_abyv': reactor_inputs['cat_abyv']/unit_qty['length']
    }
    # Calculate site density
    site_den = 0.
    for phase in phases:
        try:
            site_den += phase['site_density']
        except KeyError:
            pass
    reactor['site_den'] = site_den*unit_qty['quantity']/unit_qty['length']**2
    reactor['mol_catalyst'] = reactor['site_den']*reactor['V']*reactor['cat_abyv']
    return reactor

def get_flow_rates(reactor, mol_frac_in, mol_frac_out, mass_frac_in,
                   mass_frac_out, mol_weight):
    """Calculates the mass and mol flow rates in and out of the reactor
    
    Parameters
    ----------
        reactor : dict of `Quantity`_ objects.
            Reactor condition inputs
        mol_frac_in : (M, N) `pd.DataFrame`_ obj
            Mole fractions in. M and N represents number of runs and species
            respectively.
        mol_frac_out : (M, N) `pd.DataFrame`_ obj
            Mole fractions out. M and N represents number of runs and species
            respectively.
        mass_frac_out : (M, N) `pd.DataFrame`_ obj
            Mass fractions in. M and N represents number of runs and species
            respectively.
        mass_frac_out : (M, N) `pd.DataFrame`_ obj
            Mass fractions out. M and N represents number of runs and species
            respectively.
        mol_weight : `pd.Series`_ obj
            Mole weight in g/mol
    Returns
    -------
        mol_flow_rates : dict
            Flow rates containing the keys:

            - in
            - in_tot
            - out
            - out_tot

        mass_flow_rates : dict
            Mass flow rates containing the keys:

            - in
            - out
            - tot
        

    .. _`pd.DataFrame` : https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html#pandas-dataframe
    .. _`pd.Series`: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.Series.html#pandas.Series
    """

    '''Calculate mole flow in'''
    mol_flow_rates = {}
    mass_flow_rates = {}
    mol_flow_rates['in_tot'] = reactor['P']*reactor['Q']/reactor['T']/c.R
    mol_flow_rates['in'] = mol_flow_rates['in_tot']('mol/s')*mol_frac_in # mol/s
    '''Calulate mass flow in/out'''
    mass_flow_rates['tot'] = mol_flow_rates['in'].mul(mol_weight).sum(axis=1, skipna=True) # g/s
    mass_flow_rates['in'] = mass_frac_in.mul(mass_flow_rates['tot'], axis=0)
    mass_flow_rates['out'] = mass_frac_out.mul(mass_flow_rates['tot'], axis=0)
    '''Calculate mole flow out'''
    mol_flow_rates['out_tot'] = mass_flow_rates['out'].mul(1./mol_weight).sum(axis=1, skipna=True) # mol/s
    mol_flow_rates['out'] = mol_frac_out.mul(mol_flow_rates['out_tot'], axis=0) # mol/s
    return (mol_flow_rates, mass_flow_rates)

def plot_contour(path_out, x, y, z, title, x_label, y_label):
    """Creates a Contour plot using Plotly
    
    Parameters
    ----------
        path_out : str
            Name of file (ending with .html extension)
        x : list
            x data to plot
        y : list
            y data to plot
        z : list
            z data to plot
        title : str
            Title of plot
        x_label : str
            x axis label
        y_label : str
            y axis label
    """
    layout={'title': {'text': title},
            'xaxis': {'title': x_label},
            'yaxis': {'title': y_label}}
    fig = go.Figure(go.Contour(x=x, y=y, z=z), layout=layout)
    fig.write_html(path_out)

def plot_density(path_out, jobs_data, desc_labels, conv_data, selec_data,
                 yield_data, reactant_name, product_name):
    """Plots a density map normalized to x axis slices

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
        kwargs : keyword arguments
            Keyword arguments for `numpy.histogram2d`_
    .. _`numpy.histogram2d`: https://numpy.org/doc/1.18/reference/generated/numpy.histogram2d.html
    """
    fig = make_subplots(rows=len(desc_labels), cols=3)
    for i, desc_label in enumerate(desc_labels, start=1):
        # Get x data
        x = jobs_data[desc_label]

        # Conversion plot
        conv_trace_name = 'Conv<br />vs.<br />{}'.format(desc_label)
        conv_trace = go.Histogram2dContour(x=x, y=conv_data,
                                           name=conv_trace_name,
                                           showscale=False)
        fig.add_trace(conv_trace, row=i, col=1)
        fig.update_xaxes(title_text='{} (eV)'.format(desc_label),
                         row=i, col=1)
        fig.update_yaxes(title_text='Conversion of {}'.format(reactant_name),
                         row=i, col=1)

        # Selectivity plot
        selec_trace_name = 'Selec<br />vs.<br />{}'.format(desc_label)
        selec_trace = go.Histogram2dContour(x=x, y=selec_data,
                                            name=selec_trace_name,
                                            showscale=False)        
        fig.add_trace(selec_trace, row=i, col=2)
        fig.update_xaxes(title_text='{} (eV)'.format(desc_label),
                         row=i, col=2)
        fig.update_yaxes(title_text='Selectivity of {}'.format(product_name),
                         row=i, col=2)

        # Yield plot
        yield_trace_name = 'Yield<br />vs.<br />{}'.format(desc_label)
        yield_trace = go.Histogram2dContour(x=x, y=yield_data,
                                            name=yield_trace_name,
                                            showscale=False)        
        fig.add_trace(yield_trace, row=i, col=3)
        fig.update_xaxes(title_text='{} (eV)'.format(desc_label),
                         row=i, col=3)
        fig.update_yaxes(title_text='Yield of {}'.format(product_name),
                         row=i, col=3)

    fig.write_html(path_out)

def get_fractions(paths, omkm_path):
    """Returns the coverages, gas fractions and mole fractions associated with
    OpenMKM runs.

    Parameters
    ----------
        paths : str
            Path to OpenMKM subfolders
        omkm_path : str
            Path to OpenMKM parent folder
    Returns
    -------
        fractions : dict
            Fractions from reactor with the following keys:

            - cov_in
            - covs_out
            - mass_in
            - mass_out
            - mol_in
            - mol_out

    """
    rel_job_paths = [omkm_path.joinpath(path).as_posix() for path in paths]
    covs_in, covs_out = get_single_fractions(paths=rel_job_paths,
                                             filename='surf_cov_ss.csv',
                                             index=paths)
    mass_frac_in, mass_frac_out = get_single_fractions(paths=rel_job_paths,
                                                       filename='gas_mass_ss.csv',
                                                       index=paths)
    mol_frac_in, mol_frac_out = get_single_fractions(paths=rel_job_paths,
                                                     filename='gas_mole_ss.csv',
                                                     index=paths)
    fractions = {
        'covs_in': covs_in,
        'covs_out': covs_out,
        'mass_in': mass_frac_in,
        'mass_out': mass_frac_out,
        'mol_in': mol_frac_in,
        'mol_out': mol_frac_out
    }
    return fractions

def get_single_fractions(paths, filename, index):
    initial_fractions = []
    final_fractions = []
    for path in paths:
        full_path = os.path.join(path, filename)
        fractions = pd.read_csv(full_path, header=0)
        del fractions['t(s)']
        initial_fractions.append(fractions.iloc[0])
        final_fractions.append(fractions.iloc[-1])
    initial_fractions = pd.DataFrame(initial_fractions, index=index)
    final_fractions = pd.DataFrame(final_fractions, index=index)
    return initial_fractions, final_fractions

