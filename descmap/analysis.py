import os

import numpy as np
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
    default_units = Units().to_cti_dict()
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

def plot_contour(path_out, x, y, z, title, x_label, y_label, hover_label,
                 x_lit=None, y_lit=None, lit_label=None, zmin_cutoff=None,
                 zmax_cutoff=None):
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
        zmin_cutoff : float, optional
            Minimum cutoff z value. The minimum of ``z`` will be used if it is
            higher than ``zmin_cutoff`` or if ``zmin_cutoff`` is not specified.
        zmax_cutoff : float, optional
            Maximum cutoff z value. The maximum of ``z`` will be used if it is
            higher than ``zmax_cutoff`` of ir ``zmax_cutoff`` is not specified.
    """
    layout={'title': {'text': title},
            'xaxis': {'title': x_label,
                      'tickformat': '0.2f',
                      'ticks': 'outside',
                      'mirror': True,
                      'showline': True},
            'yaxis': {'title': y_label,
                      'tickformat': '0.2f',
                      'ticks': 'outside',
                      'mirror': True,
                      'showline': True,
                      'linewidth': 2.},
            'legend': {'x': 0., 'y': 1}}

    # Process zmin_cutoff
    temp_zmin = np.floor(z.min())
    if zmin_cutoff is None:
        zmin_cutoff = temp_zmin
    elif zmin_cutoff < temp_zmin:
        zmin_cutoff = temp_zmin
    # Process zmax_cutoff
    temp_zmax = np.ceil(z.max())
    if zmax_cutoff is None:
        zmax_cutoff = temp_zmax
    elif zmax_cutoff > temp_zmax:
        zmax_cutoff = temp_zmax

    fig = go.Figure(go.Contour(x=x, y=y, z=z, hovertext=hover_label,
                               zmin=zmin_cutoff, zmax=zmax_cutoff),
                    layout=layout) 


    if x_lit is not None:
        fig.add_trace(go.Scatter(x=x_lit, y=y_lit, hovertext=lit_label,
                                 mode='markers+text',
                                 text=lit_label,
                                 name='Literature',
                                 textposition="bottom center",
                                 legendgroup='Literature',
                                 showlegend=True,
                                 textfont={'size': 15,
                                           'color': 'black'},
                                 marker={'size': 10,
                                         'color': 'white',
                                         'line': {'width': 2,
                                                  'color': 'black'}}))
    # fig.update_layout(legend_orientation="h")
    fig.write_html(path_out)
    fig.write_image(path_out.replace('html', 'png'),
                    scale=10, width=6, height=8)

def plot_density(path_out, jobs_data, desc_labels, conv_data, selec_data,
                 yield_data, reactant_name, product_name, hover_label,
                 lit_data=None):
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
    fig = make_subplots(cols=len(desc_labels), rows=3)
    for i, desc_label in enumerate(desc_labels, start=1):
        # Get x data
        x = jobs_data[desc_label]
        QoI_data = (conv_data, selec_data, yield_data)
        QoI_labels = ('{} Conv'.format(reactant_name),
                      '{} Selectivity'.format(product_name),
                      '{} Yield'.format(product_name))

        for j, (y_data, y_label) in enumerate(zip(QoI_data, QoI_labels),
                                              start=1):
            # Since legends are grouped, only display one toggle option
            if i == 1 and j == 1:
                showlegend = True
            else:
                showlegend = False

            trace_name = '{}<br />vs.<br />{}'.format(y_label, desc_label)
            hist_trace = go.Histogram2dContour(x=x, y=y_data,
                                               name=trace_name,
                                               showscale=False,
                                               colorscale='Blues')
            scatter_trace = go.Scatter(x=x, y=y_data, mode='markers',
                                       hovertext=hover_label,
                                       name='MKM Data',
                                       marker={'size': 3,
                                               'color': 'rgba(0., 0., 0., 0.5)'},
                                       legendgroup='MKM Data',
                                       showlegend=showlegend)
            fig.add_traces(data=[hist_trace, scatter_trace],
                           cols=(i, i), rows=(j, j))

            # Update axes
            fig.update_xaxes(title_text='{} (eV)'.format(desc_label),
                             col=i, row=j, ticks='outside', tickformat='.2f',
                             showline=True, mirror=True, linewidth=1.,
                             linecolor='black')
            fig.update_yaxes(tickformat='.2f', showline=True, mirror=True,
                             row=j, col=i, linewidth=1., linecolor='black')

        # Add Y Label
        if i == 1:
            for j, y_label in enumerate(QoI_labels, start=1):
                fig.update_yaxes(title_text=y_label,
                                 col=i, row=j, ticks='outside')

        # Add literature data
        # if lit_data is not None:
        #     x_lit = lit_data[desc_label]
        #     for j in range(1, 4):
        #         for name, x_lit_point in x_lit.iteritems():
        #             fig.add_shape({'type': 'line',
        #                            'xref': 'x',
        #                            'yref': 'paper',
        #                            'x0': x_lit_point,
        #                            'x1': x_lit_point,
        #                            'y0': ,
        #                            'y1': 1,
        #                            }, col=i, row=j)

    fig.write_html(path_out)
    fig.write_image(path_out.replace('html', 'png'),
                    scale=10, width=8, height=14)

def plot_1d_volcano(path_out, x, y, title, x_label, y_label, hover_label,
                    x_lit=None, lit_label=None, ymin_cutoff=None,
                    ymax_cutoff=None):
    """Creates a Contour plot using Plotly
    
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
    """
    layout={'title': {'text': title},
            'xaxis': {'title': x_label,
                      'tickformat': '0.2f',
                      'ticks': 'outside',
                      'mirror': True,
                      'showline': True},
            'yaxis': {'title': y_label,
                      'tickformat': '0.2f',
                      'ticks': 'outside',
                      'mirror': True,
                      'showline': True,
                      'linewidth': 2.},
            'legend': {'x': 0., 'y': 1}}

    # Process zmin_cutoff
    temp_ymin = np.floor(y.min())
    if ymin_cutoff is None:
        ymin_cutoff = temp_ymin
    elif ymin_cutoff < temp_ymin:
        ymin_cutoff = temp_ymin
    # Process zmax_cutoff
    temp_ymax = np.ceil(y.max())
    if ymax_cutoff is None:
        ymax_cutoff = temp_ymax
    elif ymax_cutoff > temp_ymax:
        ymax_cutoff = temp_ymax

    fig = go.Figure(go.Scatter(x=x, y=y, hovertext=hover_label, mode='lines'),
                    layout=layout)
    fig.update_yaxes(range=[ymin_cutoff, ymax_cutoff])

    # if x_lit is not None:
    #     fig.add_trace(go.Scatter(x=x_lit, y=y_lit, hovertext=lit_label,
    #                              mode='markers+text',
    #                              text=lit_label,
    #                              name='Literature',
    #                              textposition="bottom center",
    #                              legendgroup='Literature',
    #                              showlegend=True,
    #                              textfont={'size': 15,
    #                                        'color': 'black'},
    #                              marker={'size': 10,
    #                                      'color': 'white',
    #                                      'line': {'width': 2,
    #                                               'color': 'black'}}))
    # fig.update_layout(legend_orientation="h")
    fig.write_html(path_out)
    fig.write_image(path_out.replace('html', 'png'),
                    scale=10, width=6, height=8)

    

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
        # Delete time and space coordinates
        for field in ('t(s)', 'z(m)'):
            try:
                del fractions[field]
            except KeyError:
                pass        
        initial_fractions.append(fractions.iloc[0])
        final_fractions.append(fractions.iloc[-1])
    initial_fractions = pd.DataFrame(initial_fractions, index=index)
    final_fractions = pd.DataFrame(final_fractions, index=index)
    return initial_fractions, final_fractions

