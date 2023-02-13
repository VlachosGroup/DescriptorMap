"""Functionality related to MKM kinetics analysis using electronic descriptors"""
import os
import numpy as np
import pandas as pd
from scipy import interpolate
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from matplotlib.patches import Ellipse
from matplotlib import transforms
from scipy.stats import t as student_t
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
                 zmax_cutoff=None, x_scaling=None, y_scaling=None,
                 y_scaling_upper=None, y_scaling_lower=None):
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
    zmin_cutoff, zmax_cutoff = _get_cutoff(data=z,
                                           min_cutoff=zmin_cutoff,
                                           max_cutoff=zmax_cutoff)
    fig = go.Figure(go.Contour(x=x, y=y, z=z, hovertext=hover_label,
                               zmin=zmin_cutoff, zmax=zmax_cutoff,
                               name=title,
                               showlegend=False,
                               colorscale="Viridis"),
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
                                           'color': 'white'},
                                 marker={'size': 10,
                                         'color': 'white',
                                         'line': {'width': 2,
                                                  'color': 'black'}}))
    # Add scaling line between descriptors
    if y_scaling is not None:
        # Add line to plot
        fig.add_trace(go.Scatter(x=x_scaling, y=y_scaling,
                                 mode='lines',
                                 name='Descriptor Correlation',
                                 legendgroup='Intervals',
                                 showlegend=False,
                                 line = {'width': 1,
                                         'color': 'white'}))

    minor_scaling_kwargs = {'mode': 'lines',
                            'name': 'Prediction Intervals',
                            'legendgroup': 'Intervals',
                            'line': {'width': 1,
                                     'color': 'gray'}}
    if y_scaling_upper is not None:
        fig.add_trace(go.Scatter(x=x_scaling, y=y_scaling_upper,
                                 showlegend=True,
                                 **minor_scaling_kwargs))

    if y_scaling_lower is not None:
        fig.add_trace(go.Scatter(x=x_scaling, y=y_scaling_lower,
                                 showlegend=False,
                                 **minor_scaling_kwargs))

        # x_linspace = np.linspace(np.min(x_lit), np.max(x_lit))
        # y_linspace = lit_lin_model.predict(X=x_linspace.reshape(-1, 1)).flatten()

        # upper_interval, lower_interval = get_prediction_interval(y_data=y_lit_model.flatten(),
        #                                                          y_model=y_lit_pred,
        #                                                          x_data=x_lit,
        #                                                          x_linspace=x_linspace,
        #                                                          alpha=0.99)
        # upper_interval = upper_interval + y_linspace
        # lower_interval = lower_interval + y_linspace


        # x_ell, y_ell = get_ellipse(x=x_lit.values,
        #                            y=y_lit.values,
        #                            x_model=x_lit_model,
        #                            y_model=y_lit_pred,
        #                            slope=lit_lin_model.coef_[0][0],
        #                            n=50)
        # fig.add_trace(go.Scatter(x=x_ell,
        #                          y=y_ell,
        #                          mode='lines',
        #                          name='Ellipse'))



    # fig.update_layout(legend_orientation="h")
    fig.write_html(path_out)
    fig.write_image(path_out.replace('html', 'png'),
                    scale=10, width=6, height=8)
    fig.write_image(path_out.replace('html', 'svg'),
                    scale=10, width=6, height=8)

def plot_density(path_out, jobs_data, desc_labels, conv_data, selec_data,
                 yield_data, reactant_name, product_name, hover_label,
                 lit_data=None, design_space_mask=None,
                 conv_min_cutoff=0., conv_max_cutoff=1.,
                 selec_min_cutoff=0., selec_max_cutoff=1.,
                 yield_min_cutoff=0., yield_max_cutoff=1.):
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
        ymin_cutoffs = (conv_min_cutoff, selec_min_cutoff, yield_min_cutoff)
        ymax_cutoffs = (conv_max_cutoff, selec_max_cutoff, yield_max_cutoff)

        for j, (y_data, y_label, ymin_cutoff, ymax_cutoff) in \
            enumerate(zip(QoI_data, QoI_labels, ymin_cutoffs, ymax_cutoffs),
                      start=1):
            # Since legends are grouped, only display one toggle option
            if i == 1 and j == 1:
                showlegend = True
            else:
                showlegend = False

            ymin_cutoff, ymax_cutoff = _get_cutoff(data=y_data,
                                                   min_cutoff=ymin_cutoff,
                                                   max_cutoff=ymax_cutoff)

            trace_name = '{}<br />vs.<br />{}'.format(y_label, desc_label)
            # hist_trace = go.Histogram2dContour(x=x, y=y_data,
            #                                    name=trace_name,
            #                                    showscale=False,
            #                                    colorscale='Blues')
            scatter_trace = go.Scatter(x=x, y=y_data, mode='markers',
                                       hovertext=hover_label,
                                       name='MKM Data',
                                       marker={'size': 3,
                                               'color': 'rgba(0., 0., 0., 0.5)'},
                                       legendgroup='MKM Data',
                                       showlegend=showlegend)

            '''Find Pareto Optimal for Design Space'''
            # Get unique values
            x_unique = sorted(np.unique(x))
            y_largest = []
            y_smallest = []
            y_largest_design = []
            y_smallest_design = []
            for x_val in x_unique:
                # Find indices corresponding to x value
                k = np.where(x_val == x)[0]
                y_largest_val = np.max(y_data[k])
                y_smallest_val = np.min(y_data[k])
                y_largest.append(y_largest_val)
                y_smallest.append(y_smallest_val)

                # Find indices inside design space 
                k = np.where(np.logical_and(x_val == x,
                                            np.array(design_space_mask) == True))[0]
                y_largest_design_val = np.max(y_data[k])
                y_smallest_design_val = np.min(y_data[k])
                y_largest_design.append(y_largest_design_val)
                y_smallest_design.append(y_smallest_design_val)

            pareto_largest_trace = go.Scatter(x=x_unique,
                                              y=y_largest,
                                              fill=None,
                                              mode='lines',
                                              line_color='#1f77b4',
                                              showlegend=False,
                                              legendgroup='QoI Range (All)',
                                              name='QoI Range (All)')

            pareto_smallest_trace = go.Scatter(x=x_unique,
                                               y=y_smallest,
                                               fill='tonexty',
                                               mode='lines',
                                               line_color='#1f77b4',
                                               showlegend=showlegend,
                                               legendgroup='QoI Range (All)',
                                               name='QoI Range (All)')

            pareto_largest_design_trace = go.Scatter(x=x_unique,
                                                     y=y_largest_design,
                                                     fill=None,
                                                     mode='lines',
                                                     legendgroup='QoI Range (Prediction Intervals)',
                                                     line_color='#ff7f0e',
                                                     showlegend=False,
                                                     name='QoI Range (Prediction Intervals)')
            pareto_smallest_design_trace = go.Scatter(x=x_unique,
                                                      y=y_smallest_design,
                                                      fill='tonexty',
                                                      mode='lines',
                                                      legendgroup='QoI Range (Prediction Intervals)',
                                                      line_color='#ff7f0e',
                                                      showlegend=showlegend,
                                                      name='QoI Range (Prediction Intervals)')

            fig.add_traces(data=[scatter_trace, pareto_largest_trace,
                                 pareto_smallest_trace, pareto_largest_design_trace,
                                 pareto_smallest_design_trace],
                           cols=(i, i, i, i, i), rows=(j, j, j, j, j))

            # Update axes
            fig.update_xaxes(title_text='{} (eV)'.format(desc_label),
                             col=i, row=j, ticks='outside', tickformat='.2f',
                             showline=True, mirror=True, linewidth=1.,
                             linecolor='black')
            fig.update_yaxes(tickformat='.2f', showline=True, mirror=True,
                             row=j, col=i, linewidth=1., linecolor='black',
                             range=[ymin_cutoff, ymax_cutoff])

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
    # fig.write_image(path_out.replace('html', 'png'),
    #                 scale=10, width=8, height=14)
    # fig.write_image(path_out.replace('html', 'svg'),
    #                 scale=10, width=8, height=14)
    fig.update_layout(width=500.*len(desc_labels), height=800.)
    fig.write_image(path_out.replace('html', 'svg'))

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

    ymin_cutoff, ymax_cutoff = _get_cutoff(data=y,
                                           min_cutoff=ymin_cutoff,
                                           max_cutoff=ymax_cutoff)
    fig = go.Figure(go.Scatter(x=x, y=y, hovertext=hover_label, mode='lines'),
                    layout=layout)
    fig.update_yaxes(range=[ymin_cutoff, ymax_cutoff])

    if x_lit is not None:
        interp_fn = interpolate.interp1d(x=x, y=y)
        y_lit = interp_fn(x_lit)
        fig.add_trace(go.Scatter(x=x_lit, y=y_lit, hovertext=lit_label,
                                 mode='markers+text',
                                 text=lit_label,
                                 name='Literature',
                                 textposition="top right",
                                 legendgroup='Literature',
                                 showlegend=True,
                                 textfont={'size': 15,
                                           'color': 'black'},
                                 marker={'size': 10,
                                         'color': 'white',
                                         'line': {'width': 2,
                                                  'color': 'black'}}))
        fig.update_layout(legend_orientation="h")
    fig.write_html(path_out)
    fig.write_image(path_out.replace('html', 'png'),
                    scale=10, width=6, height=8)
    fig.write_image(path_out.replace('html', 'svg'),
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

def get_ellipse(x, y, slope, x_model, y_model, n_std=1., n=50):
    """Calculates the x, y coordinates of an ellipse based on the data

    Parameters
    ----------
        x : (N,) nd.ndarray
            X coordinates
        y : (N,) nd.ndarray
            Y coordinates
        n : int
            Number of data points to fit ellipse
    Returns
    -------
        x_ell : (M,) nd.ndarray
            X coordinates corresponding to ellipse perimeter
        y_ell : (M,) np.ndarray
            Y coordinates corresponding to ellipse perimeter
    """
    # Transform data
    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    x_ell_radius = np.sqrt(1 + pearson)
    y_ell_radius = np.sqrt(1 - pearson)

    scale_factor = np.sqrt((x_model[0]-x_model[1])**2 + (y_model[0]-y_model[1])**2)/(x_ell_radius*2)


    x_scale = np.sqrt(cov[0, 0]) * n_std
    y_scale = np.sqrt(cov[1, 1]) * n_std

    t = np.linspace(0., 2.*np.pi)
    x_ell = x_ell_radius*np.cos(t)
    # x_ell = scale_factor*np.cos(t)
    y_ell = y_ell_radius*np.sin(t)
    # r_ell = np.sqrt((x_ell**2 + y_ell**2))

    # Scale data
    x_ell = scale_factor*x_ell
    y_ell = scale_factor*y_ell
    # x_ell = x_scale*x_ell
    # y_ell = y_scale*y_ell

    # Rotate data
    angle = np.arctan(slope)
    u_ell, v_ell = _rotate(theta=angle, x=x_ell, y=y_ell)
    # u_ell = x_ell*np.cos(angle) - y_ell*np.sin(angle)
    # v_ell = x_ell*np.sin(angle) + y_ell*np.cos(angle)
    # u_ell = x_ell
    # v_ell = y_ell

    # Translate data
    x_center = (x_model[0] + x_model[1])/2
    y_center = (y_model[0] + y_model[1])/2
    u_ell = u_ell + x_center
    v_ell = v_ell + y_center

    return (u_ell, v_ell)

def get_prediction_interval(y_data, y_model, x_data, x_linspace, y_linspace, alpha=0.95):

    n = len(y_data) - 2
    y_sum_sqr_err = np.sum(mean_squared_error(y_data,
                                              y_model,
                                              multioutput='raw_values'))
    stdev = np.sqrt(y_sum_sqr_err/(n - 2))

    x_mean = np.mean(x_data)
    x_var = np.sum([(x-x_mean)**2 for x in x_data])

    # Sample the student t distribution
    conf = 1 - alpha
    ppf_lookup = 1 - conf/2
    t_dist = student_t(df=n-2)
    z_score = t_dist.ppf(ppf_lookup)
    pred_interval = z_score*stdev*np.sqrt((1. + 1./n + (x_linspace - x_mean)**2/x_var))
    lower_pred_interval = y_linspace - pred_interval
    upper_pred_interval = y_linspace + pred_interval

    return (lower_pred_interval, upper_pred_interval)

def _rotate(theta, x, y):
    rot_mat = np.array([[np.cos(theta), -np.sin(theta)],
                        [np.sin(theta),  np.cos(theta)]])
    orig_coor = np.vstack((x, y))
    new_coor = np.dot(rot_mat, orig_coor)
    return new_coor[0, :], new_coor[1, :]

def _get_cutoff(data, min_cutoff=None, max_cutoff=None):
    """Helper method to calculate cutoff values.

    Parameters
    ----------
        data : pandas.Series object
            Data to evaluate
        min_cutoff : float, optional
            Minimum cutoff value. If not specified, the minimum axes value
            will be based on ``data``.
        max_cutoff : float, optional
            Maximum cutoff value. If not specified, the maximum axes value will
            be based on ``data``.
    Returns
    -------
        min_cutoff : float
            Minimum cutoff value
        max_cutoff : float
            Maximum cutoff value
    """
    # Process min cutoff
    temp_min = np.floor(data.min())
    if min_cutoff is None:
        min_cutoff = temp_min
    elif min_cutoff < temp_min:
        min_cutoff = temp_min
    # Process max_cutoff
    temp_max = np.ceil(data.max())
    if max_cutoff is None:
        max_cutoff = temp_max
    elif max_cutoff > temp_max:
        max_cutoff = temp_max
    return (min_cutoff, max_cutoff)