"""Functionality to plot interactive graphs"""
import os
from pathlib import Path
import chart_studio.plotly as py
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
import matplotlib
import matplotlib.pyplot as plt

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

def _get_cutoff(data, min_cutoff=None, max_cutoff=None):
    """Helper method to calculate cutoff values for plotting.

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