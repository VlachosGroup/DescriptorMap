import os
import json
import yaml
from pathlib import Path

import numpy as np
import pandas as pd
from plotly.subplots import make_subplots
from sklearn.linear_model import LinearRegression
import plotly.graph_objects as go

from pmutt.io.excel import read_excel
from pmutt import get_molecular_weight
from vunits import constants as c
from vunits.quantity import Quantity
from descmap.analysis.openmkm import *
from descmap.analysis.plot import *

if __name__ == "__main__":
    # Change folder to script folder
    try:
        os.chdir(os.path.dirname(__file__))
    except FileNotFoundError:
        pass

    '''Paths'''
    in_path = Path(r'../inputs.json')
    excel_path = Path(r'../inputs.xlsx')
    omkm_path = Path(r'../omkm')

    '''Read inputs'''
    with open(in_path, 'r') as f_ptr:
        inputs = json.load(f_ptr)

    '''Read relevant species'''
    reactant_name = inputs['analysis']['Reactant']
    product_name = inputs['analysis']['Product']
    desc_names = [descriptor['name'] for descriptor in inputs['descriptors']]
    n_desc = len(desc_names)
    desc_labels = ['E_{}'.format(name) for name in desc_names]
    labels = ['{} (eV)'.format(label) for label in desc_labels]

    '''Calculate molecular weights'''
    mol_weight = pd.Series({name: get_molecular_weight(comp) \
                            for name, comp in inputs['composition'].items()})

    '''Process simulations'''
    jobs_data = pd.read_excel('../job_summary.xlsx',
                              dtype={'Path': str})
    jobs_data.set_index('Path', inplace=True)
    jobs_labels = ['Path: {}'.format(job_name) for job_name in jobs_data.index]

    '''Get coverages, mole fractions, and mass fractions'''
    fractions = get_fractions(paths=list(jobs_data.index), omkm_path=omkm_path)

    '''Read literature values'''
    try:
        lit_data = pd.read_excel(io=excel_path, sheet_name='lit', skiprows=[1])
    except:
        lit_data = None
        x_lit = None
        y_lit = None
        lit_label = None
        x_scaling = None
        y_scaling = None
        y_upper_interval = None
        y_lower_interval = None
    else:
        x_lit = lit_data[desc_labels[0]]
        lit_label = lit_data['Label']
        if n_desc == 2:
            y_lit = lit_data[desc_labels[1]]

            '''Calculate Scaling Between Descriptors'''
            scaling_model = LinearRegression(fit_intercept=True)
            scaling_model.fit(X=x_lit.to_numpy().reshape(-1, 1),
                              y=y_lit.to_numpy().reshape(-1, 1))
            # Isolate contor map data that lies between literature values
            x_lit_min = np.min(x_lit)
            x_lit_max = np.max(x_lit)
            x_scaling_i = list(np.where(np.logical_and(jobs_data[desc_labels[0]] >= x_lit_min,
                                                    jobs_data[desc_labels[0]] <= x_lit_max))[0])
            x_scaling = np.unique(jobs_data[desc_labels[0]][x_scaling_i])
            x_scaling = np.concatenate([[x_lit_min], x_scaling, [x_lit_max]])
            y_scaling = scaling_model.predict(X=x_scaling.reshape(-1, 1)).flatten()

            '''Calculate prediction intervals'''
            y_lit_model = scaling_model.predict(X=x_lit.to_numpy().reshape(-1, 1)).flatten()
            y_lower_interval, y_upper_interval = get_prediction_interval(y_data=y_lit_model,
                                                                        y_model=y_lit,
                                                                        x_data=x_lit,
                                                                        x_linspace=x_scaling,
                                                                        y_linspace=y_scaling,
                                                                        alpha=0.99)

            '''Label as within design space for Pareto Optimal plots'''
            design_space_mask = []
            for i, (x, y) in enumerate(zip(jobs_data[desc_labels[0]],
                                        jobs_data[desc_labels[1]])):
                result = False
                # Check if x value falls within interpolated region
                if x >= x_lit_min and x <= x_lit_max:
                    # Get index from prediction interval
                    j = np.where(x == x_scaling)[0][0]
                    # Check if y value falls within prediction intervals
                    if y >= y_lower_interval[j] and y <= y_upper_interval[j]:
                        result = True
                design_space_mask.append(result)

    '''Create contour maps for coverages'''
    cov_path = Path('./cov')
    cov_path.mkdir(parents=True, exist_ok=True)
    for species_name, cov_out in fractions['covs_out'].iteritems():
        # Prepare plot arguments
        species_path = Path(species_name).with_suffix('.html')
        species_path = cov_path.joinpath(species_path)
        title = '{} Coverages (ML)'.format(species_name)
        if n_desc == 1:
            plot_1d_volcano(path_out=species_path.as_posix(),
                            x=jobs_data[desc_labels[0]],
                            y=cov_out,
                            title=title,
                            x_label=desc_labels[0],
                            y_label=title,
                            hover_label=jobs_labels,
                            x_lit=x_lit,
                            lit_label=lit_label,
                            ymin_cutoff=0., ymax_cutoff=1.)
        elif n_desc == 2:
            plot_contour(path_out=species_path.as_posix(),
                         x=jobs_data[desc_labels[0]],
                         y=jobs_data[desc_labels[1]],
                         z=cov_out,
                         title=title,
                         x_label=desc_labels[0],
                         y_label=desc_labels[1],
                         hover_label=jobs_labels,
                         x_lit=x_lit,
                         y_lit=y_lit,
                         lit_label=lit_label,
                         zmin_cutoff=0.,
                         zmax_cutoff=1.,
                         x_scaling=x_scaling,
                         y_scaling=y_scaling,
                         y_scaling_upper=y_upper_interval,
                         y_scaling_lower=y_lower_interval)
        else:
            err_msg = 'Currently cannot support {} descriptors'.format(n_desc)
            raise ValueError(err_msg)

    '''Create contour maps for mass fractions'''
    mass_frac_path = Path('./mass_frac')
    mass_frac_path.mkdir(parents=True, exist_ok=True)
    for species_name, species_mass_frac in fractions['mass_out'].iteritems():
        # Prepare plot arguments
        species_path = Path(species_name).with_suffix('.html')
        species_path = mass_frac_path.joinpath(species_path)
        title = '{} Mass Fraction'.format(species_name)
        if n_desc == 1:
            plot_1d_volcano(path_out=species_path.as_posix(),
                            x=jobs_data[desc_labels[0]],
                            y=species_mass_frac,
                            title=title,
                            x_label=desc_labels[0],
                            y_label=title,
                            hover_label=jobs_labels,
                            x_lit=x_lit,
                            lit_label=lit_label,
                            ymin_cutoff=0., ymax_cutoff=1.)
        elif n_desc == 2:
            plot_contour(path_out=species_path.as_posix(),
                         x=jobs_data[desc_labels[0]],
                         y=jobs_data[desc_labels[1]],
                         z=species_mass_frac,
                         title=title,
                         x_label=desc_labels[0],
                         y_label=desc_labels[1],
                         hover_label=jobs_labels,
                         x_lit=x_lit,
                         y_lit=y_lit,
                         lit_label=lit_label,
                         zmin_cutoff=0.,
                         zmax_cutoff=1.,
                         x_scaling=x_scaling,
                         y_scaling=y_scaling,
                         y_scaling_upper=y_upper_interval,
                         y_scaling_lower=y_lower_interval)
        else:
            err_msg = 'Currently cannot support {} descriptors'.format(n_desc)
            raise ValueError(err_msg)

    '''Create contour maps for mol fractions'''
    mol_frac_path = Path('./mol_frac')
    mol_frac_path.mkdir(parents=True, exist_ok=True)
    for species_name, species_mol_frac in fractions['mol_out'].iteritems():
        # Prepare plot arguments
        species_path = Path(species_name).with_suffix('.html')
        species_path = mol_frac_path.joinpath(species_path)
        title = '{} Mole Fraction'.format(species_name)

        if n_desc == 1:
            plot_1d_volcano(path_out=species_path.as_posix(),
                            x=jobs_data[desc_labels[0]],
                            y=species_mol_frac,
                            title=title,
                            x_label=desc_labels[0],
                            y_label=title,
                            hover_label=jobs_labels,
                            x_lit=x_lit,
                            lit_label=lit_label,
                            ymin_cutoff=0., ymax_cutoff=1.)
        elif n_desc == 2:
            plot_contour(path_out=species_path.as_posix(),
                         x=jobs_data[desc_labels[0]],
                         y=jobs_data[desc_labels[1]],
                         z=species_mol_frac,
                         title=title,
                         x_label=desc_labels[0],
                         y_label=desc_labels[1],
                         hover_label=jobs_labels,
                         x_lit=x_lit,
                         y_lit=y_lit,
                         lit_label=lit_label,
                         zmin_cutoff=0.,
                         zmax_cutoff=1.,
                         x_scaling=x_scaling,
                         y_scaling=y_scaling,
                         y_scaling_upper=y_upper_interval,
                         y_scaling_lower=y_lower_interval)
        else:
            err_msg = 'Currently cannot support {} descriptors'.format(n_desc)
            raise ValueError(err_msg)

    '''Create Quantity objects from reactor inputs'''
    reactor = get_reactor_inputs(reactor_inputs=inputs['reactor'],
                                 phases=inputs['phases'],
                                 units=inputs['units'])
    '''Calculate flow rates'''
    mol_flow_rates, mass_flow_rates = get_flow_rates(
        reactor=reactor,
        mol_frac_in=fractions['mol_in'],
        mol_frac_out=fractions['mol_out'],
        mass_frac_in=fractions['mass_in'],
        mass_frac_out=fractions['mass_out'],
        mol_weight=mol_weight)
    mol_gen = mol_flow_rates['out'] - mol_flow_rates['in']
    mol_cons = - mol_gen

    '''Calculate Quantities of interest'''
    conv_df = mol_cons[reactant_name]/mol_flow_rates['in'][reactant_name]
    selec_df = mol_gen[product_name]/mol_cons[reactant_name]/inputs['analysis']['Selectivity Ratio']
    yield_df = conv_df * selec_df
    tof_df = mol_cons[reactant_name]/reactor['mol_catalyst']('mol') # 1/s

    '''Rename quantities of interest'''
    conv_df.rename('Conversion of {}'.format(reactant_name), inplace=True)
    selec_df.rename('Selectivity of {}'.format(product_name), inplace=True)
    yield_df.rename('Yield of {}'.format(product_name), inplace=True)
    tof_df.rename('TOF of {}'.format(reactant_name), inplace=True)

    '''Plot Quantities of Interest'''
    conv_title = 'Conversion of {}'.format(reactant_name)
    selec_title = 'Selectivity of {}'.format(product_name)
    yield_title = 'Yield of {}'.format(product_name)
    tof_title = 'log10(TOF) of {}'.format(reactant_name)
    if n_desc == 1:
        plot_1d_volcano(path_out='./conv.html',
                        x=jobs_data[desc_labels[0]],
                        y=conv_df,
                        title=conv_title,
                        x_label=desc_labels[0],
                        y_label=conv_title,
                        hover_label=jobs_labels,
                        x_lit=x_lit,
                        lit_label=lit_label,
                        ymin_cutoff=-0.1, ymax_cutoff=1.1)

        plot_1d_volcano(path_out='./selectivity.html',
                        x=jobs_data[desc_labels[0]],
                        y=selec_df,
                        title=selec_title,
                        x_label=desc_labels[0],
                        y_label=selec_title,
                        hover_label=jobs_labels,
                        x_lit=x_lit,
                        lit_label=lit_label,
                        ymin_cutoff=-0.1, ymax_cutoff=1.1)

        plot_1d_volcano(path_out='./yield.html',
                        x=jobs_data[desc_labels[0]],
                        y=yield_df,
                        title=yield_title,
                        x_label=desc_labels[0],
                        y_label=yield_title,
                        hover_label=jobs_labels,
                        x_lit=x_lit,
                        lit_label=lit_label,
                        ymin_cutoff=-0.1, ymax_cutoff=0.4)

        plot_1d_volcano(path_out='./tof.html',
                        x=jobs_data[desc_labels[0]],
                        y=np.log10(tof_df),
                        title=tof_title,
                        x_label=desc_labels[0],
                        y_label=tof_title,
                        hover_label=jobs_labels,
                        x_lit=x_lit,
                        lit_label=lit_label)
    elif n_desc == 2:
        plot_contour(path_out='./conv.html',
                        x=jobs_data[desc_labels[0]],
                        y=jobs_data[desc_labels[1]],
                        z=conv_df,
                        title=conv_title,
                        x_label=desc_labels[0],
                        y_label=desc_labels[1],
                        hover_label=jobs_labels,
                        x_lit=x_lit,
                        y_lit=y_lit,
                        lit_label=lit_label,
                        zmin_cutoff=0.,
                        zmax_cutoff=1.,
                        x_scaling=x_scaling,
                        y_scaling=y_scaling,
                        y_scaling_upper=y_upper_interval,
                        y_scaling_lower=y_lower_interval)

        plot_contour(path_out='./selectivity.html',
                        x=jobs_data[desc_labels[0]],
                        y=jobs_data[desc_labels[1]],
                        z=selec_df,
                        title=selec_title,
                        x_label=desc_labels[0],
                        y_label=desc_labels[1],
                        hover_label=jobs_labels,
                        x_lit=x_lit,
                        y_lit=y_lit,
                        lit_label=lit_label,
                        zmin_cutoff=0.,
                        zmax_cutoff=1.,
                        x_scaling=x_scaling,
                        y_scaling=y_scaling,
                        y_scaling_upper=y_upper_interval,
                        y_scaling_lower=y_lower_interval)

        plot_contour(path_out='./yield.html',
                     x=jobs_data[desc_labels[0]],
                     y=jobs_data[desc_labels[1]],
                     z=yield_df,
                     title=yield_title,
                     x_label=desc_labels[0],
                     y_label=desc_labels[1],
                     hover_label=jobs_labels,
                     x_lit=x_lit,
                     y_lit=y_lit,
                     lit_label=lit_label,
                     zmin_cutoff=0.,
                     zmax_cutoff=1.,
                     x_scaling=x_scaling,
                     y_scaling=y_scaling,
                     y_scaling_upper=y_upper_interval,
                     y_scaling_lower=y_lower_interval)

        plot_contour(path_out='./tof.html',
                     x=jobs_data[desc_labels[0]],
                     y=jobs_data[desc_labels[1]],
                     z=np.log10(tof_df),
                     title=tof_title,
                     x_label=desc_labels[0],
                     y_label=desc_labels[1],
                     hover_label=jobs_labels,
                     x_lit=x_lit,
                     y_lit=y_lit,
                     lit_label=lit_label,
                     x_scaling=x_scaling,
                     y_scaling=y_scaling,
                     y_scaling_upper=y_upper_interval,
                     y_scaling_lower=y_lower_interval)

        '''Density Plots'''
        plot_density(path_out='density.html',
                    jobs_data=jobs_data,
                    desc_labels=desc_labels,
                    conv_data=conv_df,
                    selec_data=selec_df,
                    yield_data=yield_df,
                    reactant_name=reactant_name,
                    product_name=product_name,
                    hover_label=jobs_labels,
                    lit_data=lit_data,
                    design_space_mask=design_space_mask)
    else:
        err_msg = 'Currently cannot support {} descriptors'.format(n_desc)
        raise ValueError(err_msg)

    complete_jobs_data = pd.concat([jobs_data, conv_df, selec_df, yield_df, tof_df], axis=1)
    complete_jobs_data.to_excel('../complete_summary.xlsx')