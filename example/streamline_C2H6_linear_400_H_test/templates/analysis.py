import os
import json
import yaml
from pathlib import Path

import numpy as np
import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go

from pmutt.io.excel import read_excel
from pmutt import get_molecular_weight
from vunits import constants as c
from vunits.quantity import Quantity
from descmap.analysis import (get_reactor_inputs, plot_contour, get_fractions,
    get_flow_rates, plot_density, plot_1d_volcano)

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
    else:
        
        x_lit = lit_data[desc_label[0]]
        # y_lit = lit_data[desc_label[1]]
        lit_label = lit_data['Label']


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
                            x=jobs_data[desc_label[0]],
                            y=cov_out,
                            title=title,
                            x_label=x_label,
                            y_label=title,
                            hover_label=jobs_labels,
                            x_lit=x_lit,
                            lit_label=lit_label,
                            ymin_cutoff=0., ymax_cutoff=1.)
        elif n_desc == 2:
            plot_contour(path_out=species_path.as_posix(),
                         x=jobs_data[desc_label[0]],
                         y=jobs_data[desc_labels[1]],
                         z=cov_out,
                         title=title,
                         x_label=x_label,
                         y_label=title,
                         hover_label=jobs_labels,
                         x_lit=x_lit,
                         y_lit=y_lit,
                         lit_label=lit_label,
                         zmin_cutoff=0.,
                         zmax_cutoff=1.)
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
                            x=jobs_data[desc_label[0]],
                            y=species_mass_frac,
                            title=title,
                            x_label=x_label,
                            y_label=title,
                            hover_label=jobs_labels,
                            x_lit=x_lit,
                            lit_label=lit_label,
                            ymin_cutoff=0., ymax_cutoff=1.)
        elif n_desc == 2:
            plot_contour(path_out=species_path.as_posix(),
                         x=jobs_data[desc_label[0]],
                         y=jobs_data[desc_labels[1]],
                         z=species_mass_frac,
                         title=title,
                         x_label=x_label,
                         y_label=title,
                         hover_label=jobs_labels,
                         x_lit=x_lit,
                         y_lit=y_lit,
                         lit_label=lit_label,
                         zmin_cutoff=0.,
                         zmax_cutoff=1.)
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
                            x=jobs_data[desc_label[0]],
                            y=species_mol_frac,
                            title=title,
                            x_label=x_label,
                            y_label=title,
                            hover_label=jobs_labels,
                            x_lit=x_lit,
                            lit_label=lit_label,
                            ymin_cutoff=0., ymax_cutoff=1.)
        elif n_desc == 2:
            plot_contour(path_out=species_path.as_posix(),
                         x=jobs_data[desc_label[0]],
                         y=jobs_data[desc_labels[1]],
                         z=species_mol_frac,
                         title=title,
                         x_label=x_label,
                         y_label=title,
                         hover_label=jobs_labels,
                         x_lit=x_lit,
                         y_lit=y_lit,
                         lit_label=lit_label,
                         zmin_cutoff=0.,
                         zmax_cutoff=1.)
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
    # Conversion
    title = 'Conversion of {}'.format(reactant_name)
    if n_desc == 1:
        plot_1d_volcano(path_out='./conv.html',
                        x=jobs_data[desc_label[0]],
                        y=conv_df,
                        title=title,
                        x_label=x_label,
                        y_label=title,
                        hover_label=jobs_labels,
                        x_lit=x_lit,
                        lit_label=lit_label,
                        ymin_cutoff=0., ymax_cutoff=1.)
    elif n_desc == 2:
        plot_contour(path_out='./conv.html',
                        x=jobs_data[desc_label[0]],
                        y=jobs_data[desc_labels[1]],
                        z=conv_df,
                        title=title,
                        x_label=x_label,
                        y_label=title,
                        hover_label=jobs_labels,
                        x_lit=x_lit,
                        y_lit=y_lit,
                        lit_label=lit_label,
                        zmin_cutoff=0.,
                        zmax_cutoff=1.)
    else:
        err_msg = 'Currently cannot support {} descriptors'.format(n_desc)
        raise ValueError(err_msg)

    # Selectivity
    title = 'Selectivity of {}'.format(product_name)
    if n_desc == 1:
        plot_1d_volcano(path_out='./selectivity.html',
                        x=jobs_data[desc_label[0]],
                        y=selec_df,
                        title=title,
                        x_label=x_label,
                        y_label=title,
                        hover_label=jobs_labels,
                        x_lit=x_lit,
                        lit_label=lit_label,
                        ymin_cutoff=0., ymax_cutoff=1.)
    elif n_desc == 2:
        plot_contour(path_out='./selectivity.html',
                        x=jobs_data[desc_label[0]],
                        y=jobs_data[desc_labels[1]],
                        z=selec_df,
                        title=title,
                        x_label=x_label,
                        y_label=title,
                        hover_label=jobs_labels,
                        x_lit=x_lit,
                        y_lit=y_lit,
                        lit_label=lit_label,
                        zmin_cutoff=0.,
                        zmax_cutoff=1.)

    # Yield
    title = 'Yield of {}'.format(product_name)
    if n_desc == 1:
        plot_1d_volcano(path_out='./yield.html',
                        x=jobs_data[desc_label[0]],
                        y=yield_df,
                        title=title,
                        x_label=x_label,
                        y_label=title,
                        hover_label=jobs_labels,
                        x_lit=x_lit,
                        lit_label=lit_label,
                        ymin_cutoff=0., ymax_cutoff=1.)
    elif n_desc == 2:
        plot_contour(path_out='./yield.html',
                        x=jobs_data[desc_label[0]],
                        y=jobs_data[desc_labels[1]],
                        z=yield_df,
                        title=title,
                        x_label=x_label,
                        y_label=title,
                        hover_label=jobs_labels,
                        x_lit=x_lit,
                        y_lit=y_lit,
                        lit_label=lit_label,
                        zmin_cutoff=0.,
                        zmax_cutoff=1.)

    # TOF
    title = 'log10(TOF) of {}'.format(product_name)
    if n_desc == 1:
        plot_1d_volcano(path_out='./tof.html',
                        x=jobs_data[desc_label[0]],
                        y=np.log10(tof_df),
                        title=title,
                        x_label=x_label,
                        y_label=title,
                        hover_label=jobs_labels,
                        x_lit=x_lit,
                        lit_label=lit_label)
    elif n_desc == 2:
        plot_contour(path_out='./tof.html',
                        x=jobs_data[desc_label[0]],
                        y=jobs_data[desc_labels[1]],
                        z=np.log10(tof_df),
                        title=title,
                        x_label=x_label,
                        y_label=title,
                        hover_label=jobs_labels,
                        x_lit=x_lit,
                        y_lit=y_lit,
                        lit_label=lit_label)

    '''Density Plots'''
    plot_density(path_out='density.html',
                 jobs_data=jobs_data,
                 desc_labels=desc_label,
                 conv_data=conv_df,
                 selec_data=selec_df,
                 yield_data=yield_df,
                 reactant_name=reactant_name,
                 product_name=product_name,
                 hover_label=jobs_labels,
                 lit_data=lit_data)

    complete_jobs_data = pd.concat([jobs_data, conv_df, selec_df, yield_df, tof_df], axis=1)
    complete_jobs_data.to_excel('../complete_summary.xlsx')