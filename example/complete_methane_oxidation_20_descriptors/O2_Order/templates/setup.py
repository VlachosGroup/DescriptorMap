import json
import os
import shutil
from copy import deepcopy
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from xlrd.biffh import XLRDError

from descmap import setup
from descmap.sampling import sampling_map
from descmap.errors import raise_invalid_sampling_method
from pmutt.statmech import StatMech
from pmutt.empirical.nasa import Nasa
from pmutt.io.excel import read_excel
from pmutt.io.omkm import organize_phases, write_cti, write_yaml
from pmutt.omkm.units import Units

from pmutt.chemkin import CatSite
from pmutt.io.thermdat import write_thermdat
from pmutt.io import chemkin as ck_io


if __name__ == "__main__":
    '''Initialization'''
    try:
        os.chdir(os.path.dirname(__file__))
    except FileNotFoundError:
        pass

    '''Read inputs'''
    in_path = Path(r'../inputs.json')
    excel_path = Path(r'../inputs.xlsx')
    mkm_path = Path(r'../omkm/')
    with open(in_path, 'r') as f_ptr:
        inputs = json.load(f_ptr)

    '''Units'''
    units = Units(**inputs['units'])

    '''Read original data'''
    orig_data = {}
    sheet_names = ('refs', 'dft_species', 'extended_lsr_species', 'ga_species',
                   'nasa_species', 'shomate_species', 'beps',
                   'reactions', 'interactions', 'cat_sites')
    for sheet_name in sheet_names:
        try:
            orig_data[sheet_name] = read_excel(io=excel_path,
                                               sheet_name=sheet_name)
        except XLRDError:
            # Assign an empty list if the sheet does not exist
            orig_data[sheet_name] = []

    '''References'''
    refs = setup.initialize_references(orig_data['refs'])

    '''Descriptor Sampling'''
    descriptor_names = [descriptor_data['name'] \
                        for descriptor_data in inputs['descriptors']]
    sampling = inputs['fields']['sampling'].lower()
    sampling_fn = sampling_map[sampling]
    descriptors_map = sampling_fn(inputs['descriptors'])

    '''Create iterator to control jobs'''
    if sampling == 'linear':
        descriptors_iter = product(*descriptors_map)
    elif sampling == 'lhs':
        descriptors_iter = zip(*descriptors_map)
    else:
        raise_invalid_sampling_method(sampling)

    folder_list = []
    summary_data = []
    summary_header = ['Path']
    write_json = True
    for descriptors in descriptors_iter:
        '''Summary header'''
        # Only need to make headers once
        if len(summary_header) == 1:
            summary_indices = ['i_{}'.format(name) for name in descriptor_names]
            summary_vals = ['{}'.format(name) for name in descriptor_names]
            summary_header.extend(summary_indices  + summary_vals)

        '''Create directory'''
        if sampling == 'linear':
            job_path_list = []
            for descriptor in descriptors:
                #folder_format = '{:0%dd}' % descriptor.field_len
                folder_name = '{}'.format(descriptor.i)
                job_path_list.append(folder_name)
            job_path = Path(*job_path_list)
        elif sampling == 'lhs':
            #folder_format = '{:0%dd}' % descriptors[0].field_len
            job_path = Path('{}'.format(descriptors[0].i))
        rel_job_path = mkm_path.joinpath(job_path)
        rel_job_path.mkdir(parents=True, exist_ok=True)
        folder_list.append(job_path.as_posix())
        print('Processing {}'.format(job_path))

        '''Create paths for Chemkin files'''
        inp_dir = rel_job_path.joinpath('INP.d')
        graph_path = rel_job_path.joinpath('graphviz.sh')
        surf_path = inp_dir.joinpath('surf.inp')
        gas_path = inp_dir.joinpath('gas.inp')
        EAs_path = inp_dir.joinpath('EAs.inp')
        T_flow_path = inp_dir.joinpath('T_flow.inp')
        tube_mole_path = inp_dir.joinpath('tube_mole.inp')
        thermdat_path = inp_dir.joinpath('thermdat')
        tube_cov_path = inp_dir.joinpath('tube_COV.inp')
        tube_path = inp_dir.joinpath('tube.inp')

        '''Record descriptor data for later saving'''
        summary_indices = [descriptor.i for descriptor in descriptors]
        summary_vals = [descriptor.val for descriptor in descriptors]
        summary_data.append([job_path.as_posix()] \
                             + summary_indices + summary_vals)

        '''Cat Sites'''
        cat_site_data = read_excel(io=excel_path, sheet_name='cat_sites')[0]
        cat_site = CatSite(**cat_site_data)

        '''Species'''
        dft_species_dict = {}
        dft_species_data = deepcopy(orig_data['dft_species'])
        '''Adjust descriptor data'''
        for ind_species_data in dft_species_data:
            '''If the species is a changing descriptor, change the value'''
            ind_name = ind_species_data['name']
            # Skip species if they are not descriptors
            if ind_name not in descriptor_names:
                continue
            # Find descriptor index
            i = descriptor_names.index(ind_name)
            ind_species_data['potentialenergy'] = descriptors[i].val

        '''Create NASA polynomials from species data'''
        statmech_species_dict = {
            ind_species_data['name']: StatMech(**ind_species_data) \
                for ind_species_data in dft_species_data}

        dft_species_dict = {
                ind_species_data['name']: Nasa.from_model(references=refs,
                                                          **ind_species_data) \
                for ind_species_data in dft_species_data}

        '''Initilize all required species objects'''
        nasa_species_dict = setup.initialize_nasa_species(
                nasas_data=deepcopy(orig_data['nasa_species']))
        shomate_species_dict = setup.initialize_shomate_species(
                shomates_data=deepcopy(orig_data['shomate_species']))
        lsr_species_dict = setup.initialize_extended_lsr_species(
                extended_lsr_data=deepcopy(orig_data['extended_lsr_species']),
                statmech_species=statmech_species_dict,
                references=refs)
        ga_species_dict = setup.initialize_ga_species(
                ga_data=deepcopy(orig_data['ga_species']),
                statmech_species=statmech_species_dict,
                descriptors=descriptors,
                ref_species={**lsr_species_dict, **statmech_species_dict})
        beps_dict = setup.initialize_bep_species(
                beps_data=deepcopy(orig_data['beps']))
        interactions = setup.initialize_interactions(
                interactions_data=deepcopy(orig_data['interactions']))

        '''Combine all species'''
        non_bep_species_dict = {**dft_species_dict,
                                **lsr_species_dict,
                                **nasa_species_dict,
                                **shomate_species_dict,
                                **ga_species_dict}
        non_bep_species_list = [x for x in non_bep_species_dict.values()]
        all_species_dict = {**non_bep_species_dict,
                            **beps_dict}

        '''Adjust catalytic sites for Chemkin'''
        for species_name in non_bep_species_dict.keys():
            non_bep_species_dict[species_name].cat_site = cat_site

        '''Reactions'''
        #TO DO: Check if returns Chemkin reactions
        reactions = setup.initialize_chemkin_reactions(deepcopy(orig_data['reactions']),
                                                       all_species_dict)

#        '''Phases'''
#        phases = organize_phases(phases_data=deepcopy(orig_data['phases']),
#                                 species=list(all_species_dict.values()),
#                                 reactions=reactions,
#                                 interactions=interactions)

        '''Create INP.d'''
        if not os.path.exists(inp_dir):
            os.mkdir(inp_dir)
        
        '''Write Chemkin files'''
        write_thermdat(filename=thermdat_path, nasa_species=non_bep_species_list)

        # ### Writing gas.inp and surf.inp
        ck_io.write_gas(filename=gas_path,
                        nasa_species=non_bep_species_list,
                        reactions=reactions,
                        act_method_name='get_G_act',
                        act_unit='kcal/mol')
        ck_io.write_surf(filename=surf_path,
                         reactions=reactions,
                         act_method_name='get_G_act',
                         ads_act_method='get_H_act',
                         act_unit='kcal/mol',
                         T=613.,
                         P=1)

        # ### Writing T_flow.inp
        # Conditions used to write files
        ck_io.write_T_flow(filename=T_flow_path, T=[613.,613.,613.,613.,613.,],
                           P=[1.,1.,1.,1.,1.], Q=[1.25,1.25,1.25,1.25,1.25], 
                           abyv=[3588,3588,3588,3588,3588])

        # ### Writing EAg.inp and EAs.inp
        # Convert T_flow inputs into list of dictionaries that can be used by write_EA.
        # In the future, this will be replaced by a function
        conditions = [{'T': 613, 'P': 1., 'Q': 1.25},
                      {'T': 613, 'P': 1., 'Q': 1.25},
                      {'T': 613, 'P': 1., 'Q': 1.25},
                      {'T': 613, 'P': 1., 'Q': 1.25},
                      {'T': 613, 'P': 1., 'Q': 1.25}]
                      
        ck_io.write_EA(filename=EAs_path,
                       reactions=reactions,
                       write_gas_phase=False,
                       act_method_name='get_GoRT_act',
                       ads_act_method='get_HoRT_act',
                       conditions=conditions)
        
        # ### Writing tube_mole.inp
        # Generating a list of conditions to input
        mole_frac_conditions = [{'CH4': 0.01, 'O2': 0.002, 'HE': 0.988, 'PT(S)': 1.00},
                                {'CH4': 0.01, 'O2': 0.004, 'HE': 0.986, 'PT(S)': 1.00},
                                {'CH4': 0.01, 'O2': 0.006, 'HE': 0.984, 'PT(S)': 1.00},
                                {'CH4': 0.01, 'O2': 0.008, 'HE': 0.982, 'PT(S)': 1.00},
                                {'CH4': 0.01, 'O2': 0.010, 'HE': 0.980, 'PT(S)': 1.00}]
        
        # Write the tube_mole.inp file
        # TO DO Check filename
        ck_io.write_tube_mole(mole_frac_conditions=mole_frac_conditions, 
                              nasa_species=non_bep_species_list, 
                              filename=tube_mole_path)

        '''Copy bep.inp and tube_cov.inp'''
        shutil.copyfile(r'../templates/tube_COV.inp', tube_cov_path)
        shutil.copyfile(r'../templates/tube.inp', tube_path)
        shutil.copyfile(r'../templates/graphviz.sh', graph_path)
#        os.chmod(chemkin_path, 0o775)

    '''Write folderlist'''
    folder_list_path = mkm_path.joinpath('folderlist.txt')
    with open(folder_list_path, 'w', newline='\n') as f_ptr:
        f_ptr.write('\n'.join(folder_list))

    '''Write summary data'''
    summary_path = in_path.parent.joinpath('job_summary.xlsx')
    summary_data = pd.DataFrame(summary_data, columns=summary_header,
                                index=folder_list)
    summary_data.set_index('Path')
    summary_data.to_excel(summary_path)
    # inputs['paths']['job_summary'] = summary_path.as_posix()
