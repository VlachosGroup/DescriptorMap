import json
import os
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

if __name__ == "__main__":
    '''Initialization'''
    try:
        os.chdir(os.path.dirname(__file__))
    except FileNotFoundError:
        pass

    '''Read inputs'''
    in_path = Path(r'../inputs.json')
    excel_path = Path(r'../inputs.xlsx')
    omkm_path = Path(r'../omkm/')
    with open(in_path, 'r') as f_ptr:
        inputs = json.load(f_ptr)

    '''Units'''
    units = Units(**inputs['units'])

    '''Read original data'''
    orig_data = {}
    sheet_names = ('refs', 'dft_species', 'extended_lsr_species', 'ga_species',
                   'nasa_species', 'shomate_species', 'beps',
                   'reactions', 'interactions', 'phases')
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
            summary_vals = ['E_{}'.format(name) for name in descriptor_names]
            summary_header.extend(summary_indices  + summary_vals)

        '''Create directory'''
        if sampling == 'linear':
            job_path_list = []
            for descriptor in descriptors:
                folder_format = '{:0%dd}' % descriptor.field_len
                folder_name = folder_format.format(descriptor.i)
                job_path_list.append(folder_name)
            job_path = Path(*job_path_list)
        elif sampling == 'lhs':
            folder_format = '{:0%dd}' % descriptors[0].field_len
            job_path = Path(folder_format.format(descriptors[0].i))
        rel_job_path = omkm_path.joinpath(job_path)
        rel_job_path.mkdir(parents=True, exist_ok=True)
        folder_list.append(job_path.as_posix())
        print('Processing {}'.format(job_path))

        '''Create paths for OpenMKM files'''
        yaml_path = rel_job_path.joinpath('reactor.yaml')
        cti_path = rel_job_path.joinpath('thermo.cti')

        '''Record descriptor data for later saving'''
        summary_indices = [descriptor.i for descriptor in descriptors]
        summary_vals = [descriptor.val for descriptor in descriptors]
        summary_data.append([job_path.as_posix()] \
                             + summary_indices + summary_vals)

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
                statmech_species=statmech_species_dict)
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
        all_species_dict = {**non_bep_species_dict,
                            **beps_dict}

        '''Reactions'''
        reactions = setup.initialize_reactions(deepcopy(orig_data['reactions']),
                                               all_species_dict)

        '''Phases'''
        phases = organize_phases(phases_data=deepcopy(orig_data['phases']),
                                 species=list(all_species_dict.values()),
                                 reactions=reactions,
                                 interactions=interactions)

        '''Write YAML File'''
        write_yaml(units=units, phases=phases, filename=yaml_path,
                   **inputs['reactor'])

        '''Write CTI File'''
        write_cti(phases=phases, species=non_bep_species_dict.values(),
                  reactions=reactions, lateral_interactions=interactions,
                  units=units, filename=cti_path, T=inputs['reactor']['T'],
                  P=inputs['reactor']['P'])

        if write_json:
            '''Input composition data for analysis'''
            inputs['composition'] = {}
            for name, ind_species in all_species_dict.items():
                elements = ind_species.elements
                # Skip species that do not have elements
                if elements is None:
                    continue
                inputs['composition'][name] = elements

            '''Update JSON file'''
            with open(in_path, 'w', newline='\n') as f_ptr:
                json.dump(inputs, f_ptr, indent=2)

            write_json = False

    '''Write folderlist'''
    folder_list_path = omkm_path.joinpath('folderlist.txt')
    with open(folder_list_path, 'w', newline='\n') as f_ptr:
        f_ptr.write('\n'.join(folder_list))

    '''Write summary data'''
    summary_path = in_path.parent.joinpath('job_summary.xlsx')
    summary_data = pd.DataFrame(summary_data, columns=summary_header,
                                index=folder_list)
    summary_data.set_index('Path')
    summary_data.to_excel(summary_path)
    # inputs['paths']['job_summary'] = summary_path.as_posix()