import numpy as np

from vunits import constants as c
from pmutt.statmech.lsr import LSR
from pmutt.empirical.nasa import Nasa
from pmutt.empirical.shomate import Shomate
from pmutt.empirical.references import Reference, References
from pmutt.reaction import Reaction
from pmutt.mixture.cov import PiecewiseCovEffect
from pmutt.omkm.reaction import SurfaceReaction, BEP
from pgradd.GroupAdd.Library import GroupLibrary
import pgradd.ThermoChem


def initialize_references(refs_data):
    if len(refs_data) > 0:
        refs = []
        for ref_data in refs_data:
            ref = Reference(**ref_data)
            refs.append(ref)
        refs = References(references=refs)
    else:
        refs = None
    return refs

def initialize_nasa_species(nasas_data):
    nasa_species_dict = {}
    for ind_nasa_species_data in nasas_data:
        name = ind_nasa_species_data['name']
        nasa_species_dict[name] = Nasa(**ind_nasa_species_data)
    return nasa_species_dict

def initialize_shomate_species(shomates_data):
    shomate_species_dict = {}
    for ind_shomate_species_data in shomates_data:
        name = ind_shomate_species_data['name']
        shomate_species_dict[name] = Shomate(**ind_shomate_species_data)
    return shomate_species_dict

def initialize_ga_species(ga_data, descriptors=None, statmech_species=None,
                          ref_species=None):
    ga_species_dict = {}
    for ind_ga_species_data in ga_data:
        name = ind_ga_species_data['name']
        # Create GA species
        ga_lib = GroupLibrary.Load(ind_ga_species_data['library'])
        ga_groups = ga_lib.GetDescriptors(ind_ga_species_data['smiles'])
        ga_species = ga_lib.Estimate(ga_groups, 'thermochem')

        # Calculate thermo properties
        T = np.linspace(ind_ga_species_data.pop('T_low'),
                        ind_ga_species_data.pop('T_high'))
        CpoR = [ga_species.get_CpoR(T=T_i) for T_i in T]
        T_ref = T[0]
        HoRT = ga_species.get_HoRT(T=T_ref)
        SoR = ga_species.get_SoR(T=T_ref)

        if 'reaction' in ind_ga_species_data:
            ga_reaction = Reaction.from_string(ind_ga_species_data['reaction'],
                                               species={**ga_species_dict,
                                                        **statmech_species})

            # If piecewise, check if descriptor falls within value
            if 'low_val' in ind_ga_species_data:
                ga_reaction_species = ga_reaction.get_species()
                for descriptor in descriptors:
                    if descriptor.name in ga_reaction_species:
                        descriptor_val = descriptor.val
                        break
                # Skip the descriptor value if out of range
                if descriptor_val < ind_ga_species_data['low_val']:
                    continue
                if descriptor_val > ind_ga_species_data['high_val']:
                    continue

            ref_name = name.replace('(S)', '(PT)')
            ga_lsr = LSR(slope=ind_ga_species_data['slope'],
                         intercept=ind_ga_species_data['intercept'],
                         reaction=ga_reaction)
            BE_m = ga_lsr.get_H(T=T_ref, units='kcal/mol')
            BE_ref = ref_species[ref_name].get_H(T=T_ref, units='kcal/mol')
            delta_HoRT = (BE_m - BE_ref)/c.R('kcal/mol/K')/T_ref
            HoRT = HoRT + delta_HoRT

        # Create NASA polynomial from data
        ga_species_dict[name] = Nasa.from_data(T=T,
                                               CpoR=CpoR,
                                               T_ref=T_ref,
                                               HoRT_ref=HoRT,
                                               SoR_ref=SoR,
                                               **ind_ga_species_data)

    return ga_species_dict

def initialize_lsr_species(lsr_data, statmech_species, references=None):
    lsr_species_dict = {}
    for ind_lsr_species_data in lsr_data:
        name = ind_lsr_species_data['name']

        # Create the reaction
        ind_lsr_species_data['reaction'] = Reaction.from_string(
                reaction_str=ind_lsr_species_data['reaction'],
                species=statmech_species)

        # Create the surface species
        try:
            ind_lsr_species_data['surf_species'] = \
                    statmech_species[ind_lsr_species_data['surf_species']]
        except KeyError:
            ind_lsr_species_data['surf_species'] = 0.

        # Create the gas species
        try:
            ind_lsr_species_data['gas_species'] = \
                        statmech_species[ind_lsr_species_data['gas_species']]
        except KeyError:
            ind_lsr_species_data['gas_species'] = 0.

        '''Create a NASA polynomial from the data'''
        lsr_species_dict[name] = Nasa.from_model(references=references,
                                                 **ind_lsr_species_data)
    return lsr_species_dict

def initialize_bep_species(beps_data):
    beps_dict = {bep_data['name']: BEP(**bep_data) for bep_data in beps_data}
    return beps_dict

def initialize_interactions(interactions_data):
    interactions = [PiecewiseCovEffect(**interaction_data) \
                    for interaction_data in interactions_data]
    return interactions

def initialize_reactions(reactions_data, species):
    reactions = [SurfaceReaction.from_string(species=species, **reaction_data)\
                 for reaction_data in reactions_data]
    return reactions