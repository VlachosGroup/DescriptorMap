import numpy as np

from vunits import constants as c
from pmutt.statmech.lsr import LSR, ExtendedLSR
from pmutt.empirical.nasa import Nasa
from pmutt.empirical.shomate import Shomate
from pmutt.empirical.references import Reference, References
from pmutt.reaction import Reaction
from pmutt.reaction import ChemkinReaction, Reactions
from pmutt.mixture.cov import PiecewiseCovEffect
from pmutt.omkm.reaction import SurfaceReaction, BEP
from pgradd.GroupAdd.Library import GroupLibrary
import pgradd.ThermoChem


def initialize_references(refs_data):
    """Initialie reference species from spreadsheet data

    Parameters
    ----------
        refs_data : list of dicts
            References data. Each element of the list has the keyword arguments
            to create a `pMuTT Reference`_ object.
    Returns
    -------
        refs : `pMuTT References`_ object
            Initialized references.

    .. _`pMuTT Reference`: https://vlachosgroup.github.io/pMuTT/api/empirical/references/pmutt.empirical.references.Reference.html
    .. _`pMuTT References`: https://vlachosgroup.github.io/pMuTT/api/empirical/references/pmutt.empirical.references.References.html
    """
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
    """Initialize NASA polynomials from spreadsheet data

    Parameters
    ----------
        nasas_data : list of dicts
            NASAs data. Each element of the list corresponds to the keyword
            arguments that can initialize a `pMuTT Nasa`_ object.
    Returns
    -------
        nasas_dict : dict of `pMuTT Nasa`_ objects
            Initialized `pMuTT Nasa`_ objects where the ``name`` attribute is
            the key.
            
    .. _`pMuTT Nasa`: https://vlachosgroup.github.io/pMuTT/api/empirical/nasa/pmutt.empirical.nasa.Nasa.html
    """
    nasa_species_dict = {}
    for ind_nasa_species_data in nasas_data:
        name = ind_nasa_species_data['name']
        nasa_species_dict[name] = Nasa(**ind_nasa_species_data)
    return nasa_species_dict

def initialize_shomate_species(shomates_data):
    """Initialize Shomate polynomials from spreadsheet data

    Parameters
    ----------
        shomates_data : list of dicts
            Shomates data. Each element of the list corresponds to the keyword
            arguments that can initialize a `pMuTT Shomate`_ object.
    Returns
    -------
        shomates_dict : dict of `pMuTT Shomate`_ objects
            Shomate species where the ``name`` attribute is the key.
            
    .. _`pMuTT Shomate`: https://vlachosgroup.github.io/pMuTT/api/empirical/shomate/pmutt.empirical.shomate.Shomate.html
    """
    shomate_species_dict = {}
    for ind_shomate_species_data in shomates_data:
        name = ind_shomate_species_data['name']
        shomate_species_dict[name] = Shomate(**ind_shomate_species_data)
    return shomate_species_dict

def initialize_ga_species(ga_data, descriptors=None, statmech_species=None,
                          ref_species=None):
    """Initialize group additivity species from spreadsheet data

    Parameters
    ----------
        ga_data : list of dicts
            Group additivity data. Each element of the list corresponds to the
            keyword arguments that can initialize a ``pGrAdd Group`` object.
        descriptors : list of str, optional
            Names of descriptors to use to adjust enthalpy. Default is None.
        statmech_species : dict of `pMuTT StatMech`_ objects, optional
            StatMech objects used to initialize LSRs. Default is None.
        ref_species : dict of `pMuTT StatMech`_ objects, optional
            Reference species to adjust the enthalpy. If the pGrAdd library
            is over a particular surface, the reference species should
            represent those species over that surface. Default is None.
    Returns
    -------
        ga_species_dict : dict of `pMuTT Nasa`_ objects
            Group additivity species converted to `pMuTT Nasa` objects where
            the ``name`` attribute is the key.
            
    .. _`pMuTT Nasa`: https://vlachosgroup.github.io/pMuTT/api/empirical/nasa/pmutt.empirical.nasa.Nasa.html
    """
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
            ga_lsr = LSR(slope=ind_ga_species_data['slope'],
                         intercept=ind_ga_species_data['intercept'],
                         reaction=ga_reaction)
        elif 'reactions' in ind_ga_species_data:
            ga_reactions = []
            for reaction_str in ind_ga_species_data['reactions']:
                reaction = Reaction.from_string(reaction_str,
                                                species={**ga_species_dict,
                                                         **statmech_species})
                ga_reactions.append(reaction)

            ga_lsr = ExtendedLSR(slopes=ind_ga_species_data['slopes'],
                                 intercept=ind_ga_species_data['intercept'],
                                 reactions=ga_reactions)

        # Add enthalpy correction if LSR is used
        try:
            ga_lsr
        except NameError:
            pass
        else:
            # If piecewise, check if descriptor falls within value
            if 'low_val' in ind_ga_species_data \
               or 'high_val' in ind_ga_species_data:
                ga_reaction_species = ga_reaction.get_species()
                # Find the desciptor and its corresponding value
                for descriptor in descriptors:
                    if descriptor.name in ga_reaction_species:
                        descriptor_val = descriptor.val
                        break
                # Skips applying adjustment to enthalpy if descriptor value
                # out of range
                if descriptor_val < ind_ga_species_data['low_val']:
                    continue
                if descriptor_val > ind_ga_species_data['high_val']:
                    continue

            # Adjust binding energy by taking different between reference
            # species and new species
            # TODO Make ref_name more general
            ref_name = name.replace('(S)', '(PT)')
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
    """Initialize LSR species from spreadsheet data

    Parameters
    ----------
        lsr_data : list of dicts
            Linear scaling relationship data. Each element of the list
            corresponds to the keyword arguments that can initialize a
            `pMuTT LSR`_ object.
        statmech_species : dict of `pMuTT StatMech`_ objects
            StatMech objects used to initialize LSRs. Default is None.
        references : `pMuTT References`_ object, optional
            References to use to adjust enthalpy when creating `pMuTT Nasa`_
            objects. Default is None
    Returns
    -------
        lsr_species_dict : dict of `pMuTT Nasa`_ objects
            Linear scaling relationships converted to `pMuTT Nasa_` objects
            where the ``name`` attribute is the key

    .. _`pMuTT LSR`: https://vlachosgroup.github.io/pMuTT/api/statmech/elec/pmutt.statmech.lsr.LSR.html
    .. _`pMuTT Nasa`: https://vlachosgroup.github.io/pMuTT/api/empirical/nasa/pmutt.empirical.nasa.Nasa.html
    .. _`pMuTT References`: https://vlachosgroup.github.io/pMuTT/api/empirical/references/pmutt.empirical.references.References.html
    """
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

def initialize_extended_lsr_species(extended_lsr_data,
                                    statmech_species,
                                    references=None):
    """Initialize extended LSR species from spreadsheet data

    Parameters
    ----------
        extended_lsr_data : list of dicts
            Extended linear scaling relationship data. Each element of the list
            corresponds to the keyword arguments that can initialize a
            `pMuTT ExtendedLSR`_ object.
        statmech_species : dict of `pMuTT StatMech`_ objects
            StatMech objects used to initialize Extended LSRs. Default is None.
        references : `pMuTT References`_ object, optional
            References to use to adjust enthalpy when creating `pMuTT Nasa`_
            objects. Default is None
    Returns
    -------
        lsr_species_dict : dict of `pMuTT Nasa`_ objects
            Linear scaling relationships converted to `pMuTT Nasa_` objects
            where the ``name`` attribute is the key.

    .. _`pMuTT ExtendedLSR`: https://vlachosgroup.github.io/pMuTT/api/statmech/elec/pmutt.statmech.lsr.ExtendedLSR.html
    .. _`pMuTT Nasa`: https://vlachosgroup.github.io/pMuTT/api/empirical/nasa/pmutt.empirical.nasa.Nasa.html
    .. _`pMuTT References`: https://vlachosgroup.github.io/pMuTT/api/empirical/references/pmutt.empirical.references.References.html
    """
    lsr_species_dict = {}
    for ind_lsr_species_data in extended_lsr_data:
        name = ind_lsr_species_data['name']

        reactions = []
        for reaction_str in ind_lsr_species_data['reactions']:
            reaction = Reaction.from_string(reaction_str=reaction_str,
                                            species=statmech_species)
            reactions.append(reaction)
        ind_lsr_species_data['reactions'] = reactions

        surf_species = []
        if 'surf_species' in ind_lsr_species_data:
            for ind_surf_species in ind_lsr_species_data['surf_species']:
                # Create the surface species
                try:
                    new_species = statmech_species[ind_surf_species]
                except KeyError:
                    new_species = 0.
                surf_species.append(new_species)
            ind_lsr_species_data['surf_species'] = surf_species

        gas_species = []
        if 'gas_species' in ind_lsr_species_data:
            for ind_gas_species in ind_lsr_species_data['gas_species']:
                # Create the gasace species
                try:
                    new_species = statmech_species[ind_gas_species]
                except KeyError:
                    new_species = 0.
                gas_species.append(new_species)
            ind_lsr_species_data['gas_species'] = gas_species

        '''Create a NASA polynomial from the data'''
        lsr_species_dict[name] = Nasa.from_model(references=references,
                                                 **ind_lsr_species_data)
    return lsr_species_dict
    

def initialize_bep_species(beps_data):
    """Initialize BEP species from spreadsheet data

    Parameters
    ----------
        beps_data : list of dict
            BEP relationship data. Each element of the list correspond to the
            keyword arguments that can initialize a `pMuTT BEP`_ object.
    Returns
    -------
        beps_dict : dict of `pMuTT BEP`_ objects
            Initialized BEP relationships where the ``name`` attribute is the
            key.

    .. _`pMuTT BEP`: https://vlachosgroup.github.io/pMuTT/api/reactions/bep/pmutt.reaction.bep.BEP.html
    """
    beps_dict = {bep_data['name']: BEP(**bep_data) for bep_data in beps_data}
    return beps_dict

def initialize_interactions(interactions_data):
    """Initialize lateral interactions from spreadsheet data

    Parameters
    ----------
        interactions_data : list of dict
            Lateral interactions data. Each element of the list correspond to
            the keyword arguments that can initialize a
            `pMuTT PiecewiseCovEffect`_ object.
    Returns
    -------
        interactions : list of `pMuTT PiecewiseCovEffect`_ objects
            Initialized lateral interactions.

    .. _`pMuTT PiecewiseCovEffect`: https://vlachosgroup.github.io/pMuTT/api/mixture/cov/pmutt.mixture.cov.PiecewiseCovEffect.html
    """
    interactions = [PiecewiseCovEffect(**interaction_data) \
                    for interaction_data in interactions_data]
    return interactions

def initialize_reactions(reactions_data, species):
    """Initialize reactions

    Parameters
    ----------
        reactions_data : list of dict
            Reactions data. Each element of the list correspond to
            the keyword arguments that can initialize a
            `pMuTT SurfaceReaction`_ object using the ``from_string`` method.
        species : dict of pMuTT Model objects
            Species corresponding to the reaction strings in ``reaction_data``
    Returns
    -------
        reactions : list of `pMuTT SurfaceReaction`_ objects
            Surface reactions that can be written to OpenMKM format

    .. _`pMuTT SurfaceReaction`: https://vlachosgroup.github.io/pMuTT/api/reactions/reactions/pmutt.omkm.reaction.SurfaceReaction.html#pmutt.omkm.reaction.SurfaceReaction
    """
    reactions = [SurfaceReaction.from_string(species=species, **reaction_data) \
                 for reaction_data in reactions_data]
    return reactions

def initialize_chemkin_reactions(reactions_data, species):
    """Initialize reactions

    Parameters
    ----------
        reactions_data : list of dict
            Reactions data. Each element of the list correspond to
            the keyword arguments that can initialize a
            `pMuTT ChemkinReaction`_ object using the ``from_string`` method.
        species : dict of pMuTT Model objects
            Species corresponding to the reaction strings in ``reaction_data``
    Returns
    -------
        reactions : list of `pMuTT ChemkinReaction`_ objects
            Surface reactions that can be written to Chemkin format

    .. _`pMuTT ChemkinReaction`: https://vlachosgroup.github.io/pMuTT/api/reactions/reactions/pmutt.reaction.ChemkinReaction.html#pmutt.reaction.ChemkinReaction
    """
    reactions = Reactions([ChemkinReaction.from_string(species=species, **reaction_data) \
                           for reaction_data in reactions_data])
    return reactions