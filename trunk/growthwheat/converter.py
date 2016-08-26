# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division
import warnings

"""
    growthwheat.converter
    ~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`growthwheat.converter` defines functions to:

        * convert a :class:`dataframe <pandas.DataFrame>` to/from GrowthWheat inputs or outputs format.
        * convert a :class:`MTG <openalea.mtg.mtg.MTG>` to/from GrowthWheat inputs or outputs format.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2015.
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""

import pandas as pd

from openalea.mtg import io, fat_mtg

import simulation

#: the columns which define the topology in the input/output dataframe
HGZ_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer']
ORGAN_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ'] # exposed organs

#: the name of the organs representing a leaf
LEAF_ORGANS_NAMES = set(['sheath', 'blade'])


def from_dataframes(hgz_inputs, organ_inputs):
    """
    Convert inputs/outputs from Pandas dataframe to Growth-Wheat format.

    :Parameters:

        - `hgz_inputs` (:class:`pandas.DataFrame`) - Hidden growing zone inputs dataframe to convert, with one line by Hidden growing zone.
        - `organ_inputs` (:class:`pandas.DataFrame`) - Exposed organ inputs dataframe to convert, with one line by organ.

    :Returns:
        The inputs in a dictionary.

    :Returns Type:
        :class:`dict` of :class:`dict`

    .. seealso:: see :attr:`simulation.Simulation.inputs` for the structure of Growth-Wheat inputs.

    """
    all_hgz_dict = {}
    all_organ_dict = {}
    hgz_L_calculation_dict = {}
    hgz_inputs_columns = hgz_inputs.columns.difference(HGZ_TOPOLOGY_COLUMNS)
    organ_inputs_columns = organ_inputs.columns.difference(ORGAN_TOPOLOGY_COLUMNS)
    sheath_inputs_grouped = organ_inputs[organ_inputs.organ == 'sheath'].groupby(ORGAN_TOPOLOGY_COLUMNS)
    sheath_inputs_grouped_all_metamers = organ_inputs[organ_inputs.organ == 'sheath'].groupby(['plant', 'axis'])

    for organ_inputs_id, organ_inputs_group in organ_inputs.groupby(ORGAN_TOPOLOGY_COLUMNS):
        # organ
        organ_inputs_series = organ_inputs_group.loc[organ_inputs_group.first_valid_index()]
        organ_inputs_dict = organ_inputs_series[organ_inputs_columns].to_dict()
        all_organ_dict[organ_inputs_id] = organ_inputs_dict

    hgz_inputs_grouped = hgz_inputs.groupby(HGZ_TOPOLOGY_COLUMNS)
    for hgz_inputs_id, hgz_inputs_group in hgz_inputs_grouped:
        # hgz
        hgz_inputs_series = hgz_inputs_group.loc[hgz_inputs_group.first_valid_index()]
        hgz_inputs_dict = hgz_inputs_series[hgz_inputs_columns].to_dict()
        all_hgz_dict[hgz_inputs_id] = hgz_inputs_dict

        # Get lengths required for the calculation of the hgz length
        previous_hgz_id = tuple(list(hgz_inputs_id[:2]) + [hgz_inputs_id[-1]-1])
        # previous hgz length
        if hgz_inputs_grouped.groups.has_key(previous_hgz_id):
            previous_hgz = hgz_inputs_grouped.get_group(previous_hgz_id)
            previous_hgz_length = previous_hgz.loc[previous_hgz.first_valid_index(), 'hgz_L']
        else:
            previous_hgz_length = 0
            warnings.warn('No previous hgz found for hgz {}.'.format(hgz_inputs_id))

        # previous sheath length
        previous_sheath_id = tuple(list(hgz_inputs_id[:2]) + [hgz_inputs_id[-1]-1] + ['sheath'])
        if sheath_inputs_grouped.groups.has_key(previous_sheath_id):
            previous_sheath = sheath_inputs_grouped.get_group(previous_sheath_id)
            previous_sheath_length = previous_sheath.loc[previous_sheath.first_valid_index(), 'length']
        else:
            previous_sheath_length = 0
            warnings.warn('No previous sheath found for hgz {}.'.format(hgz_inputs_id))

        hgz_L_calculation_dict[hgz_inputs_id] = (previous_hgz_length, previous_sheath_length) #TODO: ajouter les entrenoeuds + mettre en dict

    return {'hgz': all_hgz_dict, 'organs': all_organ_dict, 'hgz_L_calculation': hgz_L_calculation_dict}

def to_dataframes(data_dict):
    """
    Convert outputs from Growth-Wheat format to Pandas dataframe.

    :Parameters:

        - `data_dict` (:class:`dict`) - The outputs in Growth-Wheat format.

    :Returns:
        One dataframe for hgz outputs and one dataframe for organ outputs.

    :Returns Type:
        :class:`tuple` of :class:`pandas.DataFrame`

    .. seealso:: see :attr:`simulation.Simulation.outputs` for the structure of Growth-Wheat outputs.

    """
    dataframes_dict = {}
    for (current_key, current_topology_columns, current_outputs_names) in (('hgz', HGZ_TOPOLOGY_COLUMNS, simulation.HGZ_OUTPUTS),
                                                                           ('organs', ORGAN_TOPOLOGY_COLUMNS, simulation.ORGAN_OUTPUTS)):
        current_data_dict = data_dict[current_key]
        current_ids_df = pd.DataFrame(current_data_dict.keys(), columns=current_topology_columns)
        current_data_df = pd.DataFrame(current_data_dict.values())
        current_df = pd.concat([current_ids_df, current_data_df], axis=1)
        current_df.sort_values(by=current_topology_columns, inplace=True)
        current_columns_sorted = current_topology_columns + current_outputs_names
        current_df = current_df.reindex_axis(current_columns_sorted, axis=1, copy=False)
        current_df.reset_index(drop=True, inplace=True)
        dataframes_dict[current_key] = current_df
    return dataframes_dict['hgz'], dataframes_dict['organs']


def from_MTG(g):
    """
    Convert a MTG to Growth-Wheat inputs.

    :Parameters:

        - g (:class:`openalea.mtg.mtg.MTG`) - A MTG which contains the inputs
          of Growth-Wheat. These inputs are: :mod:`simulation.HGZ_INPUTS` and :mod:`simulation.ORGAN_INPUTS`.

    :Returns:
        The inputs of Growth-Wheat.

    :Returns Type:
        :class:`dict` of :class:`dict`

    .. seealso:: see :attr:`simulation.Simulation.inputs` for the structure of Growth-Wheat inputs.

    """
    all_hgz_dict = {}
    all_organ_dict = {}
    hgz_L_calculation_dict = {}

    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        for axis_vid in g.components_iter(plant_vid):
            axis_label = g.label(axis_vid)
            axis_id = (plant_index, axis_label)
            for metamer_vid in g.components_iter(axis_vid):
                metamer_index = int(g.index(metamer_vid))
                growthwheat_hgz_data_from_mtg_organs_data = {}
                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in LEAF_ORGANS_NAMES: continue

                    organ_inputs_from_mtg = g.get_vertex_property(organ_vid)

                    if organ_label == 'sheath':
                        growthwheat_hgz_data_from_mtg_organs_data['leaf_Wlig'] = organ_inputs_from_mtg['diameter']
                    elif organ_label == 'blade':
                        growthwheat_hgz_data_from_mtg_organs_data['lamina_Lmax'] = organ_inputs_from_mtg['shape_mature_length']
                        growthwheat_hgz_data_from_mtg_organs_data['leaf_Wmax'] = organ_inputs_from_mtg['shape_max_width']

                    organ_id = (plant_index, axis_label, metamer_index, organ_label)
                    organ_inputs_from_mtg = g.get_vertex_property(organ_vid)
                    organ_inputs_dict = {}
                    is_valid_organ = True
                    for organ_input_name in simulation.ORGAN_INPUTS:

                        if organ_input_name in organ_inputs_from_mtg:
                            # use the input from the MTG
                            organ_inputs_dict[organ_input_name] = organ_inputs_from_mtg[organ_input_name]
                        else:
                            # use the input from the dataframe
                            if organ_input_name in organ_inputs_group_series:
                                organ_inputs_dict[organ_input_name] = organ_inputs_group_series[organ_input_name]
                            else:
                                is_valid_organ = False
                                break
                        if is_valid_organ:
                            all_organ_dict[organ_id] = organ_inputs_dict

                metamer_properties = g.get_vertex_property(metamer_vid)
                if 'hgz' in metamer_properties:
                    previous_metamer_vid = g.parent(metamer_vid)
                    hgz_id = (plant_index, axis_label, metamer_index)
                    hgz_inputs_from_mtg = metamer_properties['hgz']
                    hgz_inputs_dict = {}

                    is_valid_hgz = True
                    for hgz_input_name in simulation.HGZ_INPUTS:
                        if hgz_input_name in hgz_inputs_from_mtg:
                            # use the input from the MTG
                            hgz_inputs_dict[hgz_input_name] = hgz_inputs_from_mtg[hgz_input_name]
                        elif hgz_input_name in growthwheat_hgz_data_from_mtg_organs_data:
                            hgz_inputs_dict[hgz_input_name] = growthwheat_hgz_data_from_mtg_organs_data[hgz_input_name]
                        else:
                            # use the input from the dataframe
                            if hgz_input_name in hgz_inputs_group_series:
                                hgz_inputs_dict[hgz_input_name] = hgz_inputs_group_series[hgz_input_name]
                            else:
                                is_valid_hgz = False
                                break
                    if is_valid_hgz:
                        all_hgz_dict[hgz_id] = hgz_inputs_dict

                    # Get lengths required for the calculation of the hgz length
                    if previous_metamer_vid is not None:
                        # previous hgz length
                        if g.get_vertex_property(previous_metamer_vid).has_key('hgz'):
                            previous_hgz_length = g.get_vertex_property(previous_metamer_vid)['hgz']['hgz_L']
                        else:
                            previous_hgz_length = 0
                            warnings.warn('No previous hgz found for hgz {}.'.format(hgz_id))
                        # previous sheath length
                        previous_metamer_components = {g.class_name(component_vid): component_vid for component_vid in g.components_at_scale(previous_metamer_vid, scale=4)}
                        if previous_metamer_components.has_key('sheath'):
                            previous_sheath_length = g.get_vertex_property(previous_metamer_components['sheath'])['length']
                        else:
                            previous_sheath_length = 0
                            warnings.warn('No previous sheath found for hgz {}.'.format(hgz_id))

                        hgz_L_calculation_dict[(plant_index, axis_label, metamer_index)] = (previous_hgz_length, previous_sheath_length) #TODO: ajouter les entrenoeuds + dico
                    else:
                        raise Exception('No previous metamer found for hgz {}.'.format(hgz_id))

    return {'hgz': all_hgz_dict, 'organs': all_organ_dict, 'hgz_L_calculation': hgz_L_calculation_dict}

def update_MTG(g, geometrical_model, inputs=None, outputs=None):
    """
    Update a MTG from Growth-Wheat inputs and outputs.

    :Parameters:
            - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to update from the inputs and outputs of Growth-Wheat.

            - `geometrical_model` (:func:`geometrical_model`) - The model which deals with geometry.
              This model must implement a method `add_metamer(mtg, plant_index, axis_label)` to add
              a metamer to a specific axis of a plant in a MTG.

            - inputs (:class:`dict` of :class:`dict`) - Growth-Wheat inputs.
            These inputs are: :mod:`simulation.HGZ_INPUTS` and :mod:`simulation.ORGAN_INPUTS`.

            - outputs (:class:`dict` of :class:`dict`) - Growth-Wheat outputs.
            These outputs are: :mod:`simulation.HGZ_OUTPUTS` and :mod:`simulation.ORGAN_OUTPUTS`.

    .. seealso::

        * see :attr:`simulation.Simulation.inputs` for the structure of Growth-Wheat inputs.
        * see :attr:`simulation.Simulation.outputs` for the structure of Growth-Wheat outputs.

    """

    # add the properties if needed
    property_names = g.property_names()
    for growthwheat_data_name in set(simulation.HGZ_INPUTS_OUTPUTS + simulation.ORGAN_INPUTS_OUTPUTS):
        if growthwheat_data_name not in property_names:
            g.add_property(growthwheat_data_name)
    if 'hgz' not in property_names:
        g.add_property('hgz')

    if inputs:
        hgzs_data_dict = inputs['hgz']
        organs_data_dict = inputs['organs']
    elif outputs:
        hgzs_data_dict = outputs['hgz']
        organs_data_dict = outputs['organs']
    else:
        raise Exception('Both inputs and outputs were found to be None')

    # add new metamer(s)
    axis_to_metamers_mapping = {}
    for metamer_id in sorted(hgzs_data_dict.iterkeys()):
        axis_id = (metamer_id[0], metamer_id[1])
        if axis_id not in axis_to_metamers_mapping:
            axis_to_metamers_mapping[axis_id] = []
        axis_to_metamers_mapping[axis_id].append(metamer_id)

    axis_to_old_metamers_mapping = {}
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        for axis_vid in g.components_iter(plant_vid):
            axis_label = g.label(axis_vid)
            metamer_ids = set([(plant_index, axis_label, int(g.index(metamer_vid))) for metamer_vid in g.components_iter(axis_vid)])
            if (plant_index, axis_label) not in axis_to_metamers_mapping: continue
            new_metamer_ids = set(axis_to_metamers_mapping[(plant_index, axis_label)]).difference(metamer_ids)
            for new_metamer_id in new_metamer_ids:
                geometrical_model.add_metamer(g, plant_index, axis_label)

    # update the properties of the MTG
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        for axis_vid in g.components_iter(plant_vid):
            axis_label = g.label(axis_vid)
            for metamer_vid in g.components_iter(axis_vid):
                metamer_index = int(g.index(metamer_vid))

                hgz_id = (plant_index, axis_label, metamer_index)
                mtg_organs_data_from_growthwheat_hgz_data = {}
                if hgz_id in hgzs_data_dict:
                    hgz_outputs_dict = hgzs_data_dict[hgz_id]
                    metamer_properties = g.get_vertex_property(metamer_vid)
                    if 'hgz' not in metamer_properties:
                        g.property('hgz')[metamer_vid] = {}

                    for hgz_data in (hgzs_data_dict[hgz_id], hgzs_data_dict[hgz_id]):
                        for hgz_data_name, hgz_data_value in hgz_data.iteritems():
                            if hgz_data_name not in ('lamina_Lmax', 'leaf_Wmax'):
                                g.property('hgz')[metamer_vid][hgz_data_name] = hgz_data_value
                            else:
                                mtg_organs_data_from_growthwheat_hgz_data[hgz_data_name] = hgz_data_value

                elif 'hgz' in g.get_vertex_property(metamer_vid):
                    # remove the 'hgz' property from this metamer
                    del g.property('hgz')[metamer_vid]

                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in LEAF_ORGANS_NAMES: continue

                    if len(mtg_organs_data_from_growthwheat_hgz_data) != 0:
                        if organ_label == 'blade':
                            g.property('shape_mature_length')[organ_vid] = mtg_organs_data_from_growthwheat_hgz_data['lamina_Lmax']
                            g.property('shape_max_width')[organ_vid] = mtg_organs_data_from_growthwheat_hgz_data['leaf_Wmax']

                    organ_id = (plant_index, axis_label, metamer_index, organ_label)
                    organ_inputs_dict = organs_data_dict.get(organ_id, {})
                    organ_outputs_dict = organs_data_dict.get(organ_id, {})

                    organ_properties = g.get_vertex_property(organ_vid)
                    organ_properties.update(organ_inputs_dict)
                    organ_properties.update(organ_outputs_dict)

                    for organ_input_name in simulation.ORGAN_INPUTS:
                        g.property(organ_input_name)[organ_vid] = organ_inputs_dict.get(organ_input_name)
                    for organ_output_name in simulation.ORGAN_OUTPUTS:
                        g.property(organ_output_name)[organ_vid] = organ_outputs_dict.get(organ_output_name)