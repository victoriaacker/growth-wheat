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
HZ_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer']
ORGAN_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ'] # exposed organs

#: the name of the organs representing a leaf
LEAF_ORGANS_NAMES = set(['sheath', 'blade'])


def from_dataframes(hz_inputs, organ_inputs):
    """
    Convert inputs/outputs from Pandas dataframe to Growth-Wheat format.

    :Parameters:

        - `hz_inputs` (:class:`pandas.DataFrame`) - Hidden zone inputs dataframe to convert, with one line by Hidden zone.
        - `organ_inputs` (:class:`pandas.DataFrame`) - Exposed organ inputs dataframe to convert, with one line by organ.

    :Returns:
        The inputs in a dictionary.

    :Returns Type:
        :class:`dict` of :class:`dict`

    .. seealso:: see :attr:`simulation.Simulation.inputs` for the structure of Growth-Wheat inputs.

    """
    all_hz_dict = {}
    all_organ_dict = {}
    hz_L_calculation_dict = {}
    hz_inputs_columns = hz_inputs.columns.difference(HZ_TOPOLOGY_COLUMNS)
    organ_inputs_columns = organ_inputs.columns.difference(ORGAN_TOPOLOGY_COLUMNS)
    sheath_inputs_grouped = organ_inputs[organ_inputs.organ == 'sheath'].groupby(ORGAN_TOPOLOGY_COLUMNS)
    sheath_inputs_grouped_all_metamers = organ_inputs[organ_inputs.organ == 'sheath'].groupby(['plant', 'axis'])

    for organ_inputs_id, organ_inputs_group in organ_inputs.groupby(ORGAN_TOPOLOGY_COLUMNS):
        # organ
        organ_inputs_series = organ_inputs_group.loc[organ_inputs_group.first_valid_index()]
        organ_inputs_dict = organ_inputs_series[organ_inputs_columns].to_dict()
        all_organ_dict[organ_inputs_id] = organ_inputs_dict

    hz_inputs_grouped = hz_inputs.groupby(HZ_TOPOLOGY_COLUMNS)
    for hz_inputs_id, hz_inputs_group in hz_inputs_grouped:
        # hz
        hz_inputs_series = hz_inputs_group.loc[hz_inputs_group.first_valid_index()]
        hz_inputs_dict = hz_inputs_series[hz_inputs_columns].to_dict()
        all_hz_dict[hz_inputs_id] = hz_inputs_dict

        # Get lengths required for the calculation of the hz length
        previous_hz_id = tuple(list(hz_inputs_id[:2]) + [hz_inputs_id[-1]-1])
        # previous hz length
        if hz_inputs_grouped.groups.has_key(previous_hz_id):
            previous_hz = hz_inputs_grouped.get_group(previous_hz_id)
            previous_hz_length = previous_hz.loc[previous_hz.first_valid_index(), 'hz_L']
        else:
            previous_hz_length = None
            warnings.warn('No previous hz found for hz {}.'.format(hz_inputs_id))

        # previous sheath length
        previous_sheath_id = tuple(list(hz_inputs_id[:2]) + [hz_inputs_id[-1]-1] + ['sheath'])
        if sheath_inputs_grouped.groups.has_key(previous_sheath_id):
            previous_sheath = sheath_inputs_grouped.get_group(previous_sheath_id)
            previous_sheath_visible_length = previous_sheath.loc[previous_sheath.first_valid_index(), 'visible_length']
            if not previous_hz_length:
                previous_sheath_final_hidden_length = previous_sheath.loc[previous_sheath.first_valid_index(), 'final_hidden_length']

        else:
            previous_sheath_visible_length = 0
            warnings.warn('No previous sheath found for hz {}.'.format(hz_inputs_id))

        hz_L_calculation_dict[hz_inputs_id] = {'previous_hz_length': previous_hz_length,
                                                 'previous_sheath_visible_length': previous_sheath_visible_length,
                                                 'previous_sheath_final_hidden_length': previous_sheath_final_hidden_length} #TODO: ajouter les entrenoeuds

    return {'hz': all_hz_dict, 'organs': all_organ_dict, 'hz_L_calculation': hz_L_calculation_dict}

def to_dataframes(data_dict):
    """
    Convert outputs from Growth-Wheat format to Pandas dataframe.

    :Parameters:

        - `data_dict` (:class:`dict`) - The outputs in Growth-Wheat format.

    :Returns:
        One dataframe for hz outputs and one dataframe for organ outputs.

    :Returns Type:
        :class:`tuple` of :class:`pandas.DataFrame`

    .. seealso:: see :attr:`simulation.Simulation.outputs` for the structure of Growth-Wheat outputs.

    """
    dataframes_dict = {}
    for (current_key, current_topology_columns, current_outputs_names) in (('hz', HZ_TOPOLOGY_COLUMNS, simulation.HZ_OUTPUTS),
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
    return dataframes_dict['hz'], dataframes_dict['organs']


def from_MTG(g):
    """
    Convert a MTG to Growth-Wheat inputs.

    :Parameters:

        - g (:class:`openalea.mtg.mtg.MTG`) - A MTG which contains the inputs
          of Growth-Wheat. These inputs are: :mod:`simulation.HZ_INPUTS` and :mod:`simulation.ORGAN_INPUTS`.

    :Returns:
        The inputs of Growth-Wheat.

    :Returns Type:
        :class:`dict` of :class:`dict`

    .. seealso:: see :attr:`simulation.Simulation.inputs` for the structure of Growth-Wheat inputs.

    """
    all_hz_dict = {}
    all_organ_dict = {}
    hz_L_calculation_dict = {}

    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        for axis_vid in g.components_iter(plant_vid):
            axis_label = g.label(axis_vid)
            axis_id = (plant_index, axis_label)
            for metamer_vid in g.components_iter(axis_vid):
                metamer_index = int(g.index(metamer_vid))
                growthwheat_hz_data_from_mtg_organs_data = {}
                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in LEAF_ORGANS_NAMES: continue

                    organ_inputs_from_mtg = g.get_vertex_property(organ_vid)
                    if organ_label == 'sheath':
                        growthwheat_hz_data_from_mtg_organs_data['leaf_Wlig'] = organ_inputs_from_mtg['diameter']
                    elif organ_label == 'blade':
                        growthwheat_hz_data_from_mtg_organs_data['lamina_Lmax'] = organ_inputs_from_mtg['shape_mature_length']
                        growthwheat_hz_data_from_mtg_organs_data['leaf_Wmax'] = organ_inputs_from_mtg['shape_max_width']

                    organ_id = (plant_index, axis_label, metamer_index, organ_label)
                    organ_inputs_from_mtg = g.get_vertex_property(organ_vid)
                    organ_inputs_dict = {}
                    is_valid_organ = True
                    for organ_input_name in simulation.ORGAN_INPUTS:

                        if organ_input_name in organ_inputs_from_mtg:
                            # use the input from the MTG
                            organ_inputs_dict[organ_input_name] = organ_inputs_from_mtg[organ_input_name]
                        else:
                            is_valid_organ = False
                            break
                        if is_valid_organ:
                            all_organ_dict[organ_id] = organ_inputs_dict

                metamer_properties = g.get_vertex_property(metamer_vid)
                if 'hz' in metamer_properties:
                    previous_metamer_vid = g.parent(metamer_vid)
                    hz_id = (plant_index, axis_label, metamer_index)
                    hz_inputs_from_mtg = metamer_properties['hz']
                    hz_inputs_dict = {}

                    is_valid_hz = True
                    for hz_input_name in simulation.HZ_INPUTS:
                        if hz_input_name in hz_inputs_from_mtg:
                            # use the input from the MTG
                            hz_inputs_dict[hz_input_name] = hz_inputs_from_mtg[hz_input_name]
                        elif hz_input_name in growthwheat_hz_data_from_mtg_organs_data:
                            hz_inputs_dict[hz_input_name] = growthwheat_hz_data_from_mtg_organs_data[hz_input_name]
                        else:
                            # use the input from the dataframe
                            if hz_input_name in hz_inputs_group_series:
                                hz_inputs_dict[hz_input_name] = hz_inputs_group_series[hz_input_name]
                            else:
                                is_valid_hz = False
                                break
                    if is_valid_hz:
                        all_hz_dict[hz_id] = hz_inputs_dict

                    # Get lengths required for the calculation of the hz length
                    if previous_metamer_vid is not None:
                        # previous hz length
                        if g.get_vertex_property(previous_metamer_vid).has_key('hz'):
                            previous_hz_length = g.get_vertex_property(previous_metamer_vid)['hz']['hz_L']
                        else:
                            previous_hz_length = None
                            warnings.warn('No previous hz found for hz {}.'.format(hz_id))
                        # previous sheath length
                        previous_metamer_components = {g.class_name(component_vid): component_vid for component_vid in g.components_at_scale(previous_metamer_vid, scale=4)}
                        if previous_metamer_components.has_key('sheath'):
                            previous_sheath_visible_length = g.get_vertex_property(previous_metamer_components['sheath'])['visible_length']
                            if not previous_hz_length:
                                previous_sheath_final_hidden_length = g.get_vertex_property(previous_metamer_components['sheath'])['final_hidden_length']
                        else:
                            previous_sheath_visible_length = 0
                            warnings.warn('No previous sheath found for hz {}.'.format(hz_id))

                        hz_L_calculation_dict[(plant_index, axis_label, metamer_index)] = {'previous_hz_length': previous_hz_length,
                                                                                            'previous_sheath_visible_length': previous_sheath_visible_length,
                                                                                            'previous_sheath_final_hidden_length': previous_sheath_final_hidden_length} #TODO: ajouter les entrenoeuds
                    else:
                        raise Exception('No previous metamer found for hz {}.'.format(hz_id))

    return {'hz': all_hz_dict, 'organs': all_organ_dict, 'hz_L_calculation': hz_L_calculation_dict}

def update_MTG(g, geometrical_model, inputs=None, outputs=None):
    """
    Update a MTG from Growth-Wheat inputs and outputs.

    :Parameters:
            - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to update from the inputs and outputs of Growth-Wheat.

            - `geometrical_model` (:func:`geometrical_model`) - The model which deals with geometry.
              This model must implement a method `add_metamer(mtg, plant_index, axis_label)` to add
              a metamer to a specific axis of a plant in a MTG.

            - inputs (:class:`dict` of :class:`dict`) - Growth-Wheat inputs.
            These inputs are: :mod:`simulation.HZ_INPUTS` and :mod:`simulation.ORGAN_INPUTS`.

            - outputs (:class:`dict` of :class:`dict`) - Growth-Wheat outputs.
            These outputs are: :mod:`simulation.HZ_OUTPUTS` and :mod:`simulation.ORGAN_OUTPUTS`.

    .. seealso::

        * see :attr:`simulation.Simulation.inputs` for the structure of Growth-Wheat inputs.
        * see :attr:`simulation.Simulation.outputs` for the structure of Growth-Wheat outputs.

    """

    # add the properties if needed
    property_names = g.property_names()
    for growthwheat_data_name in set(simulation.HZ_INPUTS_OUTPUTS + simulation.ORGAN_INPUTS_OUTPUTS):
        if growthwheat_data_name not in property_names:
            g.add_property(growthwheat_data_name)
    if 'hz' not in property_names:
        g.add_property('hz')

    if inputs:
        hzs_data_dict = inputs['hz']
        hzs_data_names = simulation.HZ_INPUTS
        organs_data_dict = inputs['organs']
        organs_data_names = simulation.ORGAN_INPUTS
    elif outputs:
        hzs_data_dict = outputs['hz']
        hzs_data_names = simulation.HZ_OUTPUTS
        organs_data_dict = outputs['organs']
        organs_data_names = simulation.ORGAN_OUTPUTS
    else:
        raise Exception('Both inputs and outputs were found to be None')

    # add new metamer(s)
    axis_to_metamers_mapping = {}
    for metamer_id in sorted(hzs_data_dict.iterkeys()):
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

                hz_id = (plant_index, axis_label, metamer_index)
                mtg_organs_data_from_growthwheat_hz_data = {}
                if hz_id in hzs_data_dict:
                    hz_outputs_dict = hzs_data_dict[hz_id]
                    metamer_properties = g.get_vertex_property(metamer_vid)
                    if 'hz' not in metamer_properties:
                        g.property('hz')[metamer_vid] = {}

                    for hz_data in (hzs_data_dict[hz_id], hzs_data_dict[hz_id]):
                        for hz_data_name, hz_data_value in hz_data.iteritems():
                            if hz_data_name not in ('lamina_Lmax', 'leaf_Wmax'):
                                g.property('hz')[metamer_vid][hz_data_name] = hz_data_value
                            else:
                                mtg_organs_data_from_growthwheat_hz_data[hz_data_name] = hz_data_value

                elif 'hz' in g.get_vertex_property(metamer_vid):
                    # remove the 'hz' property from this metamer
                    del g.property('hz')[metamer_vid]

                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in LEAF_ORGANS_NAMES: continue

                    if len(mtg_organs_data_from_growthwheat_hz_data) != 0:
                        if organ_label == 'blade':
                            g.property('shape_mature_length')[organ_vid] = mtg_organs_data_from_growthwheat_hz_data['lamina_Lmax']
                            g.property('shape_max_width')[organ_vid] = mtg_organs_data_from_growthwheat_hz_data['leaf_Wmax']

                    organ_id = (plant_index, axis_label, metamer_index, organ_label)
                    organ_data_dict = organs_data_dict.get(organ_id, {})
                    if organ_data_dict:
                        for organ_data_name in organs_data_names:
                            g.property(organ_data_name)[organ_vid] = organ_data_dict.get(organ_data_name)