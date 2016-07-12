# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    growthwheat.converter
    ~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`growthwheat.converter` defines functions to:

        * convert a :class:`dataframe <pandas.DataFrame>` to/from GrowthWheat inputs or outputs format.
        * convert a :class:`MTG <openalea.mtg.mtg.MTG>` to/from GrowthWheat inputs or outputs format.

    Both dataframes and MTG follow AdelWheat naming convention.

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
ELEMENT_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ'] # exposed elements


def from_dataframes(hgz_inputs_df, element_inputs_df):
    """
    Convert inputs/outputs from Pandas dataframe to Growth-Wheat format.
    The column names of the dataframe respect the naming convention of AdelWheat.

    :Parameters:

        - `hgz_inputs_df` (:class:`pandas.DataFrame`) - Hidden growing zone inputs dataframe to convert, with one line by Hidden growing zone.
        - `element_inputs_df` (:class:`pandas.DataFrame`) - Exposed element inputs dataframe to convert, with one line by element.

    :Returns:
        The inputs in a dictionary.

    :Returns Type:
        :class:`dict` of :class:`dict`

    .. seealso:: see :attr:`simulation.Simulation.inputs` for the structure of Growth-Wheat inputs.

    """
    all_hgz_dict = {}
    all_element_dict = {}
    hgz_to_prev_sheath_dict = {}
    hgz_inputs_columns = hgz_inputs_df.columns.difference(HGZ_TOPOLOGY_COLUMNS)
    element_inputs_columns = element_inputs_df.columns.difference(ELEMENT_TOPOLOGY_COLUMNS)
    sheath_inputs_grouped = element_inputs_df[element_inputs_df.organ == 'sheath'].groupby(ELEMENT_TOPOLOGY_COLUMNS)
    sheath_inputs_grouped_all_metamers = element_inputs_df[element_inputs_df.organ == 'sheath'].groupby(['plant', 'axis'])

    for element_inputs_id, element_inputs_group in element_inputs_df.groupby(ELEMENT_TOPOLOGY_COLUMNS):
        # element
        element_inputs_series = element_inputs_group.loc[element_inputs_group.first_valid_index()]
        element_inputs_dict = element_inputs_series[element_inputs_columns].to_dict()
        all_element_dict[element_inputs_id] = element_inputs_dict

    for hgz_inputs_id, hgz_inputs_group in hgz_inputs_df.groupby(HGZ_TOPOLOGY_COLUMNS):
        # hgz
        hgz_inputs_series = hgz_inputs_group.loc[hgz_inputs_group.first_valid_index()]
        hgz_inputs_dict = hgz_inputs_series[hgz_inputs_columns].to_dict()
        all_hgz_dict[hgz_inputs_id] = hgz_inputs_dict

        # Length of the previous sheaths
        index_metamer = hgz_inputs_id[-1]
        sheath_inputs_group_all_metamers = sheath_inputs_grouped_all_metamers.get_group(tuple(hgz_inputs_id[:2])) #: all sheaths of the hgz
        for mid in reversed(range(1, index_metamer)):
            previous_sheath_id = tuple(list(hgz_inputs_id[:2]) + [mid, 'sheath'])
            if sheath_inputs_grouped.groups.has_key(previous_sheath_id):
                previous_sheath = sheath_inputs_grouped.get_group(previous_sheath_id)
                if previous_sheath['is_growing'].iloc[0]: #: the previous sheath is growing
                    previous_growing_sheaths = sheath_inputs_group_all_metamers.loc[sheath_inputs_group_all_metamers.is_growing] #: all previous growing sheaths
                    if len(previous_growing_sheaths) > 1:
                        raise Warning('Several previous growing sheaths found.') #: In wheat there is usually a single growing sheath emerged
                    previous_growing_sheaths_L = previous_growing_sheaths.length.sum() #: Sum of the previous growing sheath lengths (although only 1 previous growing sheath is expected, see Warning above)
                    previous_mature_sheath_index = previous_growing_sheaths.index[0] - 1 #: Found the first previous mature sheath
                    if previous_mature_sheath_index in sheath_inputs_group_all_metamers.index:
                        previous_mature_sheath_L = sheath_inputs_group_all_metamers.loc[previous_mature_sheath_index, 'length']
                    else:
                        previous_mature_sheath_L = 0
                        raise Warning('No previous mature sheath found.')
                else:
                    previous_growing_sheaths_L = 0
                    previous_mature_sheath_L = previous_sheath['length'].iloc[0]
                break
        else:
            raise Exception('Error: no previous sheath found.')

        hgz_to_prev_sheath_dict[hgz_inputs_id] = (previous_growing_sheaths_L, previous_mature_sheath_L) #TODO: distinguer emerge de total dans noms variables

    return {'hgz': all_hgz_dict, 'elements': all_element_dict, 'previous_sheaths_L': hgz_to_prev_sheath_dict}

def to_dataframes(data_dict):
    """
    TODO

    """
    dataframes_dict = {}
    for (current_key, current_topology_columns, current_outputs_names) in (('hgz', HGZ_TOPOLOGY_COLUMNS, simulation.HGZ_OUTPUTS),
                                                                           ('elements', ELEMENT_TOPOLOGY_COLUMNS, simulation.ELEMENT_OUTPUTS)):
        current_data_dict = data_dict[current_key]
        current_ids_df = pd.DataFrame(current_data_dict.keys(), columns=current_topology_columns)
        current_data_df = pd.DataFrame(current_data_dict.values())
        current_df = pd.concat([current_ids_df, current_data_df], axis=1)
        current_df.sort_index(by=current_topology_columns, inplace=True)
        current_columns_sorted = current_topology_columns + current_outputs_names
        current_df = current_df.reindex_axis(current_columns_sorted, axis=1, copy=False)
        current_df.reset_index(drop=True, inplace=True)
        dataframes_dict[current_key] = current_df
    return dataframes_dict['hgz'], dataframes_dict['elements']


def from_MTG(g, inputs):
    """
    Convert a MTG to Growth-Wheat inputs.
    Use data in `inputs` if `g` is incomplete.
    The property names in the MTG respect the naming convention of AdelWheat.

    :Parameters:

            - g (:class:`openalea.mtg.mtg.MTG`) - A MTG which contains the inputs
              of Growth-Wheat. These inputs are: :mod:`simulation.GROWTHWHEAT_INPUTS`.

            - `inputs` (:class:`pandas.DataFrame`) - Inputs dataframe, with one line by element.

    :Returns:
        The inputs of Growth-Wheat.

    :Returns Type:
        :class:`dict`

    .. seealso:: see :attr:`simulation.Simulation.inputs` for the structure of Growth-Wheat inputs.

    """
    all_inputs_dict = {}

    inputs_grouped = inputs.groupby(HGZ_TOPOLOGY_COLUMNS)

    # traverse the MTG recursively from top ...
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        for axis_vid in g.components_iter(plant_vid):
            axis_label = g.label(axis_vid)
            for axis_component_vid in g.components_iter(axis_vid):
                if not g.label(axis_component_vid).startswith('metamer'): continue
                metamer_vid = axis_component_vid
                metamer_index = int(g.index(metamer_vid))
                # create a new hidden growing zone
                metamer_has_hgz = False
                hgz_id = (plant_index, axis_label, metamer_index)
                if hgz_id in inputs_grouped.groups:
                    hgz_inputs_group = inputs_grouped.get_group(hgz_id)
                    hgz_inputs_group_series = hgz_inputs_group.loc[hgz_inputs_group.first_valid_index(), hgz_inputs_group.columns.intersection(simulation.HGZ_INPUTS)]
                else:
                    hgz_inputs_group_series = pd.Series()
                hgz_inputs = {}
                is_valid_hgz = True
                for organ_vid in g.components_iter(metamer_vid):
                    if not g.label(organ_vid).startswith('HiddenGrowingZone'): continue
                    hgz_vid = organ_vid
                    vertex_properties = g.get_vertex_property(hgz_vid)
                    for hgz_input_name in simulation.HGZ_INPUTS:
                        if hgz_input_name in vertex_properties:
                            # use the properties of the vertex
                            hgz_inputs[hgz_input_name] = vertex_properties[hgz_input_name]
                        else:
                            # use the value in hgz_inputs_group_series
                            if hgz_input_name in hgz_inputs_group_series:
                                hgz_inputs[hgz_input_name] = hgz_inputs_group_series[hgz_input_name]
                            else:
                                is_valid_hgz = False
                                break
                    if is_valid_hgz:
                        all_inputs_dict[hgz_id] = hgz_inputs
                        metamer_has_hgz = True
                    break
                if not metamer_has_hgz:
                    if set(hgz_inputs_group_series.index).issuperset(simulation.HGZ_INPUTS):
                        hgz_inputs.update(hgz_inputs_group_series.to_dict())
                        all_inputs_dict[hgz_id] = hgz_inputs
    return all_inputs_dict


def update_MTG(outputs, g):
    """
    Update a MTG from Growth-Wheat inputs and outputs.
    The property names in the MTG respect the naming convention of AdelWheat.

    :Parameters:
            - outputs (:class:`dict` of :class:`dict`) - Growth-Wheat outputs.
            These outputs are: :mod:`simulation.GROWTHWHEAT_OUTPUTS`.

            - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to update from the inputs and outputs of Growth-Wheat.

    .. seealso:: see :attr:`simulation.Simulation.inputs` and :attr:`simulation.Simulation.outputs` for the structure of Growth-Wheat inputs and outputs.

    """

    # add the properties if needed
    property_names = g.property_names()
    for growthwheat_output_name in simulation.HGZ_INPUTS_OUTPUTS:
        if growthwheat_output_name not in property_names:
            g.add_property(growthwheat_output_name)

    # traverse the MTG recursively from top ...
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        for axis_vid in g.components_iter(plant_vid):
            axis_label = g.label(axis_vid)
            for axis_component_vid in g.components_iter(axis_vid):
                if not g.label(axis_component_vid).startswith('metamer'): continue
                metamer_vid = axis_component_vid
                metamer_index = int(g.index(metamer_vid))

                hgz_id = (plant_index, axis_label, metamer_index)
                if hgz_id not in outputs: continue

                metamer_has_hgz = False
                metamer_components_iter = g.components_iter(metamer_vid)
                first_metamer_component_vid = next(metamer_components_iter)
                first_metamer_component_label = g.label(first_metamer_component_vid)
                if first_metamer_component_label == 'HiddenGrowingZone':
                    hgz_vid = first_metamer_component_vid
                else:
                    hgz_vid = insert_parent_at_all_scales(g, first_metamer_component_vid, label='HiddenGrowingZone')[0]
                    g = fat_mtg(g)
                for organ_property_name in simulation.HGZ_INPUTS_OUTPUTS:
                    g.property(organ_property_name)[hgz_vid] = outputs[hgz_id][organ_property_name]


##########################################################################################################
####### TODO: move the following to openalea.mtg #########################################################
##########################################################################################################

def insert_parent_at_all_scales(g, parent_id, edge_type='+', label='roots'):
    added_vertices = []
    scale = g.scale(parent_id)
    max_scale = g.max_scale()
    edge_types = g.properties()['edge_type']

    # Add a parent
    if g.parent(parent_id) is None:
        vid = g.insert_parent(parent_id, label=label, edge_type='/')
        edge_types[parent_id] = edge_type
    else:
        return added_vertices

    # Update complex and components
    cid = g.complex(parent_id)
    croots = g._components[cid]
    if parent_id in croots:
        i = croots.index(parent_id)
        g._components[cid][i] = vid
        g._complex[vid] = cid
        print 'ADDED', vid
    else:
        print 'ERROR'

    added_vertices.append(vid)

    while scale+1 <= max_scale:
        pid = g.component_roots(parent_id)[0]
        component_id = g.add_component(vid)
        vid = g.insert_parent(pid, parent_id=component_id, label=label, edge_type='/')
        edge_types[pid] = edge_type
        scale += 1
        added_vertices.append(vid)
        parent_id = pid
        print 'ADDED', vid

    return added_vertices