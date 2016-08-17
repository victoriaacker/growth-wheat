# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

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
ELEMENT_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ'] # exposed elements

#: the name of the organs representing a leaf
LEAF_ORGANS_NAMES = set(['sheath', 'blade'])


def from_dataframes(hgz_inputs, element_inputs):
    """
    Convert inputs/outputs from Pandas dataframe to Growth-Wheat format.

    :Parameters:

        - `hgz_inputs` (:class:`pandas.DataFrame`) - Hidden growing zone inputs dataframe to convert, with one line by Hidden growing zone.
        - `element_inputs` (:class:`pandas.DataFrame`) - Exposed element inputs dataframe to convert, with one line by element.

    :Returns:
        The inputs in a dictionary.

    :Returns Type:
        :class:`dict` of :class:`dict`

    .. seealso:: see :attr:`simulation.Simulation.inputs` for the structure of Growth-Wheat inputs.

    """
    all_hgz_dict = {}
    all_element_dict = {}
    hgz_to_prev_sheath_dict = {}
    hgz_inputs_columns = hgz_inputs.columns.difference(HGZ_TOPOLOGY_COLUMNS)
    element_inputs_columns = element_inputs.columns.difference(ELEMENT_TOPOLOGY_COLUMNS)
    sheath_inputs_grouped = element_inputs[element_inputs.organ == 'sheath'].groupby(ELEMENT_TOPOLOGY_COLUMNS)
    sheath_inputs_grouped_all_metamers = element_inputs[element_inputs.organ == 'sheath'].groupby(['plant', 'axis'])

    for element_inputs_id, element_inputs_group in element_inputs.groupby(ELEMENT_TOPOLOGY_COLUMNS):
        # element
        element_inputs_series = element_inputs_group.loc[element_inputs_group.first_valid_index()]
        element_inputs_dict = element_inputs_series[element_inputs_columns].to_dict()
        all_element_dict[element_inputs_id] = element_inputs_dict

    for hgz_inputs_id, hgz_inputs_group in hgz_inputs.groupby(HGZ_TOPOLOGY_COLUMNS):
        # hgz
        hgz_inputs_series = hgz_inputs_group.loc[hgz_inputs_group.first_valid_index()]
        hgz_inputs_dict = hgz_inputs_series[hgz_inputs_columns].to_dict()
        all_hgz_dict[hgz_inputs_id] = hgz_inputs_dict

        # Length of the previous sheaths
        index_metamer = hgz_inputs_id[-1]
        axis_id = hgz_inputs_id[:2]
        sheath_inputs_group_all_metamers = sheath_inputs_grouped_all_metamers.get_group(tuple(axis_id)) #: all sheaths of the axis (i.e. of all metamers of the axis)
        for mid in reversed(range(1, index_metamer)):
            previous_sheath_id = tuple(list(axis_id) + [mid, 'sheath'])
            if sheath_inputs_grouped.groups.has_key(previous_sheath_id):
                previous_sheath = sheath_inputs_grouped.get_group(previous_sheath_id)
                if previous_sheath.loc[previous_sheath.first_valid_index(), 'is_growing']: #: the previous sheath is growing
                    previous_growing_sheaths = sheath_inputs_group_all_metamers.loc[(sheath_inputs_group_all_metamers.metamer <= mid) & sheath_inputs_group_all_metamers.is_growing] #: all previous growing sheaths
                    if len(previous_growing_sheaths) > 1:
                        raise Warning('Several previous growing sheaths found.') #: In wheat there is usually a single growing sheath emerged
                    previous_growing_sheaths_L = previous_growing_sheaths.length.sum() #: Sum of the previous growing sheath lengths (although only 1 previous growing sheath is expected, see Warning above)
                    previous_mature_sheath_index = previous_growing_sheaths.first_valid_index() - 1 #: Found the first previous mature sheath
                    if previous_mature_sheath_index in sheath_inputs_group_all_metamers.index:
                        previous_mature_sheath_L = sheath_inputs_group_all_metamers.loc[previous_mature_sheath_index, 'length']
                    else:
                        previous_mature_sheath_L = 0
                        raise Warning('No previous mature sheath found.')
                else:
                    previous_growing_sheaths_L = 0
                    previous_mature_sheath_L = previous_sheath.loc[previous_sheath.first_valid_index(), 'length']
                break
        else:
            raise Exception('Error: no previous sheath found.')

        hgz_to_prev_sheath_dict[hgz_inputs_id] = (previous_growing_sheaths_L, previous_mature_sheath_L) #TODO: distinguer emerge de total dans noms variables

    return {'hgz': all_hgz_dict, 'elements': all_element_dict, 'previous_sheaths_L': hgz_to_prev_sheath_dict}


def to_dataframes(data_dict):
    """
    Convert outputs from Growth-Wheat format to Pandas dataframe.

    :Parameters:

        - `data_dict` (:class:`dict`) - The outputs in Growth-Wheat format.

    :Returns:
        One dataframe for hgz outputs and one dataframe for element outputs.

    :Returns Type:
        :class:`tuple` of :class:`pandas.DataFrame`

    .. seealso:: see :attr:`simulation.Simulation.outputs` for the structure of Growth-Wheat outputs.

    """
    dataframes_dict = {}
    for (current_key, current_topology_columns, current_outputs_names) in (('hgz', HGZ_TOPOLOGY_COLUMNS, simulation.HGZ_OUTPUTS),
                                                                           ('elements', ELEMENT_TOPOLOGY_COLUMNS, simulation.ELEMENT_OUTPUTS)):
        current_data_dict = data_dict[current_key]
        current_ids_df = pd.DataFrame(current_data_dict.keys(), columns=current_topology_columns)
        current_data_df = pd.DataFrame(current_data_dict.values())
        current_df = pd.concat([current_ids_df, current_data_df], axis=1)
        current_df.sort_values(by=current_topology_columns, inplace=True)
        current_columns_sorted = current_topology_columns + current_outputs_names
        current_df = current_df.reindex_axis(current_columns_sorted, axis=1, copy=False)
        current_df.reset_index(drop=True, inplace=True)
        dataframes_dict[current_key] = current_df
    return dataframes_dict['hgz'], dataframes_dict['elements']


def from_MTG(g, hgz_inputs=None, element_inputs=None):
    """
    Convert a MTG to Growth-Wheat inputs.
    Use data in `hgz_inputs` and `element_inputs` if `g` is incomplete.
    
    :Parameters:

        - g (:class:`openalea.mtg.mtg.MTG`) - A MTG which contains the inputs
          of Growth-Wheat. These inputs are: :mod:`simulation.HGZ_INPUTS` and :mod:`simulation.ELEMENT_INPUTS`.

        - `hgz_inputs` (:class:`pandas.DataFrame`) - hidden growing zones dataframe, with one line by hidden growing zone.

        - `element_inputs` (:class:`pandas.DataFrame`) - elements dataframe, with one line by element.

    :Returns:
        The inputs of Growth-Wheat.

    :Returns Type:
        :class:`dict` of :class:`dict`

    .. seealso:: see :attr:`simulation.Simulation.inputs` for the structure of Growth-Wheat inputs.

    """
    all_hgz_dict = {}
    all_element_dict = {}
    hgz_to_prev_sheath_dict = {}
    
    if hgz_inputs is None:
        hgz_inputs = pd.DataFrame(columns=HGZ_TOPOLOGY_COLUMNS)
    if element_inputs is None:
        element_inputs = pd.DataFrame(columns=ELEMENT_TOPOLOGY_COLUMNS)
        
    hgz_inputs_grouped = hgz_inputs.groupby(HGZ_TOPOLOGY_COLUMNS)
    element_inputs_grouped = element_inputs.groupby(ELEMENT_TOPOLOGY_COLUMNS)
    sheath_inputs_grouped = element_inputs[element_inputs.organ == 'sheath'].groupby(ELEMENT_TOPOLOGY_COLUMNS)
    sheath_inputs_grouped_all_metamers = element_inputs[element_inputs.organ == 'sheath'].groupby(['plant', 'axis'])

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
                        
                    for element_vid in g.components_iter(organ_vid):
                        element_label = g.label(element_vid)
                        if element_label == 'HiddenElement': continue
                        element_id = (plant_index, axis_label, metamer_index, organ_label)
                        if element_id in element_inputs_grouped.groups:
                            element_inputs_group = element_inputs_grouped.get_group(element_id)
                            element_inputs_group_series = element_inputs_group.loc[element_inputs_group.first_valid_index(), simulation.ELEMENT_INPUTS]
                        else:
                            element_inputs_group_series = pd.Series()
                        element_inputs_from_mtg = g.get_vertex_property(element_vid)
                        element_inputs_dict = {}
                        is_valid_element = True
                        for element_input_name in simulation.ELEMENT_INPUTS:
                            
                            if element_input_name in element_inputs_from_mtg:
                                # use the input from the MTG
                                element_inputs_dict[element_input_name] = element_inputs_from_mtg[element_input_name]
                            else:
                                # use the input from the dataframe
                                if element_input_name in element_inputs_group_series:
                                    element_inputs_dict[element_input_name] = element_inputs_group_series[element_input_name]
                                else:
                                    is_valid_element = False
                                    break
                        if is_valid_element:
                            all_element_dict[element_id] = element_inputs_dict
                
                
                # hgz inputs and length of previous sheaths
                previous_metamer_vid = g.parent(metamer_vid)
                
                hgz_id = (plant_index, axis_label, metamer_index)
                if hgz_id in hgz_inputs_grouped.groups:
                    hgz_inputs_group = hgz_inputs_grouped.get_group(hgz_id)
                    hgz_inputs_group_series = hgz_inputs_group.loc[hgz_inputs_group.first_valid_index(), simulation.HGZ_INPUTS]
                else:
                    hgz_inputs_group_series = pd.Series()
                
                metamer_properties = g.get_vertex_property(metamer_vid)
                if 'hgz' in metamer_properties:
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
                        
                    # Length of the previous sheaths
                    if previous_metamer_vid is not None:
                        for previous_metamer_organ_vid in g.components_iter(previous_metamer_vid):
                            if g.label(previous_metamer_organ_vid) == 'sheath':
                                previous_sheath_vid = previous_metamer_organ_vid
                                break
                        else:
                            raise Exception('Error: no previous sheath found.')
                            
                        previous_sheath_properties = g.get_vertex_property(previous_sheath_vid)
                        if previous_sheath_properties['is_growing']:
                            previous_growing_sheath_ancestors_vid = []
                            for previous_metamer_ancestor_vid in g.Ancestors(previous_metamer_vid):
                                for previous_metamer_ancestor_organ_vid in g.components_iter(previous_metamer_ancestor_vid):
                                    if g.label(previous_metamer_ancestor_organ_vid) == 'sheath':
                                        previous_sheath_ancestor_vid = previous_metamer_ancestor_organ_vid
                                        break
                                else:
                                    raise Exception('Error: no previous sheath found.')
                                previous_sheath_ancestor_properties = g.get_vertex_property(previous_sheath_ancestor_vid)
                                if previous_sheath_ancestor_properties['is_growing']:
                                    previous_growing_sheath_ancestors_vid.append(previous_sheath_ancestor_vid)
                            if len(previous_growing_sheath_ancestors_vid) > 1:
                                raise Warning('Several previous growing sheaths found.') #: In wheat there is usually a single growing sheath emerged
                            previous_growing_sheaths_L = 0.0
                            for previous_growing_sheath_ancestor_vid in previous_growing_sheath_ancestors_vid:
                                previous_growing_sheaths_L += g.get_vertex_property(previous_growing_sheath_ancestor_vid)['length']
                                
                            previous_mature_sheath_vid = g.parent(previous_growing_sheath_ancestors_vid[-1])
                            if previous_mature_sheath_vid is not None:
                                previous_mature_sheath_L = g.get_vertex_property(previous_mature_sheath_vid)['length']
                            else:
                                previous_mature_sheath_L = 0
                                raise Warning('No previous mature sheath found.')
                        else:
                            previous_growing_sheaths_L = 0
                            previous_mature_sheath_L = previous_sheath_properties['length']
                            
                        hgz_to_prev_sheath_dict[(plant_index, axis_label, metamer_index)] = (previous_growing_sheaths_L, previous_mature_sheath_L) #TODO: distinguer emerge de total dans noms variables   
                        
                else: # use the dataframe
                    
                    # hgz inputs
                    hgz_inputs_group_dict = hgz_inputs_group_series.to_dict()
                    if set(hgz_inputs_group_dict).issuperset(simulation.HGZ_INPUTS):
                        all_hgz_dict[hgz_id] = hgz_inputs_group_dict
                
                    # Length of the previous sheaths
                    if previous_metamer_vid is not None:
                        if axis_id in sheath_inputs_grouped_all_metamers.groups:
                            sheath_inputs_group_all_metamers = sheath_inputs_grouped_all_metamers.get_group(axis_id) #: all sheaths of the axis (i.e. of all metamers of the axis)
                            for previous_metamer_index in reversed(range(1, metamer_index)):
                                previous_sheath_id = tuple(list(axis_id) + [previous_metamer_index, 'sheath'])
                                if sheath_inputs_grouped.groups.has_key(previous_sheath_id):
                                    previous_sheath = sheath_inputs_grouped.get_group(previous_sheath_id)
                                    if previous_sheath.loc[previous_sheath.first_valid_index(), 'is_growing']: #: the previous sheath is growing
                                        previous_growing_sheaths = sheath_inputs_group_all_metamers.loc[(sheath_inputs_group_all_metamers.metamer <= previous_metamer_index) & sheath_inputs_group_all_metamers.is_growing] #: all previous growing sheaths
                                        if len(previous_growing_sheaths) > 1:
                                            raise Warning('Several previous growing sheaths found.') #: In wheat there is usually a single growing sheath emerged
                                        previous_growing_sheaths_L = previous_growing_sheaths.length.sum() #: Sum of the previous growing sheath lengths (although only 1 previous growing sheath is expected, see Warning above)
                                        previous_mature_sheath_index = previous_growing_sheaths.first_valid_index() - 1 #: Found the first previous mature sheath
                                        if previous_mature_sheath_index in sheath_inputs_group_all_metamers.index:
                                            previous_mature_sheath_L = sheath_inputs_group_all_metamers.loc[previous_mature_sheath_index, 'length']
                                        else:
                                            previous_mature_sheath_L = 0
                                            raise Warning('No previous mature sheath found.')
                                    else:
                                        previous_growing_sheaths_L = 0
                                        previous_mature_sheath_L = previous_sheath.loc[previous_sheath.first_valid_index(), 'length']
                                    break
                            else:
                                raise Exception('Error: no previous sheath found.')
                    
                            hgz_to_prev_sheath_dict[(plant_index, axis_label, metamer_index)] = (previous_growing_sheaths_L, previous_mature_sheath_L) #TODO: distinguer emerge de total dans noms variables
                
    return {'hgz': all_hgz_dict, 'elements': all_element_dict, 'previous_sheaths_L': hgz_to_prev_sheath_dict}

def update_MTG(inputs, outputs, g, geometrical_model):
    """
    Update a MTG from Growth-Wheat inputs and outputs.

    :Parameters:
            - inputs (:class:`dict` of :class:`dict`) - Growth-Wheat inputs.
            These inputs are: :mod:`simulation.HGZ_INPUTS` and :mod:`simulation.ELEMENT_INPUTS`.
    
            - outputs (:class:`dict` of :class:`dict`) - Growth-Wheat outputs.
            These outputs are: :mod:`simulation.HGZ_OUTPUTS` and :mod:`simulation.ELEMENT_OUTPUTS`.

            - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to update from the inputs and outputs of Growth-Wheat.
            
            - `geometrical_model` (:func:`geometrical_model`) - The model which deals with geometry. 
              This model must implement a method `add_metamer(mtg, plant_index, axis_label)` to add 
              a metamer to a specific axis of a plant in a MTG. 
        
    .. seealso:: 
        
        * see :attr:`simulation.Simulation.inputs` for the structure of Growth-Wheat inputs.
        * see :attr:`simulation.Simulation.outputs` for the structure of Growth-Wheat outputs.

    """

    # add the properties if needed
    property_names = g.property_names()
    for growthwheat_data_name in set(simulation.HGZ_INPUTS_OUTPUTS + simulation.ELEMENT_INPUTS_OUTPUTS):
        if growthwheat_data_name not in property_names:
            g.add_property(growthwheat_data_name)
    if 'hgz' not in property_names:
        g.add_property('hgz')
    
    hgzs_inputs_dict = inputs['hgz']
    elements_inputs_dict = inputs['elements']
    hgzs_outputs_dict = outputs['hgz']
    elements_outputs_dict = outputs['elements']
    
    # add new metamer(s)
    axis_to_metamers_mapping = {}
    for metamer_id in sorted(hgzs_outputs_dict.iterkeys()):
        axis_id = (metamer_id[0], metamer_id[1])
        if axis_id not in axis_to_metamers_mapping:
            axis_to_metamers_mapping[axis_id] = []
        axis_to_metamers_mapping[axis_id].append(metamer_id)
        
    axis_to_old_metamers_mapping = {}
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        for axis_vid in g.components_iter(plant_vid):
            axis_label = g.label(axis_vid)
            metamer_ids = set([(plant_index, axis_label, g.index(metamer_vid)) for metamer_vid in g.components_iter(axis_vid)])
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
                if hgz_id in hgzs_outputs_dict:
                    hgz_outputs_dict = hgzs_outputs_dict[hgz_id]
                    metamer_properties = g.get_vertex_property(metamer_vid)
                    if 'hgz' not in metamer_properties:
                        g.property('hgz')[metamer_vid] = {}
                    
                    for hgz_data in (hgzs_inputs_dict[hgz_id], hgzs_outputs_dict[hgz_id]):
                        for hgz_data_name, hgz_data_value in hgz_data.iteritems():
                            if hgz_data_name not in ('leaf_Wlig', 'lamina_Lmax', 'leaf_Wmax'):
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
                        if organ_label == 'sheath':
                            g.property('diameter')[organ_vid] = mtg_organs_data_from_growthwheat_hgz_data['leaf_Wlig']
                        elif organ_label == 'blade':
                            g.property('shape_mature_length')[organ_vid] = mtg_organs_data_from_growthwheat_hgz_data['lamina_Lmax']
                            g.property('shape_max_width')[organ_vid] = mtg_organs_data_from_growthwheat_hgz_data['leaf_Wmax']
                            
                    exposed_element_id = (plant_index, axis_label, metamer_index, organ_label)
                    element_inputs_dict = elements_inputs_dict.get(exposed_element_id, {})
                    element_outputs_dict = elements_outputs_dict.get(exposed_element_id, {})
                    
                    for element_vid in g.components_iter(organ_vid):
                        element_label = g.label(element_vid)
                        if element_label == 'HiddenElement': continue
                        
                        for element_input_name in simulation.ELEMENT_INPUTS:
                            g.property(element_input_name)[element_vid] = element_inputs_dict.get(element_input_name)
                        for element_output_name in simulation.ELEMENT_OUTPUTS:
                            g.property(element_output_name)[element_vid] = element_outputs_dict.get(element_output_name)
                
                            