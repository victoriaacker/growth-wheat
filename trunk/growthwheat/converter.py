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
HIDDENZONE_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer']
ORGAN_TOPOLOGY_COLUMNS = ['plant', 'axis', 'metamer', 'organ'] # visible part of organs
ROOT_TOPOLOGY_COLUMNS = ['plant', 'axis']

#: the name of the organs representing a leaf
LEAF_ORGANS_NAMES = set(['sheath', 'blade'])


def from_dataframes(hiddenzone_inputs, organ_inputs, root_inputs):
    """
    Convert inputs/outputs from Pandas dataframe to Growth-Wheat format.

    :Parameters:

        - `hiddenzone_inputs` (:class:`pandas.DataFrame`) - Hidden zone inputs dataframe to convert, with one line by hidden zone.
        - `organ_inputs` (:class:`pandas.DataFrame`) - Exposed organ inputs dataframe to convert, with one line by organ.
        - `root_inputs` (:class:`pandas.DataFrame`) - Root inputs dataframe to convert, with one line by root.

    :Returns:
        The inputs in a dictionary.

    :Returns Type:
        :class:`dict` of :class:`dict`

    .. seealso:: see :attr:`simulation.Simulation.inputs` for the structure of Growth-Wheat inputs.

    """
    all_hiddenzone_dict = {}
    all_organ_dict = {}
    all_root_dict = {}
    hiddenzone_L_calculation_dict = {}
    hiddenzone_inputs_columns = hiddenzone_inputs.columns.difference(HIDDENZONE_TOPOLOGY_COLUMNS)
    organ_inputs_columns = organ_inputs.columns.difference(ORGAN_TOPOLOGY_COLUMNS)
    root_inputs_columns = root_inputs.columns.difference(ROOT_TOPOLOGY_COLUMNS)
    sheath_inputs_grouped = organ_inputs[organ_inputs.organ == 'sheath'].groupby(ORGAN_TOPOLOGY_COLUMNS)
    sheath_inputs_grouped_all_metamers = organ_inputs[organ_inputs.organ == 'sheath'].groupby(['plant', 'axis'])

    for organ_inputs_id, organ_inputs_group in organ_inputs.groupby(ORGAN_TOPOLOGY_COLUMNS):
        # organ
        organ_inputs_series = organ_inputs_group.loc[organ_inputs_group.first_valid_index()]
        organ_inputs_dict = organ_inputs_series[organ_inputs_columns].to_dict()
        all_organ_dict[organ_inputs_id] = organ_inputs_dict

    hiddenzone_inputs_grouped = hiddenzone_inputs.groupby(HIDDENZONE_TOPOLOGY_COLUMNS)
    for hiddenzone_inputs_id, hiddenzone_inputs_group in hiddenzone_inputs_grouped:
        # hiddenzone
        hiddenzone_inputs_series = hiddenzone_inputs_group.loc[hiddenzone_inputs_group.first_valid_index()]
        hiddenzone_inputs_dict = hiddenzone_inputs_series[hiddenzone_inputs_columns].to_dict()
        all_hiddenzone_dict[hiddenzone_inputs_id] = hiddenzone_inputs_dict

    for root_inputs_id, root_inputs_group in root_inputs.groupby(ROOT_TOPOLOGY_COLUMNS):
        # root
        root_inputs_series = root_inputs_group.loc[root_inputs_group.first_valid_index()]
        root_inputs_dict = root_inputs_series[root_inputs_columns].to_dict()
        all_root_dict[root_inputs_id] = root_inputs_dict

    return {'hiddenzone': all_hiddenzone_dict, 'organs': all_organ_dict, 'roots': all_root_dict}

def to_dataframes(data_dict):
    """
    Convert outputs from Growth-Wheat format to Pandas dataframe.

    :Parameters:

        - `data_dict` (:class:`dict`) - The outputs in Growth-Wheat format.

    :Returns:
        One dataframe for hiddenzone outputs and one dataframe for organ outputs.

    :Returns Type:
        :class:`tuple` of :class:`pandas.DataFrame`

    .. seealso:: see :attr:`simulation.Simulation.outputs` for the structure of Growth-Wheat outputs.

    """
    dataframes_dict = {}
    for (current_key, current_topology_columns, current_outputs_names) in (('hiddenzone', HIDDENZONE_TOPOLOGY_COLUMNS, simulation.HIDDENZONE_OUTPUTS),
                                                                           ('organs', ORGAN_TOPOLOGY_COLUMNS, simulation.ORGAN_OUTPUTS),
                                                                           ('roots', ROOT_TOPOLOGY_COLUMNS, simulation.ROOT_OUTPUTS)):
        current_data_dict = data_dict[current_key]
        current_ids_df = pd.DataFrame(current_data_dict.keys(), columns=current_topology_columns)
        current_data_df = pd.DataFrame(current_data_dict.values())
        current_df = pd.concat([current_ids_df, current_data_df], axis=1)
        current_df.sort_values(by=current_topology_columns, inplace=True)
        current_columns_sorted = current_topology_columns + current_outputs_names
        current_df = current_df.reindex_axis(current_columns_sorted, axis=1, copy=False)
        current_df.reset_index(drop=True, inplace=True)
        dataframes_dict[current_key] = current_df
    return dataframes_dict['hiddenzone'], dataframes_dict['organs'], dataframes_dict['roots']


def from_MTG(g):
    """
    Convert a MTG to Growth-Wheat inputs.

    :Parameters:

        - g (:class:`openalea.mtg.mtg.MTG`) - A MTG which contains the inputs
          of Growth-Wheat. These inputs are: :mod:`simulation.HIDDENZONE_INPUTS` and :mod:`simulation.ORGAN_INPUTS`.

    :Returns:
        The inputs of Growth-Wheat.

    :Returns Type:
        :class:`dict` of :class:`dict`

    .. seealso:: see :attr:`simulation.Simulation.inputs` for the structure of Growth-Wheat inputs.

    """
    all_hiddenzone_dict = {}
    all_organ_dict = {}
    all_root_dict = {}

    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        for axis_vid in g.components_iter(plant_vid):
            axis_label = g.label(axis_vid)

            # Roots
            axis_properties = g.get_vertex_property(axis_vid)
            if 'roots' in axis_properties:
                roots_id = (plant_index, axis_label)
                root_inputs_from_mtg = axis_properties['roots']
                roots_inputs_dict = {}

                is_valid_roots = True
                for roots_input_name in simulation.ROOT_INPUTS:
                    if roots_input_name in root_inputs_from_mtg:
                        # use the input from the MTG
                        roots_inputs_dict[roots_input_name] = root_inputs_from_mtg[roots_input_name]
                    else:
                        is_valid_roots = False
                        break
                if is_valid_roots:
                    all_root_dict[roots_id] = roots_inputs_dict

            for metamer_vid in g.components_iter(axis_vid):
                metamer_index = int(g.index(metamer_vid))
                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in LEAF_ORGANS_NAMES: continue

                    organ_inputs_from_mtg = g.get_vertex_property(organ_vid)

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
                    # temp
                    for element_vid in g.components_iter(organ_vid):
                        if g.get_vertex_property(element_vid)['label'] in ('StemElement', 'LeafElement1'):
                            organ_inputs_dict['green_area'] = g.get_vertex_property(element_vid)['green_area']

                    if is_valid_organ:
                        all_organ_dict[organ_id] = organ_inputs_dict

                metamer_properties = g.get_vertex_property(metamer_vid)
                if 'hiddenzone' in metamer_properties:
                    hiddenzone_id = (plant_index, axis_label, metamer_index)
                    hiddenzone_inputs_from_mtg = metamer_properties['hiddenzone']
                    hiddenzone_inputs_dict = {}

                    is_valid_hiddenzone = True
                    for hiddenzone_input_name in simulation.HIDDENZONE_INPUTS:
                        if hiddenzone_input_name in hiddenzone_inputs_from_mtg:
                            # use the input from the MTG
                            hiddenzone_inputs_dict[hiddenzone_input_name] = hiddenzone_inputs_from_mtg[hiddenzone_input_name]
                        else:
                            is_valid_hiddenzone = False
                            break
                    if is_valid_hiddenzone:
                        all_hiddenzone_dict[hiddenzone_id] = hiddenzone_inputs_dict

    return {'hiddenzone': all_hiddenzone_dict, 'organs': all_organ_dict, 'roots': all_root_dict}

def update_MTG(g, inputs=None, outputs=None):
    """
    Update a MTG from Growth-Wheat inputs and outputs.

    :Parameters:
            - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to update from the inputs and outputs of Growth-Wheat.



            - inputs (:class:`dict` of :class:`dict`) - Growth-Wheat inputs.
            These inputs are: :mod:`simulation.HIDDENZONE_INPUTS` and :mod:`simulation.ORGAN_INPUTS`.

            - outputs (:class:`dict` of :class:`dict`) - Growth-Wheat outputs.
            These outputs are: :mod:`simulation.HIDDENZONE_OUTPUTS` and :mod:`simulation.ORGAN_OUTPUTS`.

    .. seealso::

        * see :attr:`simulation.Simulation.inputs` for the structure of Growth-Wheat inputs.
        * see :attr:`simulation.Simulation.outputs` for the structure of Growth-Wheat outputs.

    """

    # add the properties if needed
    property_names = g.property_names()
    for growthwheat_data_name in set(simulation.HIDDENZONE_INPUTS_OUTPUTS + simulation.ORGAN_INPUTS_OUTPUTS):
        if growthwheat_data_name not in property_names:
            g.add_property(growthwheat_data_name)

    if inputs:
        hiddenzones_data_dict = inputs['hiddenzone']
        hiddenzones_data_names = simulation.HIDDENZONE_INPUTS
        organs_data_dict = inputs['organs']
        organs_data_names = simulation.ORGAN_INPUTS
        roots_data_dict = inputs['roots']
        roots_data_names = simulation.ROOT_INPUTS
    elif outputs:
        hiddenzones_data_dict = outputs['hiddenzone']
        hiddenzones_data_names = simulation.HIDDENZONE_OUTPUTS
        organs_data_dict = outputs['organs']
        organs_data_names = simulation.ORGAN_OUTPUTS
        roots_data_dict = outputs['roots']
        roots_data_names = simulation.ROOT_OUTPUTS
    else:
        raise Exception('Both inputs and outputs were found to be None')

    # update the properties of the MTG
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        for axis_vid in g.components_iter(plant_vid):
            axis_label = g.label(axis_vid)

            root_id = (plant_index, axis_label)
            if root_id in roots_data_dict:
                root_outputs_dict = roots_data_dict[root_id]
                axis_properties = g.get_vertex_property(axis_vid)
                if 'roots' not in axis_properties:
                    g.property('roots')[axis_vid] = {}
                for root_data in (roots_data_dict[root_id], roots_data_dict[root_id]):
                    for root_data_name, root_data_value in root_data.iteritems():
                            g.property('roots')[axis_vid][root_data_name] = root_data_value

            for metamer_vid in g.components_iter(axis_vid):
                metamer_index = int(g.index(metamer_vid))

                hiddenzone_id = (plant_index, axis_label, metamer_index)
                if hiddenzone_id in hiddenzones_data_dict:
                    hiddenzone_outputs_dict = hiddenzones_data_dict[hiddenzone_id]
                    metamer_properties = g.get_vertex_property(metamer_vid)
                    if 'hiddenzone' not in metamer_properties:
                        g.property('hiddenzone')[metamer_vid] = {}

                    for hiddenzone_data in (hiddenzones_data_dict[hiddenzone_id], hiddenzones_data_dict[hiddenzone_id]):
                        for hiddenzone_data_name, hiddenzone_data_value in hiddenzone_data.iteritems():
                                g.property('hiddenzone')[metamer_vid][hiddenzone_data_name] = hiddenzone_data_value

                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in LEAF_ORGANS_NAMES: continue

                    organ_id = (plant_index, axis_label, metamer_index, organ_label)
                    organ_data_dict = organs_data_dict.get(organ_id, {})
                    if organ_data_dict:
                        for organ_data_name in organs_data_names:
                            g.property(organ_data_name)[organ_vid] = organ_data_dict.get(organ_data_name)

                            # Write properties in element scale #: TODO: voir ce qui ce passe avec +sieurs élements
                            for element_vid in g.components_iter(organ_vid):
                                if g.get_vertex_property(element_vid)['label'] in ('StemElement', 'LeafElement1'):
                                    g.property(organ_data_name)[element_vid] = organ_data_dict.get(organ_data_name)