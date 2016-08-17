# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    growthwheat.simulation
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`growthwheat.simulation` is the front-end to run the CN-Wheat model.

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

import numpy as np
import pandas as pd

import model, parameters
import logging
import copy


#: the inputs needed by GrowthWheat
HGZ_INPUTS = ['leaf_is_growing', 'hgz_L', 'leaf_L', 'leaf_Lmax', 'leaf_Lem_prev', 'lamina_Lmax', 'sheath_Lmax', 'leaf_Lwmax', 'leaf_Wmax', 'leaf_Wlig', 'SSLW', 'SSSW', 't_prev_leaf_emerged', 'leaf_is_emerged', 'sucrose', 'amino_acids', 'fructan', 'hgz_mstruct']
ORGAN_INPUTS = ['length', 'width', 'area', 'is_growing', 'mstruct', 'sucrose', 'amino_acids', 'fructan']

#: the outputs computed by GrowthWheat
HGZ_OUTPUTS = ['leaf_is_growing', 'hgz_L', 'leaf_L', 'leaf_Lmax', 'leaf_Lem_prev', 'lamina_Lmax', 'sheath_Lmax', 'leaf_Lwmax', 'leaf_Wmax', 'leaf_Wlig', 'SSLW', 'SSSW', 't_prev_leaf_emerged', 'leaf_is_emerged', 'sucrose', 'amino_acids', 'fructan', 'hgz_mstruct']
ORGAN_OUTPUTS = ['length', 'width', 'area', 'is_growing', 'mstruct', 'sucrose', 'amino_acids', 'fructan']

#: the inputs and outputs of GrowthWheat.
HGZ_INPUTS_OUTPUTS = sorted(set(HGZ_INPUTS + HGZ_OUTPUTS))
ORGAN_INPUTS_OUTPUTS = sorted(set(ORGAN_INPUTS + ORGAN_OUTPUTS))


class SimulationError(Exception): pass
class SimulationRunError(SimulationError): pass


class Simulation(object):
    """The Simulation class permits to initialize and run a simulation.
    """

    def __init__(self, delta_t=1):

        #: The inputs of growth-Wheat.
        #:
        #: `inputs` is a dictionary of dictionaries:
        #:     {'hgz': {(plant_index, axis_label, metamer_index): {hgz_input_name: hgz_input_value, ...}, ...},
        #:      'organs': {(plant_index, axis_label, metamer_index, organ_label): {organ_input_name: organ_input_value, ...}, ...},
        #:      'previous_sheaths_L': {(plant_index, axis_label, metamer_index): (previous_growing_sheath_L, previous_mature_sheath_L), ...}}
        #: See :TODO?
        #: for more information about the inputs.
        self.inputs = {}

        #: The outputs of growth-Wheat.
        #:
        #: `outputs` is a dictionary of dictionaries:
        #:     {'hgz': {(plant_index, axis_label, metamer_index): {hgz_input_name: hgz_input_value, ...}, ...},
        #:      'organs': {(plant_index, axis_label, metamer_index, organ_label): {organ_input_name: organ_input_value, ...}, ...},}
        #: See :TODO?
        #: for more information about the inputs.
        self.outputs = {}

        #: the delta t of the simulation (in seconds)
        self.delta_t = delta_t


    def initialize(self, inputs):
        """
        Initialize :attr:`inputs` from `inputs`.

        :Parameters:

            - `inputs` (:class:`dict`)
              `inputs` must be a dictionary with the same structure as :attr:`inputs`.
        """
        self.inputs.clear()
        self.inputs.update(inputs)


    def run(self):
        # Copy the inputs into the output dict
        self.outputs.update({inputs_type: copy.deepcopy(all_inputs) for inputs_type, all_inputs in self.inputs.iteritems() if inputs_type in set(['hgz', 'organs'])})

        # Hidden growing zones
        all_hgz_inputs = self.inputs['hgz']
        all_hgz_outputs = self.outputs['hgz']

        # organs
        all_organs_inputs = self.inputs['organs']
        all_organs_outputs = self.outputs['organs']

        # Previous sheaths
        all_prev_sheath_inputs = self.inputs['previous_sheaths_L']

        for hgz_id, hgz_inputs in all_hgz_inputs.iteritems():
            curr_hgz_outputs = all_hgz_outputs[hgz_id]

            # Initialisation of the exports towards the growing lamina or sheath
            delta_lamina_mstruct, delta_sheath_mstruct, export_sucrose, export_amino_acids = 0, 0, 0, 0

            # Found previous hidden growing zone
            prev_hgz_id = tuple(list(hgz_id[:2]) + [hgz_id[2] - 1])
            if prev_hgz_id in all_hgz_inputs:
                prev_leaf_emerged = all_hgz_inputs[prev_hgz_id]['leaf_is_emerged']
            else:
                prev_leaf_emerged = True

            # Found previous sheath
            previous_growing_sheath_L = all_prev_sheath_inputs[hgz_id][0]
            previous_mature_sheath_L = all_prev_sheath_inputs[hgz_id][1]

            # Update hidden growing zone length
            curr_hgz_outputs['hgz_L'] = model.calculate_hgz_length(previous_growing_sheath_L, previous_mature_sheath_L)
            delta_hgz_L = hgz_inputs['hgz_L'] - curr_hgz_outputs['hgz_L']

            if not prev_leaf_emerged: #: Before the emergence of the previous leaf. Exponential-like growth.
                ## delta leaf length
                delta_leaf_L = model.calculate_deltaL_preE(hgz_inputs['sucrose'], hgz_inputs['leaf_L'], hgz_inputs['amino_acids'], hgz_inputs['hgz_mstruct'], self.delta_t)
                ## delta hidden growing zone mstruct (=delta leaf mstruct at this stage)
                delta_hgz_mstruct = model.calculate_delta_mstruct_preE(hgz_inputs['leaf_L'], delta_leaf_L)

            else: #: After the emergence of the previous leaf.
                if hgz_inputs['t_prev_leaf_emerged'] == 0:
                    # Stores values once the previous leaf has emerged
                    curr_hgz_outputs['leaf_Lem_prev'] = hgz_inputs['leaf_L']                                                                       # Leaf length at the time of the emergence of the previous leaf
                    curr_hgz_outputs['leaf_Lmax'] = model.calculate_leaf_Lmax(curr_hgz_outputs['leaf_Lem_prev'])                                   # Final leaf length
                    sheath_lamina_ratio = model.calculate_SL_ratio(hgz_id[2])                                                                      # Sheath:Lamina final length ratio
                    curr_hgz_outputs['lamina_Lmax'] = model.calculate_lamina_Lmax(curr_hgz_outputs['leaf_Lmax'], sheath_lamina_ratio)              # Final lamina length
                    curr_hgz_outputs['sheath_Lmax'] = model.calculate_sheath_Lmax(curr_hgz_outputs['leaf_Lmax'], curr_hgz_outputs['lamina_Lmax'])  # Final sheath length
                    curr_hgz_outputs['leaf_Lwmax'] = model.calculate_leaf_Lwmax(curr_hgz_outputs['lamina_Lmax'])                                   # Position of the maximal width along the leaf
                    curr_hgz_outputs['leaf_Wmax'] = model.calculate_leaf_Wmax(curr_hgz_outputs['lamina_Lmax'], hgz_inputs['fructan'], hgz_inputs['hgz_mstruct'])        # Maximal leaf width
                    curr_hgz_outputs['leaf_Wlig'] = model.calculate_leaf_Wlig(curr_hgz_outputs['leaf_Wmax'])                                       # Lamina width at ligule position, also used to determine sheath width
                    curr_hgz_outputs['SSLW'] = model.calculate_SSLW(hgz_inputs['fructan'], hgz_inputs['hgz_mstruct'] )                             # Structural Specific Lamina Weight
                    curr_hgz_outputs['SSSW'] = model.calculate_SSSW(curr_hgz_outputs['SSLW'])                                                      # Structural Specific Sheath Weight

                ## delta leaf length
                delta_leaf_L = model.calculate_deltaL_postE(hgz_inputs['leaf_L'], curr_hgz_outputs['leaf_Lmax'], hgz_inputs['sucrose'], self.delta_t)

                lamina_id = hgz_id + tuple(['blade'])
                #: Lamina has not emerged
                if not curr_hgz_outputs['leaf_is_emerged']:
                    ## delta hgz mstruct (=delta leaf mstruct at this stage)
                    delta_hgz_mstruct = model.calculate_delta_mstruct_preE(hgz_inputs['leaf_L'], delta_leaf_L)

                    #: Test of leaf emergence against hidden growing zone length. Assumes that a leaf cannot emerge before the previous one # TODO: besoin correction pour savoir à quel pas de temps exact??
                    curr_hgz_outputs['leaf_is_emerged'] = model.calculate_leaf_emergence(hgz_inputs['leaf_L'], curr_hgz_outputs['hgz_L'])
                    if curr_hgz_outputs['leaf_is_emerged']: # Initialise lamina outputs
                        all_organs_outputs[lamina_id] = dict.fromkeys(ORGAN_OUTPUTS, 0)
                        all_organs_outputs[lamina_id]['is_growing'] = True

                #: Lamina has emerged and is growing
                elif curr_hgz_outputs['leaf_is_emerged'] and all_organs_inputs[lamina_id]['is_growing']:
                    curr_organ_outputs = all_organs_outputs[lamina_id]
                    ## Length of emerged lamina
                    lamina_L = model.calculate_lamina_L(hgz_inputs['leaf_L'], curr_hgz_outputs['hgz_L'])
                    curr_organ_outputs['length'] = lamina_L

                    ## Width of emerged lamina
                    lamina_W = model.calculate_lamina_W(lamina_L, curr_hgz_outputs['leaf_Lwmax'], curr_hgz_outputs['leaf_Wmax'], curr_hgz_outputs['leaf_Wlig'], curr_hgz_outputs['lamina_Lmax'])
                    curr_organ_outputs['width'] = lamina_W

                    ## Delta of lamina width
                    # TODO: CE TEST POURRAIT ETRE REMPLACE PAR delta_lamina_W = lamina_W(t) - lamina_W(t-1)
                    delta_lamina_W  = model.calculate_delta_lamina_W(lamina_L, curr_hgz_outputs['leaf_Lwmax'], curr_hgz_outputs['leaf_Wmax'], delta_leaf_L, curr_hgz_outputs['lamina_Lmax'], curr_hgz_outputs['leaf_Wlig'])

                    ## Delta of lamina area
                    delta_lamina_area = model.calculate_delta_lamina_area(delta_leaf_L, all_organs_inputs[lamina_id]['width'], delta_lamina_W) # TODO: OK?
                    ## Delta mstruct of the emerged lamina
                    delta_lamina_mstruct = model.calculate_delta_mstruct_postE(delta_lamina_area, curr_hgz_outputs['SSLW'])
                    ## Export of sucrose from hgz towards emerged lamina
                    export_sucrose = model.calculate_export_sucrose(delta_lamina_mstruct, hgz_inputs['sucrose'], hgz_inputs['hgz_mstruct'])
                    ## Export of amino acids from hgz towards emerged lamina
                    export_amino_acids = model.calculate_export_amino_acids(delta_lamina_mstruct, hgz_inputs['amino_acids'], hgz_inputs['hgz_mstruct'])

                    ## Delta mstruct of the hgz
                    delta_hgz_mstruct = model.calculate_delta_mstruct_preE(hgz_inputs['hgz_L'], delta_hgz_L)

                    # Test end of growth
                    if lamina_L >= curr_hgz_outputs['lamina_Lmax']:
                        curr_organ_outputs['is_growing'] = False
                        # Initialise sheath outputs
                        sheath_id = hgz_id + tuple(['sheath'])
                        all_organs_outputs[sheath_id] = dict.fromkeys(ORGAN_OUTPUTS, 0)
                        all_organs_outputs[sheath_id]['is_growing'] = True

                    # Update of lamina outputs
                    curr_organ_outputs['area'] += delta_lamina_area
                    curr_organ_outputs['mstruct'] += delta_lamina_mstruct
                    curr_organ_outputs['sucrose'] += export_sucrose
                    curr_organ_outputs['amino_acids'] += export_amino_acids
                    self.outputs['organs'][lamina_id] = curr_organ_outputs

                # Mature lamina, growing sheath
                elif curr_hgz_outputs['leaf_is_emerged'] and not all_organs_inputs[lamina_id]['is_growing']:
                    sheath_id = hgz_id + tuple(['sheath'])
                    curr_organ_outputs = all_organs_outputs[sheath_id]

                    ## Length of emerged sheath
                    lamina_L = self.outputs['organs'][hgz_id + tuple(['lamina'])]['length']
                    curr_organ_outputs['length'] = model.calculate_sheath_L(hgz_inputs['leaf_L'], curr_hgz_outputs['hgz_L'], lamina_L)

                    ## Sheath width
                    curr_organ_outputs['width'] = curr_hgz_outputs['leaf_Wlig'] # TODO: attention la largeur n'est pas vraiment égale à Wlig

                    ## Sheath area
                    curr_organ_outputs['area'] = model.calculate_sheath_area(curr_organ_outputs['length'], curr_hgz_outputs['leaf_Wlig']) # TODO: developee ou projetee??, faire 1 relation pour avoir le diametre apparent (différent de Wlig) et donc l'interception
                    delta_sheath_area = curr_organ_outputs['area'] - all_organs_inputs[sheath_id]['area']

                    ## Delta mstruct of the emerged sheath
                    delta_sheath_mstruct = model.calculate_delta_mstruct_postE(delta_sheath_area, curr_hgz_outputs['SSSW'])
                    ## Export of sucrose from hgz towards emerged sheath
                    export_sucrose = model.calculate_export_sucrose(delta_sheath_mstruct, hgz_inputs['sucrose'], hgz_inputs['hgz_mstruct'])
                    ## Export of amino acids from hgz towards emerged sheath
                    export_amino_acids = model.calculate_export_amino_acids(delta_sheath_mstruct, hgz_inputs['amino_acids'], hgz_inputs['hgz_mstruct'])

                    ## Delta mstruct of the hgz
                    delta_hgz_mstruct = model.calculate_delta_mstruct_preE(hgz_inputs['hgz_L'], delta_hgz_L)

                    #: Test end of growth
                    if hgz_inputs['leaf_L'] >= curr_hgz_outputs['leaf_Lmax']: #TODO:  hgz_inputs['leaf_L'] ou  hgz_outputs['leaf_L']
                        curr_organ_outputs['is_growing'] = False
                        curr_hgz_outputs['leaf_is_growing'] = False

                    # Update of sheath outputs
                    curr_organ_outputs['mstruct'] += delta_sheath_mstruct
                    curr_organ_outputs['sucrose'] += export_sucrose
                    curr_organ_outputs['amino_acids'] += export_amino_acids
                    self.outputs['organs'][sheath_id] = curr_organ_outputs

                ## Increment t_prev_leaf_emerged
                curr_hgz_outputs['t_prev_leaf_emerged'] = model.calculate_t_prev_leaf_emerged(hgz_inputs['t_prev_leaf_emerged'], self.delta_t)


            # End of leaf calculations
            sucrose_consumption_mstruct = model.calculate_s_mstruct_sucrose(delta_hgz_mstruct, delta_lamina_mstruct, delta_sheath_mstruct) # Consumption of sucrose due to mstruct growth
            Respiration = model.Respiration(delta_hgz_mstruct + delta_lamina_mstruct + delta_sheath_mstruct)   # Respiration growth, temporary
            AA_consumption_mstruct = model.calculate_s_mstruct_amino_acids(delta_hgz_mstruct, delta_lamina_mstruct, delta_sheath_mstruct)  # Consumption of amino acids due to mstruct growth
            # Update of leaf outputs, TODO: attention aux valeurs negatives
            curr_hgz_outputs['leaf_L'] = min(curr_hgz_outputs['leaf_Lmax'], (hgz_inputs['leaf_L'] + delta_leaf_L))
            curr_hgz_outputs['hgz_mstruct'] = hgz_inputs['hgz_mstruct'] + delta_hgz_mstruct
            curr_hgz_outputs['sucrose'] = hgz_inputs['sucrose'] - (sucrose_consumption_mstruct + Respiration + export_sucrose)
            curr_hgz_outputs['amino_acids'] = hgz_inputs['amino_acids'] - (AA_consumption_mstruct + export_amino_acids)

            if curr_hgz_outputs['leaf_is_growing']:
                self.outputs['hgz'][hgz_id] = curr_hgz_outputs
            else: # End of leaf growth, hgz compartments allocated to the sheath.
                # TODO: separer enclosed/exposed
                del self.outputs['hgz'][hgz_id]
                all_organs_outputs[sheath_id]['mstruct'] = curr_hgz_outputs['mstruct']
                all_organs_outputs[sheath_id]['sucrose'] = curr_hgz_outputs['sucrose']
                all_organs_outputs[sheath_id]['amino_acids'] = curr_hgz_outputs['amino_acids']
                all_organs_outputs[sheath_id]['fructan'] = curr_hgz_outputs['fructan']