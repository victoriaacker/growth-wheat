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

import model
import logging
import warnings
import copy

#: the inputs needed by GrowthWheat
HZ_INPUTS = ['leaf_is_growing', 'hz_L', 'leaf_L', 'leaf_Lmax', 'leaf_Lem_prev', 'lamina_Lmax', 'sheath_Lmax', 'leaf_Wmax', 'SSLW', 'SSSW', 'leaf_is_emerged', 'sucrose', 'amino_acids', 'fructan', 'hz_mstruct']
ORGAN_INPUTS = ['visible_length', 'is_growing', 'final_hidden_length', 'length']

#: the outputs computed by GrowthWheat
HZ_OUTPUTS = ['leaf_is_growing', 'hz_L', 'leaf_L', 'delta_leaf_L', 'leaf_Lmax', 'leaf_Lem_prev', 'lamina_Lmax', 'sheath_Lmax', 'leaf_Wmax', 'SSLW', 'SSSW', 'leaf_is_emerged', 'sucrose', 'amino_acids', 'fructan', 'hz_mstruct']
ORGAN_OUTPUTS = ['visible_length', 'is_growing', 'final_hidden_length', 'length']

#: the inputs and outputs of GrowthWheat.
HZ_INPUTS_OUTPUTS = sorted(set(HZ_INPUTS + HZ_OUTPUTS))
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
        #:     {'hz': {(plant_index, axis_label, metamer_index): {hz_input_name: hz_input_value, ...}, ...},
        #:      'organs': {(plant_index, axis_label, metamer_index, organ_label): {organ_input_name: organ_input_value, ...}, ...},
        #:      'hz_L_calculation': {(plant_index, axis_label, metamer_index): {'previous_hz_length': previous_hz_length,
        #:                                                                       'previous_sheath_visible_length': previous_sheath_visible_length,
        #:                                                                       'previous_sheath_final_hidden_length': previous_sheath_final_hidden_length}, ...}}
        #: See :TODO?
        #: for more information about the inputs.
        self.inputs = {}

        #: The outputs of growth-Wheat.
        #:
        #: `outputs` is a dictionary of dictionaries:
        #:     {'hz': {(plant_index, axis_label, metamer_index): {hz_input_name: hz_input_value, ...}, ...},
        #:      'organs': {(plant_index, axis_label, metamer_index, organ_label): {organ_input_name: organ_input_value, ...}, ...},
        #:      'hz_L_calculation': {(plant_index, axis_label, metamer_index): {'previous_hz_length': previous_hz_length,
        #:                                                                       'previous_sheath_visible_length': previous_sheath_visible_length,
        #:                                                                       'previous_sheath_final_hidden_length': previous_sheath_final_hidden_length}, ...}}
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
        self.outputs.update({inputs_type: copy.deepcopy(all_inputs) for inputs_type, all_inputs in self.inputs.iteritems() if inputs_type in set(['hz', 'organs', 'hz_L_calculation'])})

        # Hidden zones
        all_hz_inputs = self.inputs['hz']
        all_hz_outputs = self.outputs['hz']

        # organs
        all_organs_inputs = self.inputs['organs']
        all_organs_outputs = self.outputs['organs']

        # Previous sheaths
        all_hz_L_calculation_inputs = self.inputs['hz_L_calculation']

        for hz_id, hz_inputs in all_hz_inputs.iteritems():
            curr_hz_outputs = all_hz_outputs[hz_id]

            # Found previous hidden zone TODO: a améliorer
            prev_hz_id = tuple(list(hz_id[:2]) + [hz_id[2] - 1])
            if prev_hz_id in all_hz_inputs:
                prev_leaf_emerged = all_hz_inputs[prev_hz_id]['leaf_is_emerged']
            else:
                prev_leaf_emerged = True

            # Update hidden zone length
            previous_hz_L = all_hz_L_calculation_inputs[hz_id]['previous_hz_length']
            previous_sheath_L = all_hz_L_calculation_inputs[hz_id]['previous_sheath_visible_length']
            previous_sheath_final_hidden_L = all_hz_L_calculation_inputs[hz_id]['previous_sheath_final_hidden_length']
            curr_hz_outputs['hz_L'] = model.calculate_hz_length(previous_hz_L, previous_sheath_L, previous_sheath_final_hidden_L)
            delta_hz_L = hz_inputs['hz_L'] - curr_hz_outputs['hz_L']

            if not prev_leaf_emerged: #: Before the emergence of the previous leaf. Exponential-like growth.
                ## delta leaf length
                delta_leaf_L = model.calculate_deltaL_preE(hz_inputs['sucrose'], hz_inputs['leaf_L'], hz_inputs['amino_acids'], hz_inputs['hz_mstruct'], self.delta_t)

            else: #: After the emergence of the previous leaf.
                ## delta leaf length
                delta_leaf_L = model.calculate_deltaL_postE(hz_inputs['leaf_L'], curr_hz_outputs['leaf_Lmax'], hz_inputs['sucrose'], self.delta_t)

                lamina_id = hz_id + tuple(['blade'])
                #: Lamina has not emerged
                if not curr_hz_outputs['leaf_is_emerged']:
                    #: Test of leaf emergence against hidden zone length. Assumes that a leaf cannot emerge before the previous one # TODO: besoin correction pour savoir à quel pas de temps exact??
                    curr_hz_outputs['leaf_is_emerged'] = model.calculate_leaf_emergence(hz_inputs['leaf_L'], curr_hz_outputs['hz_L'])
                    if curr_hz_outputs['leaf_is_emerged']: # Initialise lamina outputs
                        all_organs_outputs[lamina_id] = dict.fromkeys(ORGAN_OUTPUTS, 0)
                        all_organs_outputs[lamina_id]['is_growing'] = True

                        # Initialise variables for the next hidden zone
                        next_hz_id = tuple(list(hz_id[:2]) + [hz_id[2] + 1])
                        if next_hz_id in all_hz_inputs:
                            next_hz_inputs = all_hz_inputs[next_hz_id]
                            next_hz_outputs = all_hz_outputs[next_hz_id]
                            next_hz_outputs['leaf_Lem_prev'] = next_hz_inputs['leaf_L']                                                                  # Leaf length at the time of the emergence of the previous leaf
                            next_hz_outputs['leaf_Lmax'] = model.calculate_leaf_Lmax(next_hz_outputs['leaf_Lem_prev'])                                   # Final leaf length
                            sheath_lamina_ratio = model.calculate_SL_ratio(next_hz_id[2])                                                                 # Sheath:Lamina final length ratio
                            next_hz_outputs['lamina_Lmax'] = model.calculate_lamina_Lmax(next_hz_outputs['leaf_Lmax'], sheath_lamina_ratio)              # Final lamina length
                            next_hz_outputs['sheath_Lmax'] = model.calculate_sheath_Lmax(next_hz_outputs['leaf_Lmax'], next_hz_outputs['lamina_Lmax'])  # Final sheath length
                            next_hz_outputs['leaf_Wmax'] = model.calculate_leaf_Wmax(next_hz_outputs['lamina_Lmax'], next_hz_inputs['fructan'], next_hz_inputs['hz_mstruct'])        # Maximal leaf width
                            next_hz_outputs['SSLW'] = model.calculate_SSLW(next_hz_inputs['fructan'], next_hz_inputs['hz_mstruct'] )                   # Structural Specific Lamina Weight
                            next_hz_outputs['SSSW'] = model.calculate_SSSW(next_hz_outputs['SSLW'])                                                      # Structural Specific Sheath Weight
                            self.outputs['hz'][next_hz_id] = next_hz_outputs
                        else:
                            warnings.warn('No next hidden zone found for hz {}.'.format(hz_id))

                #: Lamina has emerged and is growing
                elif curr_hz_outputs['leaf_is_emerged'] and all_organs_inputs[lamina_id]['is_growing']:
                    curr_organ_outputs = all_organs_outputs[lamina_id]
                    ## Length of emerged lamina
                    lamina_L = min(model.calculate_lamina_L(hz_inputs['leaf_L'], curr_hz_outputs['hz_L']), curr_hz_outputs['lamina_Lmax'])
                    curr_organ_outputs['visible_length'] = lamina_L

                    # Test end of growth
                    if lamina_L >= curr_hz_outputs['lamina_Lmax']:
                        curr_organ_outputs['is_growing'] = False
                        curr_organ_outputs['final_hidden_length'] = 0
                        curr_organ_outputs['length'] = lamina_L
                        # Initialise sheath outputs
                        sheath_id = hz_id + tuple(['sheath'])
                        all_organs_outputs[sheath_id] = dict.fromkeys(ORGAN_OUTPUTS, 0)
                        all_organs_outputs[sheath_id]['is_growing'] = True

                    # Update of lamina outputs
                    self.outputs['organs'][lamina_id] = curr_organ_outputs

                # Mature lamina, growing sheath
                else:
                    sheath_id = hz_id + tuple(['sheath'])
                    curr_organ_outputs = all_organs_outputs[sheath_id]

                    ## Length of emerged sheath
                    lamina_L = self.outputs['organs'][hz_id + tuple(['blade'])]['visible_length']
                    sheath_L = min(model.calculate_sheath_L(hz_inputs['leaf_L'], curr_hz_outputs['hz_L'], lamina_L), curr_hz_outputs['sheath_Lmax'])
                    curr_organ_outputs['visible_length'] = sheath_L

                    #: Test end of growth
                    if hz_inputs['leaf_L'] >= curr_hz_outputs['leaf_Lmax']: #TODO:  hz_inputs['leaf_L'] ou  hz_outputs['leaf_L']
                        curr_organ_outputs['final_hidden_length'] = curr_hz_outputs['hz_L']
                        curr_organ_outputs['length'] = sheath_L + curr_organ_outputs['final_hidden_length'] #: Length of the mature sheath = visible length + hz length
                        curr_organ_outputs['is_growing'] = False
                        curr_hz_outputs['leaf_is_growing'] = False

                    # Update of sheath outputs
                    self.outputs['organs'][sheath_id] = curr_organ_outputs

            # Update of leaf outputs, TODO: attention aux valeurs negatives
            curr_hz_outputs['leaf_L'] = np.nanmin([curr_hz_outputs['leaf_Lmax'], (hz_inputs['leaf_L'] + delta_leaf_L)])
            curr_hz_outputs['delta_leaf_L'] = np.nanmin([delta_leaf_L, (curr_hz_outputs['leaf_Lmax'] - hz_inputs['leaf_L'])])

            if curr_hz_outputs['leaf_is_growing']:
                self.outputs['hz'][hz_id] = curr_hz_outputs
            else: # End of leaf growth
                del self.outputs['hz'][hz_id]