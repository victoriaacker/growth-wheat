# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    growthwheat.simulation
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`growthwheat.simulation`.

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

import model
import copy

from respiwheat.model import RespirationModel

#: the inputs needed by GrowthWheat
HIDDENZONE_INPUTS = ['leaf_is_growing','internode_is_growing','leaf_L', 'delta_leaf_L', 'internode_L','delta_internode_L', 'leaf_dist_to_emerge', 'delta_leaf_dist_to_emerge','internode_dist_to_emerge','delta_internode_dist_to_emerge', 'SSLW', 'SSSW', 'leaf_is_emerged','SSINW','internode_is_visible', 'sucrose', 'amino_acids', 'fructan','proteins','leaf_enclosed_mstruct','leaf_enclosed_Nstruct','internode_mstruct', 'internode_Nstruct','mstruct']
ORGAN_INPUTS = ['is_growing', 'mstruct', 'green_area', 'sucrose', 'amino_acids', 'fructan','proteins','Nstruct']
ROOT_INPUTS = ['sucrose', 'amino_acids', 'mstruct', 'Nstruct']

#: the outputs computed by GrowthWheat
HIDDENZONE_OUTPUTS = ['sucrose', 'amino_acids','fructan','proteins', 'leaf_enclosed_mstruct', 'leaf_enclosed_Nstruct', 'internode_mstruct', 'internode_Nstruct','mstruct', 'Respi_growth', 'sucrose_consumption_mstruct', 'AA_consumption_mstruct']
ORGAN_OUTPUTS = ['sucrose', 'amino_acids','fructan','proteins', 'mstruct', 'Nstruct']
ROOT_OUTPUTS = ['sucrose', 'amino_acids', 'mstruct', 'Nstruct', 'Respi_growth', 'rate_mstruct_growth']

#: the inputs and outputs of GrowthWheat.
HIDDENZONE_INPUTS_OUTPUTS = sorted(set(HIDDENZONE_INPUTS + HIDDENZONE_OUTPUTS))
ORGAN_INPUTS_OUTPUTS = sorted(set(ORGAN_INPUTS + ORGAN_OUTPUTS))
ROOT_INPUTS_OUTPUTS = sorted(set(ROOT_INPUTS + ROOT_OUTPUTS))


class SimulationError(Exception): pass
class SimulationRunError(SimulationError): pass


class Simulation(object):
    """The Simulation class permits to initialize and run a simulation.
    """

    def __init__(self, delta_t=1):

        #: The inputs of growth-Wheat.
        #:
        #: `inputs` is a dictionary of dictionaries:
        #:     {'hiddenzone': {(plant_index, axis_label, metamer_index): {hiddenzone_input_name: hiddenzone_input_value, ...}, ...},
        #:      'organs': {(plant_index, axis_label, metamer_index, organ_label): {organ_input_name: organ_input_value, ...}, ...},
        #:      'roots': {(plant_index, axis_label): {root_input_name: root_input_value, ...}, ...}}
        #: See :TODO?
        #: for more information about the inputs.
        self.inputs = {}

        #: The outputs of growth-Wheat.
        #:
        #: `outputs` is a dictionary of dictionaries:
        #:     {'hiddenzone': {(plant_index, axis_label, metamer_index): {hiddenzone_input_name: hiddenzone_input_value, ...}, ...},
        #:      'organs': {(plant_index, axis_label, metamer_index, organ_label): {organ_input_name: organ_input_value, ...}, ...}
        #:      'roots': {(plant_index, axis_label): {root_input_name: root_input_value, ...}, ...}}
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
        self.outputs.update({inputs_type: copy.deepcopy(all_inputs) for inputs_type, all_inputs in self.inputs.iteritems() if inputs_type in set(['hiddenzone', 'organs', 'roots'])})

        # Hidden growing zones
        all_hiddenzone_inputs = self.inputs['hiddenzone']
        all_hiddenzone_outputs = self.outputs['hiddenzone']

        # organs
        all_organs_inputs = self.inputs['organs']
        all_organs_outputs = self.outputs['organs']

        # roots
        all_roots_inputs = self.inputs['roots']
        all_roots_outputs = self.outputs['roots']

        # hidden zones and organs
        for hiddenzone_id, hiddenzone_inputs in all_hiddenzone_inputs.iteritems():

            curr_hiddenzone_outputs = all_hiddenzone_outputs[hiddenzone_id]

            # Initialisation of the exports towards the growing lamina or sheath
            delta_lamina_mstruct, delta_sheath_mstruct, delta_lamina_Nstruct, delta_sheath_Nstruct, export_sucrose, export_amino_acids = 0, 0, 0, 0, 0, 0

            # Delta Growth internode

            if not hiddenzone_inputs['internode_is_visible']: #: Internode is not visible
               ## delta mstruct of the internode
               delta_internode_enclosed_mstruct = model.calculate_delta_internode_enclosed_mstruct(hiddenzone_inputs['internode_L'], hiddenzone_inputs['delta_internode_L'])
               ## delta Nstruct of the internode
               delta_internode_enclosed_Nstruct = model.calculate_delta_Nstruct(delta_internode_enclosed_mstruct)
            else : #: Internode is visible
                ## delta mstruct of the enclosed internode
                delta_internode_enclosed_mstruct = model.calculate_delta_internode_enclosed_mstruct(hiddenzone_inputs['internode_dist_to_emerge'], hiddenzone_inputs['delta_internode_dist_to_emerge'])
                ## delta Nstruct of the enclosed internode
                delta_internode_enclosed_Nstruct = model.calculate_delta_Nstruct(delta_internode_enclosed_mstruct)

                internode_id = hiddenzone_id + tuple(['internode'])
                curr_organ_inputs = all_organs_inputs[internode_id]
                curr_organ_outputs = all_organs_outputs[internode_id]
                ## Delta mstruct of the emerged internode
                delta_internode_mstruct = model.calculate_delta_emerged_tissue_mstruct(hiddenzone_inputs['SSINW'], curr_organ_inputs['mstruct'], curr_organ_inputs['green_area'])
                ## Delta Nstruct of the emerged internode
                delta_internode_Nstruct = model.calculate_delta_Nstruct(delta_internode_mstruct)
                ## Export of sucrose from hiddenzone towards emerged internode
                export_sucrose = model.calculate_export_sucrose(delta_internode_mstruct, hiddenzone_inputs['sucrose'], hiddenzone_inputs['mstruct'])
                ## Export of amino acids from hiddenzone towards emerged internode
                export_amino_acids = model.calculate_export_amino_acids(delta_internode_mstruct, hiddenzone_inputs['amino_acids'], hiddenzone_inputs['mstruct'])

               # Update of internode outputs
                curr_organ_outputs['mstruct'] += delta_internode_mstruct
                curr_organ_outputs['Nstruct'] += delta_internode_Nstruct
                curr_organ_outputs['sucrose'] += export_sucrose
                curr_organ_outputs['amino_acids'] += export_amino_acids
                self.outputs['organs'][internode_id] = curr_organ_outputs

            # Delta Growth leaf

            if not hiddenzone_inputs['leaf_is_emerged']: #: leaf has not emerged
                ## delta mstruct of the hidden leaf
                delta_leaf_enclosed_mstruct = model.calculate_delta_leaf_enclosed_mstruct(hiddenzone_inputs['leaf_L'], hiddenzone_inputs['delta_leaf_L'])
                ## delta Nstruct of the hidden leaf
                delta_leaf_enclosed_Nstruct = model.calculate_delta_Nstruct(delta_leaf_enclosed_mstruct)

            else: #: leaf has emerged
                ## delta mstruct of the enclosed leaf (which length is assumed to equal the length of the pseudostem)
                delta_leaf_enclosed_mstruct = model.calculate_delta_leaf_enclosed_mstruct(hiddenzone_inputs['leaf_dist_to_emerge'], hiddenzone_inputs['delta_leaf_dist_to_emerge'])
                ## delta Nstruct of the enclosed en leaf
                delta_leaf_enclosed_Nstruct = model.calculate_delta_Nstruct(delta_leaf_enclosed_mstruct)

                share_leaf_getting_visible_mstruct, share_leaf_getting_visible_Nstruct = 0, 0
                if delta_leaf_enclosed_mstruct < 0: #: Pseudostem length decreased because internode elongation. We assume that internode n starts to elongate when sheath n-1 is mature.
                   delta_leaf_enclosed_mstruct, delta_leaf_enclosed_Nstruct = 0, 0
                   share_leaf_getting_visible_mstruct = abs(delta_leaf_enclosed_mstruct) / hiddenzone_inputs['mstruct']
                   share_leaf_getting_visible_Nstruct = abs(delta_leaf_enclosed_Nstruct) / hiddenzone_inputs['leaf_enclosed_Nstruct']

                lamina_id = hiddenzone_id + tuple(['blade'])
                #: Lamina is growing
                if all_organs_inputs[lamina_id]['is_growing']:
                    curr_organ_inputs = all_organs_inputs[lamina_id]
                    curr_organ_outputs = all_organs_outputs[lamina_id]
                    ## Delta mstruct of the emerged lamina
                    delta_lamina_mstruct = model.calculate_delta_emerged_tissue_mstruct(hiddenzone_inputs['SSLW'], curr_organ_inputs['mstruct'], curr_organ_inputs['green_area'])
                    ## Delta Nstruct of the emerged lamina
                    delta_lamina_Nstruct = model.calculate_delta_Nstruct(delta_lamina_mstruct)
                    ## Export of sucrose from hiddenzone towards emerged lamina
                    export_sucrose = model.calculate_export_sucrose(delta_lamina_mstruct, hiddenzone_inputs['sucrose'], hiddenzone_inputs['mstruct'])
                    ## Export of amino acids from hiddenzone towards emerged lamina
                    export_amino_acids = model.calculate_export_amino_acids(delta_lamina_mstruct, hiddenzone_inputs['amino_acids'], hiddenzone_inputs['mstruct'])

                    # Update of lamina outputs
                    curr_organ_outputs['mstruct'] += delta_lamina_mstruct + share_leaf_getting_visible_mstruct * hiddenzone_inputs['mstruct']
                    curr_organ_outputs['Nstruct'] += delta_lamina_Nstruct + share_leaf_getting_visible_Nstruct * hiddenzone_inputs['leaf_enclosed_Nstruct']
                    curr_organ_outputs['sucrose'] += export_sucrose + share_leaf_getting_visible_mstruct * hiddenzone_inputs['sucrose']
                    curr_organ_outputs['amino_acids'] += export_amino_acids + share_leaf_getting_visible_mstruct * hiddenzone_inputs['amino_acids']
                    curr_organ_outputs['fructan'] += share_leaf_getting_visible_mstruct * hiddenzone_inputs['fructan']
                    curr_organ_outputs['proteins'] += share_leaf_getting_visible_mstruct * hiddenzone_inputs['proteins']
                    self.outputs['organs'][lamina_id] = curr_organ_outputs

                else: #: Mature lamina, growing sheath
                    sheath_id = hiddenzone_id + tuple(['sheath'])
                    curr_organ_inputs = all_organs_inputs[sheath_id]
                    curr_organ_outputs = all_organs_outputs[sheath_id]
                    ## Delta mstruct of the emerged sheath
                    delta_sheath_mstruct = model.calculate_delta_emerged_tissue_mstruct(hiddenzone_inputs['SSSW'], curr_organ_inputs['mstruct'], curr_organ_inputs['green_area'])
                    ## Delta Nstruct of the emerged sheath
                    delta_sheath_Nstruct = model.calculate_delta_Nstruct(delta_sheath_mstruct)
                    ## Export of sucrose from hiddenzone towards emerged sheath
                    export_sucrose = model.calculate_export_sucrose(delta_sheath_mstruct, hiddenzone_inputs['sucrose'], hiddenzone_inputs['mstruct'])
                    ## Export of amino acids from hiddenzone towards emerged sheath
                    export_amino_acids = model.calculate_export_amino_acids(delta_sheath_mstruct, hiddenzone_inputs['amino_acids'], hiddenzone_inputs['mstruct'])

                    # Update of sheath outputs
                    curr_organ_outputs['mstruct'] += delta_sheath_mstruct + share_leaf_getting_visible_mstruct * hiddenzone_inputs['mstruct']
                    curr_organ_outputs['Nstruct'] += delta_sheath_Nstruct + share_leaf_getting_visible_Nstruct * hiddenzone_inputs['leaf_enclosed_Nstruct']
                    curr_organ_outputs['sucrose'] += export_sucrose + share_leaf_getting_visible_mstruct * hiddenzone_inputs['sucrose']
                    curr_organ_outputs['amino_acids'] += export_amino_acids + share_leaf_getting_visible_mstruct * hiddenzone_inputs['amino_acids']
                    curr_organ_outputs['fructan'] += share_leaf_getting_visible_mstruct * hiddenzone_inputs['fructan']
                    curr_organ_outputs['proteins'] += share_leaf_getting_visible_mstruct * hiddenzone_inputs['proteins']
                    self.outputs['organs'][sheath_id] = curr_organ_outputs


            # CN consumption due to mstruct/Nstruct growth of the enclosed leaf and of the internode
            curr_hiddenzone_outputs['AA_consumption_mstruct'] = model.calculate_s_Nstruct_amino_acids( (delta_leaf_enclosed_Nstruct + delta_internode_enclosed_Nstruct), delta_lamina_Nstruct, delta_sheath_Nstruct)  #: Consumption of amino acids due to mstruct growth (µmol N)
            curr_hiddenzone_outputs['sucrose_consumption_mstruct'] = model.calculate_s_mstruct_sucrose( (delta_leaf_enclosed_mstruct + delta_internode_enclosed_mstruct), delta_lamina_mstruct, delta_sheath_mstruct, curr_hiddenzone_outputs['AA_consumption_mstruct']) #: Consumption of sucrose due to mstruct growth (µmol C)
            curr_hiddenzone_outputs['Respi_growth'] = RespirationModel.R_growth(curr_hiddenzone_outputs['sucrose_consumption_mstruct'])                                      #: Respiration growth (µmol C)

            # Update of outputs
            curr_hiddenzone_outputs['leaf_enclosed_mstruct'] += delta_leaf_enclosed_mstruct
            curr_hiddenzone_outputs['leaf_enclosed_Nstruct'] += delta_leaf_enclosed_Nstruct
            curr_hiddenzone_outputs['internode_mstruct'] += delta_internode_enclosed_mstruct
            curr_hiddenzone_outputs['internode_Nstruct'] += delta_internode_enclosed_Nstruct
            curr_hiddenzone_outputs['sucrose'] -= (curr_hiddenzone_outputs['sucrose_consumption_mstruct'] + curr_hiddenzone_outputs['Respi_growth'] + export_sucrose) # TODO: Add checks for negative sucrose value
            curr_hiddenzone_outputs['amino_acids'] -= (curr_hiddenzone_outputs['AA_consumption_mstruct'] + export_amino_acids) # TODO: Add checks for negative AA value
            self.outputs['hiddenzone'][hiddenzone_id] = curr_hiddenzone_outputs

            # Remobilisation at the end of leaf elongation
            if not hiddenzone_inputs['leaf_is_growing'] and hiddenzone_inputs['delta_leaf_L'] > 0:
               share = curr_hiddenzone_outputs['leaf_mstruct'] / curr_hiddenzone_outputs['mstruct']

               ## Add to sheath
               sheath_id = hiddenzone_id + tuple(['sheath'])
               curr_organ_outputs = self.outputs['organs'][sheath_id]
               curr_organ_outputs['mstruct'] += curr_hiddenzone_outputs['leaf_enclosed_mstruct']
               curr_organ_outputs['Nstruct'] += curr_hiddenzone_outputs['leaf_enclosed_Nstruct']
               curr_organ_outputs['sucrose'] += curr_hiddenzone_outputs['sucrose'] * share
               curr_organ_outputs['amino_acids'] += curr_hiddenzone_outputs['amino_acids'] * share
               curr_organ_outputs['fructan'] += curr_hiddenzone_outputs['fructan'] * share
               curr_organ_outputs['proteins'] += curr_hiddenzone_outputs['proteins'] * share
               self.outputs['organs'][sheath_id] = curr_organ_outputs

               ## Remove in hiddenzone
               curr_hiddenzone_outputs = self.outputs['hiddenzone'][hiddenzone_id]
               curr_hiddenzone_outputs['leaf_enclosed_mstruct'] = 0
               curr_hiddenzone_outputs['leaf_enclosed_Nstruct'] = 0
               curr_hiddenzone_outputs['sucrose'] -= curr_hiddenzone_outputs['sucrose'] * share
               curr_hiddenzone_outputs['amino_acids'] -= curr_hiddenzone_outputs['amino_acids'] * share
               curr_hiddenzone_outputs['fructan'] -= curr_hiddenzone_outputs['fructan'] * share
               curr_hiddenzone_outputs['proteins'] -= curr_hiddenzone_outputs['proteins'] * share
               self.outputs['hiddenzone'][hiddenzone_id] = curr_hiddenzone_outputs

            # Remobilisation at the end of internode elongation
            if not hiddenzone_inputs['internode_is_growing'] and hiddenzone_inputs['delta_internode_L'] > 0: # Internodes stop to elongate after leaves

               ## Add to internode
               internode_id = hiddenzone_id + tuple(['internode'])
               curr_organ_outputs = self.outputs['organs'][internode_id] # TODO: Initialise internode organ outputs when starting its elongation
               curr_organ_outputs['mstruct'] += curr_hiddenzone_outputs['internode_enclosed_mstruct']
               curr_organ_outputs['Nstruct'] += curr_hiddenzone_outputs['internode_enclosed_Nstruct']
               curr_organ_outputs['sucrose'] += curr_hiddenzone_outputs['sucrose']
               curr_organ_outputs['amino_acids'] += curr_hiddenzone_outputs['amino_acids']
               curr_organ_outputs['fructan'] += curr_hiddenzone_outputs['fructan'] * share
               curr_organ_outputs['proteins'] += curr_hiddenzone_outputs['proteins'] * share
               self.outputs['organs'][internode_id] = curr_organ_outputs


        # Roots
        for root_id, root_inputs in all_roots_inputs.iteritems():
            curr_root_outputs = all_roots_outputs[root_id]
            # Growth
            mstruct_C_growth, mstruct_growth, Nstruct_growth, Nstruct_N_growth = model.calculate_roots_mstruct_growth(root_inputs['sucrose'], root_inputs['amino_acids'], root_inputs['mstruct'], self.delta_t)
            # Respiration growth
            curr_root_outputs['Respi_growth'] = RespirationModel.R_growth(mstruct_C_growth)
            # Update of root outputs
            curr_root_outputs['mstruct'] += mstruct_growth
            curr_root_outputs['sucrose'] -= (mstruct_C_growth + curr_root_outputs['Respi_growth'])
            curr_root_outputs['Nstruct'] += Nstruct_growth
            curr_root_outputs['amino_acids'] -= (Nstruct_N_growth)
            curr_root_outputs['rate_mstruct_growth'] = mstruct_growth / self.delta_t
            self.outputs['roots'][root_id] = curr_root_outputs