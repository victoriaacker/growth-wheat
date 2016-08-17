# -*- coding: latin-1 -*-
"""
    main
    ~~~~

    An example to show how to:
        
        * initialize and run the model Growth-Wheat,
        * format the outputs of Growth-Wheat.

    You must first install :mod:`growthwheat` and its dependencies 
    before running this script with the command `python`.

    :copyright: Copyright 2014-2016 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2016.
    
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""

import os

import numpy as np
import pandas as pd

from growthwheat import model, simulation, converter

INPUTS_DIRPATH = 'inputs'
HGZ_INPUTS_FILENAME = 'hgz_inputs.csv'
ORGAN_INPUTS_FILENAME = 'organ_inputs.csv'

OUTPUTS_DIRPATH = 'outputs'
HGZ_OUTPUTS_FILENAME = 'hgz_outputs.csv'
ORGAN_OUTPUTS_FILENAME = 'organ_outputs.csv'

OUTPUTS_PRECISION = 10

if __name__ == '__main__':

    # Create population
    simulation_ = simulation.Simulation(delta_t=3600)
    # read inputs from Pandas dataframe
    hgz_inputs_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, HGZ_INPUTS_FILENAME))
    organ_inputs_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, ORGAN_INPUTS_FILENAME))
    # convert the dataframe to simulation inputs format
    inputs = converter.from_dataframes(hgz_inputs_df, organ_inputs_df)
    # initialize the simulation with the inputs
    simulation_.initialize(inputs)
    # run the simulation
    simulation_.run()
    # convert the outputs to Pandas dataframe
    hgz_outputs_df, organ_outputs_df = converter.to_dataframes(simulation_.outputs)
    # write the dataframe to CSV
    hgz_outputs_df.to_csv(os.path.join(OUTPUTS_DIRPATH, HGZ_OUTPUTS_FILENAME), index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION))
    organ_outputs_df.to_csv(os.path.join(OUTPUTS_DIRPATH, ORGAN_OUTPUTS_FILENAME), index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION))