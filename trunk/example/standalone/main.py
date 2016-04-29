# -*- coding: latin-1 -*-
"""
    test_growthwheat
    ~~~~~~~~~~~~~~~~

    Test the Growth-Wheat model.

    You must first install :mod:`growthwheat` (and add it to your PYTHONPATH)
    before running this script with the command `python`.

    To get a coverage report of :mod:`growthwheat`, use the command:
    `nosetests --with-coverage --cover-package=growthwheat test_growthwheat.py`.

    CSV files must contain only ASCII characters and ',' as separator.

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

import os

import numpy as np
import pandas as pd

from growthwheat import model, simulation, converter

INPUTS_DIRPATH = 'inputs'
HGZ_INPUTS_FILENAME = 'hgz_inputs.csv'
ELEMENT_INPUTS_FILENAME = 'exposed_element_inputs.csv'

OUTPUTS_DIRPATH = 'outputs'
HGZ_OUTPUTS_FILENAME = 'hgz_outputs.csv'
ELEMENT_OUTPUTS_FILENAME = 'exposed_element_outputs.csv'

OUTPUTS_PRECISION = 10

if __name__ == '__main__':
    start_time = 0
    stop_time = 1
    growthwheat_ts = 1
    number_of_output_steps = 2
    # Create population
    simulation_ = simulation.Simulation(delta_t=3600)
    # read inputs from Pandas dataframe
    hgz_inputs_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, HGZ_INPUTS_FILENAME))
    element_inputs_df = pd.read_csv(os.path.join(INPUTS_DIRPATH, ELEMENT_INPUTS_FILENAME))
    # convert the dataframe to simulation inputs format
    inputs = converter.from_dataframe(hgz_inputs_df, element_inputs_df)
    # initialize the simulation with the inputs
    simulation_.initialize(inputs)
    # run the simulation
    simulation_.run()
    # convert the outputs to Pandas dataframe
    hgz_outputs_df, element_outputs_df = converter.to_dataframe(simulation_.outputs)
    # write the dataframe to CSV
    hgz_outputs_df.to_csv(os.path.join(OUTPUTS_DIRPATH, HGZ_OUTPUTS_FILENAME), index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION))
    element_outputs_df.to_csv(os.path.join(OUTPUTS_DIRPATH, ELEMENT_OUTPUTS_FILENAME), index=False, na_rep='NA', float_format='%.{}f'.format(OUTPUTS_PRECISION))