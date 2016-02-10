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
from scipy.interpolate import interp1d

#### LOGGING###
import logging
import logging.config
import json
LOGGING_CONFIG_FILEPATH = os.path.join('..', 'logging.json')

LOGGING_LEVEL = logging.INFO # can be one of: DEBUG, INFO, WARNING, ERROR, CRITICAL

def setup_logging(config_filepath='logging.json', level=logging.INFO,
                  log_model=True, log_compartments=True, log_derivatives=True):
    """Setup logging configuration.
    """
    if os.path.exists(config_filepath):
        with open(config_filepath, 'r') as f:
            config = json.load(f)
        logging.config.dictConfig(config)
    else:
        logging.basicConfig()
    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    logging.getLogger('growthwheat.model').disabled = not log_model # set to False to log messages from cnwheat.model
    logging.getLogger('growthwheat.compartments').disabled = not log_compartments # set to False to log the compartments
    logging.getLogger('growthwheat.derivatives').disabled = not log_derivatives # set to False to log the derivatives

setup_logging(LOGGING_CONFIG_FILEPATH, LOGGING_LEVEL, log_model=False, log_compartments=False, log_derivatives=False)
########

from growthwheat import model
from growthwheat import simulation
from growthwheat import parameters

INPUTS_DIRPATH = 'inputs'

PHOTOSYNTHESE_LIMBE_FILENAME = 'photosynthese_limbe.csv'
PHOTOSYNTHESE_GAINE_FILENAME = 'photosynthese_gaine.csv'

OUTPUTS_DIRPATH = 'outputs'

INPUTS_OUTPUTS_PRECISION = 6

def test_run():

    # croissance
    L = parameters.Main.L0
    Le = parameters.Main.L0
    Lem = 0.0
    M_em = 0.0
    Mstruc_tot = 6.39e-8
    S_photoS = 0.0
    S_tot = 0.0
    W_photoS_bis = 0.0
    Wbis = 0.0
    x = 0.0

    # Limbe_emerge
    Camidon_em = 0.0
    Cfruc_em = 0.0
    Csuc_em = 0.0
    Ctri_em = 0.0
    Mstruc_em = 0.0
    Naa_em = 0.0
    Proteines_em = 0.0

    # phloeme
    Csuc_phlo = 583
    Naa_phlo = 61

    # Zone_cachee
    C_respi_croi = 0.0
    Camidon_croi = 0.0
    Csuc_pool_croi = 1E-3
    Ctri_croi = 0.0
    Fruc_pool_croi = 0.0
    Fruc_pool_E = 3.5E-4
    Mstruc_croi = 6.39e-8
    Mstruc_E = 1e-6
    Naa_pool_croi = 2E-4
    Prot_pool_croi = 0.0

    start_time = 0
    stop_time = 800
    growthwheat_ts = 1
    number_of_output_steps = 2

    # Photosynthese df
    ## Tableau de valeurs de photosynthèse surfacique gaine. En µmol de C/m²/s
    photosynthese_gaine_df = pd.read_csv(os.path.join('inputs', 'photosynthese_gaine.csv'), index_col='t')
    new_index = np.linspace(start_time, stop_time+1, stop_time*number_of_output_steps + 3)
    photosynthese_gaine_df = photosynthese_gaine_df.loc[:stop_time+1].reindex(new_index)
    photosynthese_gaine_df.interpolate(inplace=True)

    ## photosynthèse surfacique limbe
    photosynthese_limbe_df = pd.read_csv(os.path.join('inputs', 'photosynthese_limbe.csv'), index_col='t')
    new_index = np.linspace(start_time, stop_time+1, stop_time*number_of_output_steps + 3)
    photosynthese_limbe_df = photosynthese_limbe_df.loc[:stop_time+1].reindex(new_index)
    photosynthese_limbe_df.interpolate(inplace=True)


    # Mstruc
    Mstruc_df = pd.read_csv(os.path.join('inputs', 'Mstruc.csv'), index_col='t')

    # Create population
    simulation_ = simulation.Simulation()

    global_croissance_df = pd.DataFrame(columns=simulation_.CROISSANCE)
    global_Zone_cachee_df = pd.DataFrame(columns=simulation_.ZONE_CACHEE)
    global_Limbe_emerge_df = pd.DataFrame(columns=simulation_.LIMBE_EMERGE)
    global_phloeme_df = pd.DataFrame(columns=simulation_.PHLOEME)

    croissance_df_list = []
    Zone_cachee_df_list = []
    Limbe_emerge_df_list = []
    phloeme_df_list = []

    croissance = model.croissance(L, Le, Lem, M_em, Mstruc_tot, S_photoS, S_tot, W_photoS_bis, Wbis, x)
    Limbe_emerge = model.Limbe_emerge(Camidon_em, Cfruc_em, Csuc_em, Ctri_em, Mstruc_em, Naa_em, Proteines_em)
    Limbe_emerge.Mstruc_em_df = interp1d(range(0, 801), Mstruc_df['Mstruc_em'])
    phloeme = model.phloeme(Csuc_phlo, Naa_phlo)
    Zone_cachee = model.Zone_cachee(C_respi_croi, Camidon_croi, Csuc_pool_croi, Ctri_croi, Fruc_pool_croi, Fruc_pool_E, Mstruc_croi, Mstruc_E, Naa_pool_croi, Prot_pool_croi)
    Zone_cachee.Mstruc_croi_df = interp1d(range(0, 801), Mstruc_df['Mstruc_croi'])

    population = [croissance, Limbe_emerge, phloeme, Zone_cachee]

    simulation_.initialize(population)

    for t_growthwheat in xrange(start_time, stop_time, growthwheat_ts):
        # Photosynthesis data
        Limbe_emerge.Photosynthesis_limbe = photosynthese_limbe_df[photosynthese_limbe_df.index==t_growthwheat]['Photosynthese_nette'].iloc[0]
        Limbe_emerge.Multi_dix_limbe = photosynthese_limbe_df[photosynthese_limbe_df.index==t_growthwheat]['multipli_respi'].iloc[0]
        Zone_cachee.Multi_dix_gaine = photosynthese_gaine_df[photosynthese_gaine_df.index==t_growthwheat]['Multi_dix'].iloc[0]

        if t_growthwheat >= croissance.x_em(croissance.xb()):
            Limbe_emerge.has_emerged = True

        simulation_.run(start_time=t_growthwheat, stop_time=t_growthwheat + growthwheat_ts, number_of_output_steps=number_of_output_steps + 1)
        croissance_postprocessing_df, Zone_cachee_postprocessing_df, Limbe_emerge_postprocessing_df, phloeme_postprocessing_df = simulation_.postprocessings()

        croissance_postprocessing_df = croissance_postprocessing_df.loc[croissance_postprocessing_df.t == t_growthwheat, :].reset_index(drop=True)
        Zone_cachee_postprocessing_df = Zone_cachee_postprocessing_df.loc[Zone_cachee_postprocessing_df.t == t_growthwheat, :].reset_index(drop=True)
        Limbe_emerge_postprocessing_df = Limbe_emerge_postprocessing_df.loc[Limbe_emerge_postprocessing_df.t == t_growthwheat, :].reset_index(drop=True)
        phloeme_postprocessing_df = phloeme_postprocessing_df.loc[phloeme_postprocessing_df.t == t_growthwheat, :].reset_index(drop=True)

        croissance_df_list.append(croissance_postprocessing_df)
        Zone_cachee_df_list.append(Zone_cachee_postprocessing_df)
        Limbe_emerge_df_list.append(Limbe_emerge_postprocessing_df)
        phloeme_df_list.append(phloeme_postprocessing_df)


    global_croissance_df = pd.concat(croissance_df_list, ignore_index=True)
    global_Zone_cachee_df = pd.concat(Zone_cachee_df_list, ignore_index=True)
    global_Limbe_emerge_df = pd.concat(Limbe_emerge_df_list, ignore_index=True)
    global_phloeme_df = pd.concat(phloeme_df_list, ignore_index=True)

    global_croissance_df.to_csv(os.path.join(OUTPUTS_DIRPATH, 'croissance.csv'), na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))
    global_Zone_cachee_df.to_csv(os.path.join(OUTPUTS_DIRPATH, 'Zone_cachee.csv'), na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))
    global_Limbe_emerge_df.to_csv(os.path.join(OUTPUTS_DIRPATH, 'Limbe_emerge.csv'), na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))
    global_phloeme_df.to_csv(os.path.join(OUTPUTS_DIRPATH, 'phloeme.csv'), na_rep='NA', index=False, float_format='%.{}f'.format(INPUTS_OUTPUTS_PRECISION))

if __name__ == '__main__':
    test_run()
