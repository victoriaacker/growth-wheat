# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    growthwheat.model
    ~~~~~~~~~~~~~

    The module :mod:`growthwheat.model` defines the equations of the kinetic of leaf growth (mass flows) according to leaf elongation. Also includes root growth.

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

import math
import parameters

def calculate_delta_mstruct_preE(leaf_L, delta_leaf_L):
    """ Relation between leaf length and the delta of mstruct in the hidden growing zone.
    Parameters alpha_mass_growth and beta_mass_growth estimated from Williams (1975) and expressed in g of dry mass.
    Parameter RATIO_MSTRUCT_DM is then used to convert in g of structural dry mass.

    :Parameters:
        - `leaf_L` (:class:`float`) - Total leaf length (m)
        - `delta_leaf_L` (:class:`float`) - delta of leaf length (m)
    :Returns:
        delta mstruct (g)
    :Returns Type:
        :class:`float`
    """
    return parameters.ALPHA * parameters.BETA * leaf_L**(parameters.BETA-1) * delta_leaf_L * parameters.RATIO_MSTRUCT_DM

def calculate_delta_mstruct_postE(SSW, previous_mstruct, area):
    """ delta mstruct of emerged tissue (lamina and sheath). Function used when a model (e.g. ADEL-Wheat) has computed the plant geometry and thus updated organ area.

    :Parameters:
        - `SSW` (:class:`float`) - Structural Specific Weight (g m-2)
        - `previous_mstruct` (:class:`float`) - mstruct at the previous time step i.e. not yet updated (g)
        - `area` (:class:`float`) - Area at the current time step, as updated by the geometrical model (m2)
    :Returns:
        delta mstruct (g)
    :Returns Type:
        :class:`float`
    """
    updated_mstruct = SSW * area
    delta_mstruct = updated_mstruct - previous_mstruct
    return delta_mstruct

def calculate_export_sucrose(delta_mstruct, sucrose, hgz_mstruct):
    """Export of sucrose from the hidden growing zone towards the emerged part of the leaf integrated over delta_t (µmol C sucrose).

    :Parameters:
        - `delta_mstruct` (:class:`float`) - Delta of structural dry mass of the emerged part of the leaf (g)
        - `sucrose` (:class:`float`) - Sucrose amount in the hidden growing zone (µmol C)
        - `hgz_mstruct` (:class:`float`) - Structural mass of the hidden growing zone (g)

    :Returns:
        Sucrose export (µmol C)
    :Returns Type:
        :class:`float`
    """
    return delta_mstruct * max(0, (sucrose / hgz_mstruct))

def calculate_export_amino_acids(delta_mstruct, amino_acids, hgz_mstruct):
    """Export of amino acids from the hidden growing zone towards the emerged part of the leaf integrated over delta_t (µmol N amino acids).

    :Parameters:
        - `delta_mstruct` (:class:`float`) - Delta of structural dry mass of the emerged part of the leaf (g)
        - `amino_acids` (:class:`float`) - Amino acids amount in the hidden growing zone (µmol N)
        - `hgz_mstruct` (:class:`float`) - Structural mass of the hidden growing zone (g)

    :Returns:
        amino acids export (µmol N)
    :Returns Type:
        :class:`float`
    """
    return delta_mstruct * max(0, (amino_acids / hgz_mstruct))

def calculate_s_mstruct_sucrose(delta_hgz_mstruct, delta_lamina_mstruct, delta_sheath_mstruct):
    """Consumption of sucrose for the calculated mstruct growth (µmol C consumed by mstruct growth)

    :Parameters:
        - `delta_hgz_mstruct` (:class:`float`) - Hidden growing zone growth of mstruct (g)
        - `delta_lamina_mstruct` (:class:`float`) - Lamina growth of mstruct (g)
        - `delta_sheath_mstruct` (:class:`float`) - Sheath growth of mstruct (g)
    :Returns:
        Sucrose consumption (µmol C)
    :Returns Type:
        :class:`float`
    """
    return (delta_hgz_mstruct + delta_lamina_mstruct + delta_sheath_mstruct) * parameters.RATIO_CN_MSTRUCT * parameters.RATIO_SUCROSE_MSTRUCT

def calculate_s_mstruct_amino_acids(delta_hgz_mstruct, delta_lamina_mstruct, delta_sheath_mstruct):
    """Consumption of amino acids for the calculated mstruct growth (µmol N consumed by mstruct growth)

    :Parameters:
        - `delta_hgz_mstruct` (:class:`float`) - Hidden growing zone growth (g)
        - `delta_lamina_mstruct` (:class:`float`) - Lamina growth (g)
        - `delta_sheath_mstruct` (:class:`float`) - Sheath growth(g)
    :Returns:
        Amino acid consumption (µmol N)
    :Returns Type:
        :class:`float`
    """
    return (delta_hgz_mstruct + delta_lamina_mstruct + delta_sheath_mstruct) * parameters.RATIO_CN_MSTRUCT * parameters.RATIO_AMINO_ACIDS_MSTRUCT

## Roots
def calculate_roots_mstruct_growth(sucrose, amino_acids, mstruct, delta_t):
    """Root structural dry mass growth integrated over delta_t

    : Parameters:
        - `sucrose` (:class:`float`) - Amount of sucrose in roots (µmol C)
        - `amino_acids` (:class:`float`) - Amount of amino acids in roots (µmol N)
        - `mstruct` (:class:`float`) - Root structural mass (g)

    : Returns:
        mstruct_C_growth (µmol C), mstruct_growth (g), Nstruct_growth (g), Nstruct_N_growth (µmol N)

    :Returns Type:
        :class:`float`
    """
    conc_sucrose = max(0, sucrose/mstruct)

    mstruct_C_growth = (conc_sucrose * parameters.VMAX_ROOTS_GROWTH) / (conc_sucrose + parameters.K_ROOTS_GROWTH) * delta_t * mstruct     #: root growth in C (µmol of C)
    mstruct_growth = (mstruct_C_growth*1E-6 * parameters.C_MOLAR_MASS) / parameters.RATIO_C_MSTRUCT_ROOTS                                 #: root growth (g of structural dry mass)
    Nstruct_growth = mstruct_growth * parameters.RATIO_N_MSTRUCT_ROOTS_                                                                   #: root growth in N (g of structural dry mass)
    Nstruct_N_growth = min(amino_acids, (Nstruct_growth / parameters.N_MOLAR_MASS)*1E6)                                                   #: root growth in nitrogen (µmol N)

    return mstruct_C_growth, mstruct_growth, Nstruct_growth, Nstruct_N_growth
