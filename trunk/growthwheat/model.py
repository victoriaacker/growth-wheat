# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    growthwheat.model
    ~~~~~~~~~~~~~

    The module :mod:`growthwheat.model` defines the equations of the CN exchanges in a population of plants.

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

def calculate_hgz_length(prev_growing_sheath_L, prev_mature_sheath_L):
    """ length of the hidden growing zone given by the previous sheaths.

    :Parameters:
        - `prev_growing_sheath_L` (:class:`float`) - Length of the previous growing sheath (mm)
        - `prev_mature_sheath_L` (:class:`float`) - Length of the previous mature sheath (mm)
    :Returns:
        Hidden growing zone length (mm)
    :Returns Type:
        :class:`float`
    """
    return prev_growing_sheath_L + prev_mature_sheath_L

def calculate_deltaL_preE(sucrose, leaf_L, amino_acids, mstruct, delta_t):
    """ delta of leaf length over delta_t as a function of sucrose and amino acids, from initiation to the emergence of the previous leaf.

    :Parameters:
        - `sucrose` (:class:`float`) - Amount of sucrose (µmol C)
        - `leaf_L` (:class:`float`) - Total leaf length (mm)
        - `amino_acids` (:class:`float`) - Amount of amino acids (µmol N)
        - `mstruct` (:class:`float`) - Structural mass (g of C and N)
    :Returns:
        delta delta_leaf_L (mm)
    :Returns Type:
        :class:`float`
    """
    if sucrose > 0:
        delta_leaf_L = leaf_L * ((sucrose / mstruct) / (parameters.Kc + (sucrose / mstruct))) * (((amino_acids/mstruct) **3) / (parameters.Kn**3 + (amino_acids / mstruct)**3)) * parameters.RERmax * delta_t
    else:
        delta_leaf_L = 0
    return delta_leaf_L

def calculate_deltaL_postE(leaf_L, leaf_Lmax, sucrose, delta_t):
    """ delta of leaf length, from the emergence of the previous leaf to the end of growth (predefined growth kinetic depending on leaf state).

    :Parameters:
        - `leaf_L` (:class:`float`) - Total leaf length (mm)
        - `leaf_Lmax` (:class:`float`) - Final leaf length (mm)
        - `sucrose` (:class:`float`) - Amount of sucrose (µmol C)
    :Returns:
        delta delta_leaf_L (mm)
    :Returns Type:
        :class:`float`
    """
    # TODO: probleme unite au moins du temps car fonction calée sur des DD. Mais il y a autre chose
    if sucrose > 0:
        delta_leaf_L = parameters.K * leaf_L * ((leaf_L - leaf_Lmax)**2)**parameters.N * delta_t
    else:
        delta_leaf_L = 0
    return delta_leaf_L

def calculate_leaf_Lmax(leaf_Lem_prev):
    """ Final leaf length.

    :Parameters:
        - `leaf_Lem_prev` (:class:`float`) - Leaf length at the emergence of the previous leaf (mm)
    :Returns:
        Final leaf length (mm)
    :Returns Type:
        :class:`float`
    """
    return leaf_Lem_prev * parameters.Y0

def calculate_SL_ratio(metamer):
    """ Sheath:Lamina final length ratio according to the rank. Parameters from Dornbush (2011).

    :Parameters:
        - `metamer` (:class:`float`) - metamer rank
    :Returns:
        Sheath:Lamina ratio (dimensionless)
    :Returns Type:
        :class:`float`
    """
    return -0.0021 * metamer**3 + 0.037 * metamer**2 - 0.1527 * metamer + 0.4962

def calculate_lamina_Lmax(leaf_Lmax, sheath_lamina_ratio):
    """ Final lamina length.

    :Parameters:
        - `leaf_Lmax` (:class:`float`) - Final leaf length (mm)
        - `sheath_lamina_ratio` (:class:`float`) - Sheath:Lamina ratio (dimensionless)

    :Returns:
        final lamina length (mm)
    :Returns Type:
        :class:`float`
    """
    return leaf_Lmax / (1 + sheath_lamina_ratio)

def calculate_sheath_Lmax(leaf_Lmax, lamina_Lmax):
    """ Final sheath length.

    :Parameters:
        - `leaf_Lmax` (:class:`float`) - Final leaf length (mm)
        - `lamina_Lmax` (:class:`float`) - Final lamina length (mm)

    :Returns:
        final sheath length (mm)
    :Returns Type:
        :class:`float`
    """
    return leaf_Lmax - lamina_Lmax

def calculate_leaf_Wmax(lamina_Lmax, fructan, mstruct):
    """ Maximal leaf width.
    0.0575 et 0.12 issu graph Dornbush

    :Parameters:
        - `lamina_Lmax` (:class:`float`) - Maximal lamina length (mm)
        - `fructan` (:class:`float`) - Fructan in the hidden growing zone at the time of the previous leaf emergence (µmol C).
        - `mstruct` (:class:`float`) - Mstruct of the hidden growing zone at the time of the previous leaf emergence (g).
    :Returns:
        maximal leaf width (mm)
    :Returns Type:
        :class:`float`
    """
    return (0.0575 * lamina_Lmax - 0.12) * (parameters.EC_wmax * 2 * parameters.Ksslw/(parameters.Ksslw + (fructan / mstruct)) + (1-parameters.EC_wmax))

def calculate_SSLW(fructan, mstruct):
    """ Structural Specific Lamina Weight.

    :Parameters:
        - `fructan` (:class:`float`) - Fructan in the hidden growing zone at the time of the previous leaf emergence (µmol C).
        - `mstruct` (:class:`float`) - Mstruct of the hidden growing zone at the time of the previous leaf emergence (g).
    :Returns:
        Structural Specific Leaf Weight (g mm-2)
    :Returns Type:
        :class:`float`
    """
    conc_fructan = fructan / mstruct
    return parameters.min_SSLW + (parameters.max_SSLW - parameters.min_SSLW) * conc_fructan/ (conc_fructan + parameters.Ksslw)

def calculate_SSSW(SSLW):
    """ Structural Specific Sheath Weight.

    :Parameters:
        - `SSLW` (:class:`float`) - Structural Specific Leaf Weight (g mm-2).
    :Returns:
        Structural Specific Sheath Weight (g mm-2)
    :Returns Type:
        :class:`float`
    """
    return SSLW * parameters.ratio_SSSW_SSLW

def calculate_leaf_emergence(leaf_L, hgz_L):
    """Calculate if a given leaf has emerged from the hidden growing zone

    :Parameters:
        - `leaf_L` (:class:`float`) - Total leaf length (mm)
        - `hgz_L` (:class:`float`) - Length of the hidden growing zone (mm)
    :Returns:
        Specifies if the leaf has emerged (True) or not (False)
    :Returns Type:
        :class:`bool`
    """
    return leaf_L > hgz_L

def calculate_lamina_L(leaf_L, hgz_L):
    """ Emerged lamina length given by the difference between leaf length and hidden growing zone length.

    :Parameters:
        - `leaf_L` (:class:`float`) - Total leaf length (mm)
        - `hgz_L` (:class:`float`) - Length of the hidden growing zone (mm)
    :Returns:
        lamina length (mm)
    :Returns Type:
        :class:`float`
    """
    lamina_L = leaf_L - hgz_L
    if lamina_L <=0:
        print "Warning: the leaf is shorther than the hgz"
    return max(0, lamina_L)

def calculate_sheath_L(leaf_L, hgz_L, lamina_L):
    """ Emerged sheath length. Assumes that leaf_L = hgz_L + sheath_L + lamina_L

    :Parameters:
        - `leaf_L` (:class:`float`) - Total leaf length (mm)
        - `hgz_L` (:class:`float`) - Length of the hidden growing zone (mm)
        - `lamina_L` (:class:`float`) - Lamina length (mm)
    :Returns:
        sheath length (mm)
    :Returns Type:
        :class:`float`
    """
    return leaf_L - hgz_L - lamina_L

def calculate_t_prev_leaf_emerged(t_prev_leaf_emerged, delta_t):
    """Increment the time spent after the emergence of the previous leaf. If sucrose is null, then 't_prev_leaf_emerged' is not incremented

    :Parameters:
        - `t_prev_leaf_emerged` (:class:`float`) - Time spent after the emergence of the previous leaf at (t-1) (hour)
    :Returns:
         Time spent after the emergence of the previous leaf at (t) (hour)
    :Returns Type:
        :class:`float`
    """
    return t_prev_leaf_emerged + (delta_t / parameters.hour_to_second_conversion_factor)