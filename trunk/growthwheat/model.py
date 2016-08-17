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
    if sucrose > 0:
        delta_leaf_L = parameters.K * leaf_L * ((leaf_L - leaf_Lmax)**2)**parameters.N * delta_t
    else:
        delta_leaf_L = 0
    return delta_leaf_L

def calculate_delta_mstruct_preE(leaf_L, delta_leaf_L):
    """ Relation between leaf length and the delta of mstruct in the hidden growing zone.
    Parameters alpha_mass_growth and beta_mass_growth estimated from Williams (1975) and expressed in g of dry mass.
    Parameter RATIO_MSTRUCT_DM is then used to convert in g of structural dry mass.
    Finally, parameter RATIO_CN_MSTRUCT is used to convert in g of structural dry mass for C and N only.

    :Parameters:
        - `leaf_L` (:class:`float`) - Total leaf length (mm)
        - `delta_leaf_L` (:class:`float`) - delta of leaf length (mm)
    :Returns:
        delta mstruct (g of mstruct expressed in C and N)
    :Returns Type:
        :class:`float`
    """
    return parameters.alpha_mass_growth * parameters.beta_mass_growth * leaf_L**(parameters.beta_mass_growth-1) * delta_leaf_L * parameters.RATIO_MSTRUCT_DM * parameters.RATIO_CN_MSTRUCT

def calculate_delta_mstruct_postE(delta_emerged_area, leaf_SSLW):
    """ Relation between the delta of an emerged area (lamina or sheath) and the delta of mstruct.
    Parameter RATIO_CN_MSTRUCT is used to convert in g of structural dry mass for C and N only.

    :Parameters:
        - `delta_emerged_area` (:class:`float`) - delta of leaf surface (mm2)
        - `leaf_SSLW` (:class:`float`) - Structural Specific Leaf Weight (g mm-2)
    :Returns:
        delta mstruct (g of mstruct expressed in C and N)
    :Returns Type:
        :class:`float`
    """
    return leaf_SSLW * delta_emerged_area * parameters.RATIO_CN_MSTRUCT

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

def calculate_leaf_Lwmax(lamina_Lmax):
    """ Position of the maximal leaf width (leaf_Wmax) along the lamina.

    :Parameters:
        - `lamina_Lmax` (:class:`float`) - Maximal lamina length (mm)
    :Returns:
        position of the maximal leaf width (mm)
    :Returns Type:
        :class:`float`
    """
    return parameters.Swmax * lamina_Lmax

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

def calculate_leaf_Wlig(leaf_Wmax):
    """ Lamina width at ligule position, also used to determine sheath width.

    :Parameters:
        - `leaf_Wmax` (:class:`float`) - Maximal leaf width (mm)
    :Returns:
        leaf width (mm)
    :Returns Type:
        :class:`float`
    """
    return parameters.Klig * leaf_Wmax

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

def calculate_leaf_W(leaf_L, leaf_Lwmax, leaf_Wmax, lamina_Lmax, leaf_Wlig, leaf_Lmax):
    """ Leaf width at leaf base (from Dornbush et al., 2011).

    :Parameters:
        - `leaf_L` (:class:`float`) - Total leaf length (mm)
        - `leaf_Lwmax` (:class:`float`) - Position of the maximal leaf width (leaf_Wmax) along the leaf (mm)
        - `leaf_Wmax` (:class:`float`) - Maximal leaf width (mm)
        - `lamina_Lmax` (:class:`float`) - Maximal lamina length (mm)
        - `leaf_Wlig` (:class:`float`) - Lamina width at ligule position (mm)
        - `leaf_Lmax` (:class:`float`) - Final leaf length (mm)
    :Returns:
        leaf width (mm)
    :Returns Type:
        :class:`float`
    """
    if leaf_L >=0 and leaf_L<= leaf_Lwmax:
        leaf_W = leaf_Wmax * (leaf_L / leaf_Lwmax)**parameters.c1
    elif leaf_L > leaf_Lwmax and leaf_L <= lamina_Lmax:
        leaf_W = leaf_Wlig + ((leaf_Wmax-leaf_Wlig)/math.log(1+parameters.c2)) * math.log(1+parameters.c2*(lamina_Lmax-leaf_L)/(lamina_Lmax-leaf_Lwmax))
    elif leaf_L >= lamina_Lmax and leaf_L<=leaf_Lmax:
        leaf_W = leaf_Wlig
    else:
        leaf_W = 0.0
    return leaf_W

def calculate_s_mstruct_sucrose(delta_hgz_mstruct, delta_lamina_mstruct, delta_sheath_mstruct):
    """Consumption of sucrose for the calculated growth of mstruct (µmol C consumed by mstruct growth)

    :Parameters:
        - `delta_hgz_mstruct` (:class:`float`) - Hidden growing zone growth of mstruct (g of CN)
        - `delta_lamina_mstruct` (:class:`float`) - Lamina growth of mstruct (g of CN)
        - `delta_sheath_mstruct` (:class:`float`) - Sheath growth of mstruct (g of CN)
    :Returns:
        Sucrose consumption (µmol C)
    :Returns Type:
        :class:`float`
    """
    return (delta_hgz_mstruct + delta_lamina_mstruct + delta_sheath_mstruct) * parameters.ratio_sucrose_mstruct

def calculate_s_mstruct_amino_acids(delta_hgz_mstruct, delta_lamina_mstruct, delta_sheath_mstruct):
    """Consumption of amino acids for the calculated growth of mstruct (µmol N consumed by mstruct growth)

    :Parameters:
        - `delta_hgz_mstruct` (:class:`float`) - Hidden growing zone growth of mstruct (g of CN)
        - `delta_lamina_mstruct` (:class:`float`) - Lamina growth of mstruct (g of CN)
        - `delta_sheath_mstruct` (:class:`float`) - Sheath growth of mstruct (g of CN)
    :Returns:
        Amino acid consumption (µmol N)
    :Returns Type:
        :class:`float`
    """
    return (delta_hgz_mstruct + delta_lamina_mstruct + delta_sheath_mstruct) * parameters.ratio_amino_acids_mstruct

def Respiration(delta_mstruct):
    '''Flow from sucrose to C_respi_croi
    '''
    return (delta_mstruct / parameters.RATIO_CN_MSTRUCT)  * parameters.taux_respi * parameters.NB_C_SUCROSE

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


def calculate_lamina_W(lamina_L, leaf_Lwmax, leaf_Wmax, leaf_Wlig, lamina_Lmax):
    """ Emerged lamina width (recomputed from Dornbush et al., 2011)

    :Parameters:
        - `lamina_L` (:class:`float`) - lamina length (mm)
        - `leaf_Lwmax` (:class:`float`) - Position of the maximal leaf width (leaf_Wmax) along the lamina (mm)
        - `leaf_Wmax` (:class:`float`) - Maximal leaf width (mm)
        - `lamina_Lmax` (:class:`float`) - Maximal lamina length (mm)
        - `leaf_Wlig` (:class:`float`) - Lamina width at ligule position (mm)

    :Returns:
        lamina width (mm)
    :Returns Type:
        :class:`float`
    """
    if lamina_L <= leaf_Lwmax:
        lamina_W = leaf_Wmax * (lamina_L / leaf_Lwmax)**parameters.c1
    else: #if lamina_L > leaf_Lwmax and lamina_L<=lamina_Lmax: normalement pas besoin de ce test car quand si lamina__L = lamina_Lmax alors limbe sera mature et donc ne passerra plus ici
        lamina_W = leaf_Wlig + ((leaf_Wmax - leaf_Wlig) / math.log(1 + parameters.c2)) * math.log(1 + parameters.c2 * (lamina_Lmax - lamina_L) / (lamina_Lmax - leaf_Lwmax))
    return lamina_W

def calculate_delta_lamina_W(lamina_L, leaf_Lwmax, leaf_Wmax, delta_leaf_L, lamina_Lmax, leaf_Wlig):
    """ delta lamina width over delta_t

    :Parameters:
        - `lamina_L` (:class:`float`) - lamina length (mm)
        - `leaf_Lwmax` (:class:`float`) - Position of the maximal leaf width (leaf_Wmax) along the lamina (mm)
        - `leaf_Wmax` (:class:`float`) - Maximal leaf width (mm)
        - `delta_leaf_L` (:class:`float`) - delta of leaf length (mm)
        - `lamina_Lmax` (:class:`float`) - Maximal lamina length (mm)
        - `leaf_Wlig` (:class:`float`) - Lamina width at ligule position (mm)

    :Returns:
        delta lamina width (mm)
    :Returns Type:
        :class:`float`
    """
    if lamina_L <= leaf_Lwmax:
        delta_lamina_W = leaf_Wmax * parameters.c1 * (1 / (leaf_Lwmax**parameters.c1)) * lamina_L**(parameters.c1 - 1) * delta_leaf_L
    else: #if lamina_L > leaf_Lwmax and lamina_L<=lamina_Lmax: normalement pas besoin de ce test car quand si lamina_L = lamina_Lmax alors limbe sera mature et donc ne passerra plus ici
        delta_lamina_W = -parameters.c2 * (1 / (lamina_Lmax - leaf_Lwmax)) * ((leaf_Wmax - leaf_Wlig) / math.log(1 + parameters.c2)) * 1 / (1 + parameters.c2 * ((lamina_Lmax - lamina_L) / (lamina_Lmax - leaf_Lwmax))) * delta_leaf_L
    return delta_lamina_W


def calculate_delta_lamina_area(delta_leaf_L, lamina_W, delta_lamina_W):
    """ delta lamina area.

    :Parameters:
        - `delta_leaf_L` (:class:`float`) - delta of leaf length (mm)
        - `lamina_W` (:class:`float`) - lamina width at the previous time step(mm)
        - `delta_lamina_W` (:class:`float`) - Delta lamina width (mm)
    :Returns:
        delta lamina area (mm2)
    :Returns Type:
        :class:`float`
    """
    return delta_leaf_L * (lamina_W + delta_lamina_W / 2) #: Assumes a trapezoid form for the delta area


def calculate_export_sucrose(delta_mstruct, sucrose, hgz_mstruct):
    """Export of sucrose from the hidden growing zone towards the emerged part of the leaf integrated over delta_t (µmol C sucrose).

    :Parameters:
        - `delta_mstruct` (:class:`float`) - Delta of structural dry mass of the emerged part of the leaf (g of CN)
        - `sucrose` (:class:`float`) - Sucrose amount in the hidden growing zone (µmol C)
        - `hgz_mstruct` (:class:`float`) - Structural mass of the hidden growing zone (g of CN)

    :Returns:
        Sucrose export (µmol C)
    :Returns Type:
        :class:`float`
    """
    return delta_mstruct * max(0, (sucrose / hgz_mstruct))

def calculate_export_amino_acids(delta_mstruct, amino_acids, hgz_mstruct):
    """Export of amino acids from the hidden growing zone towards the emerged part of the leaf integrated over delta_t (µmol N amino acids).

    :Parameters:
        - `delta_mstruct` (:class:`float`) - Delta of structural dry mass of the emerged part of the leaf (g of CN)
        - `amino_acids` (:class:`float`) - Amino acids amount in the hidden growing zone (µmol N)
        - `hgz_mstruct` (:class:`float`) - Structural mass of the hidden growing zone (g of CN)

    :Returns:
        amino acids export (µmol N)
    :Returns Type:
        :class:`float`
    """
    return delta_mstruct * max(0, (amino_acids / hgz_mstruct))

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

def calculate_sheath_area(sheath_L, leaf_Wlig):
    """ Emerged sheath area.

    :Parameters:
        - `sheath_L` (:class:`float`) - Sheath length (mm)
        - `leaf_Wlig` (:class:`float`) - Lamina width at ligule position (mm)
    :Returns:
        Sheath area (mm2)
    :Returns Type:
        :class:`float`
    """
    return sheath_L * leaf_Wlig * math.pi

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
