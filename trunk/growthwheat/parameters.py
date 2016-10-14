# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    growthwheat.parameters
    ~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`growthwheat.parameters` defines the constant parameters.

    :copyright: Copyright 2015 INRA-ECOSYS, see AUTHORS.
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
C_MOLAR_MASS = 12                     #: Carbon molar mass (g mol-1)
N_MOLAR_MASS = 14                     #: Nitrogen molar mass (g mol-1)

## Shoot
ALPHA = 1.537e-02                     #: Parameter of the relation between leaf mass and leaf length (g m^(-BETA))
BETA = 1.28                           #: Parameter of the relation between leaf mass and leaf length (dimensionless)
RATIO_SUCROSE_MSTRUCT =  0.384        #: Mass of C (under carbohydrate form, g) in 1 g of mstruct (Penning de Vries, Witlage and Kremer, 1978)
RATIO_AMINO_ACIDS_MSTRUCT  = 0.0322   #: Mass of N (under amino acid/protein form, g) in 1 g of mstruct (Penning de Vries, Witlage and Kremer, 1978)
AMINO_ACIDS_C_RATIO = 3.67            #: Mean number of mol of C in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
AMINO_ACIDS_N_RATIO = 1.17            #: Mean number of mol of N in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
RATIO_MSTRUCT_DM = 0.8                #: Ratio mstruc/dry matter (dimensionless)

## Roots
VMAX_ROOTS_GROWTH = 0.015             #: Maximal rate of root structural dry matter growth (µmol C s-1 g-1 MS)
K_ROOTS_GROWTH = 1250                 #: Affinity coefficient of root structural dry matter growth (µmol C g-1 MS)
RATIO_C_MSTRUCT_ROOTS = 0.384         #: Mean contribution of carbon to root structural dry mass (g C g-1 Mstruct)
RATIO_N_MSTRUCT_ROOTS_ = 0.02         #: Mean contribution of nitrogen to root structural dry mass (g N g-1 Mstruct)
