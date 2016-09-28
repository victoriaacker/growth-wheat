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

ALPHA = 1.537e-02                     #: Parameter of the relation between leaf mass and leaf length (g m^(-BETA))
BETA = 1.28                           #: Parameter of the relation between leaf mass and leaf length (dimensionless)
RATIO_SUCROSE_MSTRUCT =  58656.667    #: Amount of C under carbohydrate form (µmol C) in 1 g of mstruct (of CN), calculation = ( (92.3e-2)/(12) - (7.7*4.15e-2)/(14*1.25) ) *1e6
RATIO_AMINO_ACIDS_MSTRUCT  = 5500.0   #: Amount of N under amino acid/protein form (µmol N) in 1 g of mstruct (of CN), calculation = (7.7e-2)/(14)* 1e6
RATIO_CN_MSTRUCT = 0.416              #: Contribution of C and N masses to mstruct
RATIO_MSTRUCT_DM = 0.8                #: Ratio mstruc/dry matter (dimensionless)

## Roots
VMAX_ROOTS_GROWTH = 0.015             #: Maximal rate of root structural dry matter growth (µmol C s-1 g-1 MS)
K_ROOTS_GROWTH = 1250                 #: Affinity coefficient of root structural dry matter growth (µmol C g-1 MS)
C_MOLAR_MASS = 12                     #: Carbon molar mass (g mol-1)
N_MOLAR_MASS = 14                     #: Nitrogen molar mass (g mol-1)
RATIO_C_MSTRUCT_ROOTS = 0.384         #: Mean contribution of carbon to root structural dry mass (g C g-1 Mstruct)
RATIO_N_MSTRUCT_ROOTS_ = 0.02         #: Mean contribution of nitrogen to root structural dry mass (g N g-1 Mstruct)
