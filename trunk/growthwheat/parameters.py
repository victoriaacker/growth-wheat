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

hour_to_second_conversion_factor = 3600. #: TODO: temporary

EC_wmax = 0.3 #: variation de + ou - 15% de maximal leaf width (SU)
Y0 = 137 #: Facteur agrandissement feuille en mode automate (SU)
K = 0.0239 #: Parameter of the growth function after previous leaf emergence
N = 0.4612 #: Parameter of the growth function after previous leaf emergence
Ksslw = 10000 #: Affinité SSLW aux fructanes
Kc = 350 #: affinité du RER au C (µmol/g)
Kn = 40 #: affinité du RER à N (µmol/g)
min_SSLW = 2.2e-05 #: g/mm²
max_SSLW = 5e-05 #: g/mm²
ratio_SSSW_SSLW = 5 # ratio gaine/limbe des matieres seches structurales spécifiques (calculé depuis les données de J. Bertheloot, 2004)
RERmax = 4e-06 #: s-1