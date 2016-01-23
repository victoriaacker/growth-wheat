# -*- coding: latin-1 -*-

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

t = 800

class Main:
    alpha_mass_growth = 1.73e-05
    beta_mass_growth = 1.74
    C_contents_roots = 0.96
    conductance = 8e-05
    conv_aa_C = 4.15
    conv_aa_N = 1.25
    conv_suc_C = 12
    correction_m = 0.416
    delta_Dproteins = 1.85e-06
    delta_Dstorage = 0.0001
    delta_t = 3600
    K_Dfructan = 100
    K_Regul_Sfructans = 0.5
    K_Saa_nitrates = 3
    K_Saa_trioses = 0.2
    K_Sfructans = 5000
    K_Sprotein = 100
    K_storage = 20
    K_sucrose = 0.66
    Kc = 350
    Kn = 40
    Ksslw = 10000
    L0 = 0.04
    Maa = 6.86e-05
    Mass_plant = 0.25
    max_SSLW = 5e-05
    MC = 1.2e-05
    min_SSLW = 2.2e-05
    MN = 1.4e-05
    Msuc = 0.000144
    N_feuille = 5
    n_Regul_Sfructans = 15
    n_Sfructans = 3
    p1 = 0
    q = 1
    Rdark = 1e-06
    RERmax = 4e-06
    respi_roots = 60
    rg_feuille = 3
    S_R_ratio = 5
    taux_respi = 890
    tE = 300
    Vmax_aa = 1
    Vmax_Dfructan = 0.035
    Vmax_Regul_Sfructans = 1
    Vmax_Sfructans = 0.2
    Vmax_Sprotein = 0.002
    Vmax_storage = 2
    Vmax_sucrose = 1

    
class croissance:
    c1 = 0.56
    c2 = 268274
    EC_wmax = 0.3
    Klig = 0.6
    n_wmax_suc = 3
    ratio_Mstruc = 0.8
    Y0 = 137


class Limbe_emerge:
    pass
    

class phloeme:
    pass


class Zone_cachee:
    pass




