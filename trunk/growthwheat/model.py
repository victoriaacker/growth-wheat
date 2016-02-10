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

import logging

class Population(object):
    """
    The class :class:`Population` defines the CN exchanges at the population scale.

    A :class:`population <Population>` must have one or several :class:`plants <Plant>`.
    """

    #PARAMETERS = parameters.PopulationParameters #: the internal parameters of the population

    def __init__(self, organs=None):
        if organs is None:
            organs = []
        self.organs = organs #: the list of organs



class Main:

    def __init__(self):
        pass

    # variables

    def Mstruc_bis(self, Mstruc_croi, Mstruc_em):
        '''
        pour vérifier que Mstruc_croi + Mstruc_em = Mstruc_tot. Ne semble pas utilise sinon
        '''
        return Mstruc_croi + Mstruc_em

    def taux_aa(self):
        '''conversion de deltaM en µmol d' aa
        '''
        return (7.7e-2)/(14*1.25)* 1e6

    def taux_suc(self):
        '''conversion de deltaM en µmol de sucrose
        '''
        return ( (92.3e-2)/(12*12) - (7.7*4.15e-2)/(14*1.25*12) ) *1e6

    def Total_unloading(self, Unloading_Csuc_phlo, Loading_Csuc_em):
        return Unloading_Csuc_phlo - Loading_Csuc_em

    # fluxes

    def Export_Csuc_em(self, Out_Csuc_em):
        '''Flow from Out_Csuc_em to In1
        '''
        return Out_Csuc_em

    def Export_Csuc_pool_croi(self, Out_Csuc_pool_croi):
        '''Flow from Out_Csuc_pool_croi to In_Csuc_em
        '''
        return Out_Csuc_pool_croi

    def Export_Mstruc_croi(self, Out_Mstruc_croi):
        '''Flow from Out_Mstruc_croi to In_Mstruc_em
        '''
        return Out_Mstruc_croi

    def Export_Naa_em(self, Out_Naa_em):
        '''Flow from Out_Naa_em to In1
        '''
        return Out_Naa_em

    def Export_Naa_pool_croi(self, Out_Naa_pool_croi):
        '''Flow from Out_Naa_pool_croi to In_Naa_em
        '''
        return Out_Naa_pool_croi

    def Loading_Csuc_em(self, Out_loading_Csuc_em):
        '''Flow from Out_loading_Csuc_em to In2
        '''
        return Out_loading_Csuc_em

    def Loading_Naa_em(self, Out_loading_Naa_em):
        '''Flow from Out_loading_Naa_em to In2
        '''
        return Out_loading_Naa_em

    def Unloading_Csuc_phlo(self, Out_Csuc_phlo):
        '''Flow from Out_Csuc_phlo to In_Csuc_pool_croi
        '''
        return Out_Csuc_phlo

    def Unloading_Naa_phlo(self, Out_Naa_phlo):
        '''Flow from Out_Naa_phlo to In_Naa_pool_croi
        '''
        return Out_Naa_phlo


class croissance:

    def __init__(self, L, Le, Lem, M_em, Mstruc_tot, S_photoS, S_tot, W_photoS_bis, Wbis, x):

        # state variables
        self.L = L
        self.Le = Le
        self.Lem = Lem
        self.M_em = M_em
        self.Mstruc_tot = Mstruc_tot
        self.S_photoS = S_photoS
        self.S_tot = S_tot
        self.W_photoS_bis = W_photoS_bis
        self.Wbis = Wbis
        self.x = x

    # variables

    def calculate_delta_x(self, t, Csuc_pool_croi, tE):
        '''temps carbone
        '''
        if Csuc_pool_croi <= 0 and tE < t:
            return 0
        return 1

    def calculate_deltaL(self, Cpool_croi, xE, x, xend, Ltot_max, Le, xmax, xb, L, Kc, Npool_croi, Kn, RERmax, delta_t):
        '''mm
        '''
        if Cpool_croi == 0:
            return 0
        if xE <= x and x <= xend and 0 < Cpool_croi:
            return (Ltot_max-Le) * ( (-1/(xend-xmax)) * ( (x-xb)/(xend-xb) ) **( (xend -xb)/(xend-xmax) ) + (1+ (xend-x)/(xend-xmax)) * (1/(xend-xmax)) * ( (x-xb)/(xend-xb)) **( (xend-xb)/(xend-xmax)-1) )
        if x < xE:
            return L * (Cpool_croi/(Kc+ Cpool_croi)) * ((Npool_croi**3)/(Kn**3+Npool_croi**3)) * RERmax*delta_t
        return 0

    def calculate_DeltaMstruc(self, x, xE, alpha_mass_growth, beta_mass_growth, L, deltaL, ratio_Mstruc, correction_m, SSLW, deltaS):
        '''g (uniquement C et N)
        '''
        if x <= xE:
            return alpha_mass_growth* beta_mass_growth * L**(beta_mass_growth-1) * deltaL * ratio_Mstruc * correction_m
        if xE < x:
            return SSLW * deltaS* correction_m
        return 0

    def calculate_deltaS(self, xE, x, xend, deltaL, W, deltaW):
        if xE <= x and x <= xend:
            return deltaL*(W + deltaW)
        return 0

    def calculate_deltaS_photoS(self, x_em, x, L, ligulation, Ltot_max, deltaL, W_photoS, deltaW_photoS):
        '''mm²
        '''
        if x_em<x and L<(ligulation * Ltot_max):
            return deltaL*(W_photoS+deltaW_photoS)
        return 0

    def calculate_deltaW(self, L, Lwmax, Wmax, c1, deltaL, L_photoS_max, c2, Wlig):
        '''mm
        '''
        if 0<=L and L<= Lwmax:
            return Wmax * c1 * (1/(Lwmax**c1)) * L**(c1-1) * deltaL
        if Lwmax<=L and L <= L_photoS_max:
            return -c2 * (1/(L_photoS_max - Lwmax) ) * ((Wmax-Wlig)/math.log(1+c2)) * 1/(1 + c2 *( (L_photoS_max - L)/(L_photoS_max-Lwmax) ) ) * deltaL
        return 0

    def calculate_deltaW_photoS(self, x_em, x, L_photoS, Lwmax, Wmax, c1, deltaL, L_photoS_max, c2, Wlig):
        if x_em<=x and 0<=L_photoS and L_photoS<=Lwmax:
            return Wmax * c1 * (1/(Lwmax**c1)) * L_photoS**(c1-1) * deltaL
        if Lwmax<L_photoS and L_photoS<L_photoS_max:
            -c2 * (1/(L_photoS_max - Lwmax) ) * ((Wmax-Wlig)/math.log(1+c2)) * 1/(1 + c2 *( (L_photoS_max - L_photoS)/(L_photoS_max-Lwmax) ) ) * deltaL
        return 0

    def calculate_L_gaine_em(self, x_em, x, L_zone_cachee, Lem):
        '''mm
        '''
        if x_em<x:
            return L_zone_cachee-Lem
        return 0

    def calculate_L_photoS(self, x, x_em, L, Ltot_max, ligulation, Lem):
        '''mm
        '''
        if x<x_em:
            return 0
        if (x_em<=x) and L< (Ltot_max*ligulation):
            return L-Lem
        if x_em<=x and (ligulation * Ltot_max)<=L:
            return ligulation *Ltot_max-Lem
        return 0

    def calculate_L_photoS_max(self, xE, x, ligulation, Ltot_max, Lem_bis):
        '''mm
        '''
        if xE<=x:
            return ligulation *Ltot_max -Lem_bis
        return 0

    def calculate_L_zone_cachee(self, L, L_photoS):
        '''mm
        '''
        return L-L_photoS

    def calculate_Lbis(self, L_zone_cachee, L_photoS):
        return L_zone_cachee + L_photoS

    def calculate_Lem_bis(self, Le, Ltot_max, xend, x_em, xmax, xE):
        return Le + (Ltot_max-Le) * (1 + (xend-x_em)/(xend-xmax) ) * ( (x_em-xE) / (xend-xE))**( (xend-xE)/(xend-xmax) )

    def calculate_ligulation(self):
        return 0.9

    def calculate_Ltot_max(self, xE, x, Le, Y0):
        '''mm
        '''
        if xE<=x:
            return Le*Y0
        return 0

    def calculate_Lwmax(self, xE, x, Swmax, L_photoS_max):
        '''mm
        '''
        if xE<=x:
            return Swmax * L_photoS_max
        return 0

    def calculate_M0(self, alpha_mass_growth, L0, beta_mass_growth):
        return alpha_mass_growth * (L0)**beta_mass_growth

    def calculate_RER(self, deltaL, L):
        return deltaL/L

    def calculate_RER_m(self, DeltaMstruc, Mstruc_tot):
        return DeltaMstruc/Mstruc_tot

    def calculate_S_gaine_em(self, x_em, x, L_gaine_em, Wlig):
        '''mm²
        '''
        if x_em<x:
            return L_gaine_em * Wlig
        return 0

    def SSLW(self, xE, x, min_SSLW, max_SSLW, Fruc_E, Ksslw):
        '''g/mm²
        '''
        if xE<x:
            return (min_SSLW + (max_SSLW - min_SSLW) * Fruc_E/ (Fruc_E + Ksslw) )
        return 0.0


    def Swmax(self, N_feuille):
        return -(N_feuille/30)+103/150

    def t_em(self, tE):
        '''h
        '''
        return tE+100

    def tb(self, tE):
        '''h
        '''
        return tE

    def tend(self, t_em):
        '''h
        '''
        return t_em+200

    def tmax(self, t_em):
        '''h
        '''
        return t_em+120

    def W(self, L, Lwmax, xE, x, Wmax, c1, L_photoS_max, Wlig, c2, Ltot_max):
        '''mm
        '''
        if 0<= L and L<= Lwmax and xE <= x:
            return Wmax *(L/Lwmax)**c1
        if Lwmax < L and L <= L_photoS_max and xE<=x:
            return Wlig + ((Wmax-Wlig)/math.log(1+c2)) * math.log(1+c2*(L_photoS_max-L)/(L_photoS_max-Lwmax))
        if L_photoS_max<=L and L<=Ltot_max:
            return Wlig
        return 0.0

    def W_photoS(self, L_photoS, Lwmax, Wmax, c1, L_photoS_max, Wlig, c2):
        if 0<L_photoS and L_photoS<=Lwmax:
            return Wmax *(L_photoS/Lwmax)**c1
        if Lwmax<L_photoS and L_photoS<=L_photoS_max:
            Wlig + ((Wmax-Wlig)/math.log(1+c2)) * math.log(1+c2*(L_photoS_max-L_photoS)/(L_photoS_max-Lwmax))
        return 0.0

    def Wlig(self, xE, x, Klig, Wmax):
        '''mm
        '''
        if xE<=x:
            return Klig * Wmax
        return 0

    def Wmax(self, xE, x, L_photoS_max, EC_wmax, Ksslw, Fruc_E):
        '''mm
        '''
        if xE<=x:
            return (0.0575 * L_photoS_max - 0.12) * (EC_wmax * 2 * Ksslw/(Ksslw + Fruc_E) + (1-EC_wmax))
        return 0.0

    def x_em(self, xb):
        return xb+100

    def xb(self):
        return 300

    def xE(self):
        return 300

    def xend(self, xb):
        return xb+300

    def xmax(self, xb):
        return xb+120

    # compartments

    def calculate_L_derivative(self, deltaL):
        '''mm
        '''
        return deltaL

    def calculate_Le_derivative(self, x, xE, deltaL):
        if x<xE:
            return deltaL
        return 0

    def calculate_Lem_derivative(self, x, x_em, deltaL):
        '''mm
        '''
        if x<=x_em:
            return deltaL
        return 0

    def calculate_M_em_derivative(self, x, x_em, DeltaMstruc):
        if x<x_em:
            return DeltaMstruc
        return 0

    def calculate_Mstruc_tot_derivative(self, DeltaMstruc):
        '''g (uniquement C et N)
        '''
        return DeltaMstruc

    def calculate_S_photoS_derivative(self, deltaS_photoS):
        '''mm²
        '''
        return deltaS_photoS

    def calculate_S_tot_derivative(self, deltaS):
        '''mm²
        '''
        return deltaS

    def calculate_W_photoS_bis_derivative(self, deltaW_photoS):
        return deltaW_photoS

    def calculate_Wbis_derivative(self, deltaW):
        return deltaW

    def calculate_x_derivative(self, delta_x):
        '''temps carbone
        '''
        return delta_x


class Limbe_emerge:

    def __init__(self, Camidon_em, Cfruc_em, Csuc_em, Ctri_em, Mstruc_em, Naa_em, Proteines_em, Photosynthesis_limbe=None, Multi_dix_limbe = None, has_emerged=False):

        self.Camidon_em = Camidon_em
        self.Cfruc_em = Cfruc_em
        self.Csuc_em = Csuc_em
        self.Ctri_em =Ctri_em
        self.Mstruc_em = Mstruc_em
        self.Naa_em = Naa_em
        self.Proteines_em = Proteines_em
        self.Photosynthesis_limbe = Photosynthesis_limbe
        self.Multi_dix_limbe = Multi_dix_limbe
        self.has_emerged = has_emerged

    # variables

    def Amidon(self, Mstruc_em, Camidon_em):
        if 0<Mstruc_em:
            return Camidon_em/Mstruc_em
        return 0

    def C_em(self, Mstruc_em, Csuc_em):
        '''µmol/g de MS émergée
        '''
        return Csuc_em/Mstruc_em

    def Fruc_em(self, Mstruc_em, Cfruc_em):
        if 0<Mstruc_em:
            return Cfruc_em/Mstruc_em
        return 0

    def Metabo_tot_em(self, Proteines_em, Naa_em, MN, conv_aa_C, conv_aa_N, MC, Csuc_em, Cfruc_em, Ctri_em, Camidon_em):
        return (Proteines_em + Naa_em) *(MN + (conv_aa_C/conv_aa_N)*MC) + (Csuc_em +  Cfruc_em + Ctri_em + Camidon_em) *MC

    def N_em(self, Mstruc_em, Naa_em):
        if 0<Mstruc_em:
            return Naa_em/Mstruc_em
        return 0

    def PhotoS(self, x_em, x, Photosynthesis, S_photoS, delta_t):
        if x_em<x:
            return Photosynthesis*1E-6*S_photoS*delta_t*0.5
        return 0

    def Prot_em(self, Mstruc_em, Proteines_em):
        if 0<Mstruc_em:
            return Proteines_em / Mstruc_em
        return 0

    def Regul_Sfructanes_limbe(self, Vmax_Regul_Sfructans, K_Regul_Sfructans, n_Regul_Sfructans, Loading_Csuc_em, Export_Csuc_em):
        return (Vmax_Regul_Sfructans*K_Regul_Sfructans**n_Regul_Sfructans)/(K_Regul_Sfructans**n_Regul_Sfructans + (max(0,Loading_Csuc_em)+max(0, Export_Csuc_em))**n_Regul_Sfructans)

    def Respi(self, x_em, x, multi_dix, S_photoS, delta_t):
        if x_em<=x:
            return max(0,multi_dix)*1E-6*S_photoS*delta_t*0.1
        return 0

    def S_amino_acids(self, Vmax_aa, Tri, K_Saa_trioses, delta_t):
        return 0.001*Vmax_aa * max(0,Tri) /(1 + K_Saa_trioses) * delta_t

    def test1(self, S_storage, S_sucrose, S_amino_acidsC_limbe, Mstruc_em):
        return (S_storage + S_sucrose+S_amino_acidsC_limbe)*Mstruc_em

    def Tri(self, Mstruc_em, Ctri_em):
        if 0<Mstruc_em:
            return (Ctri_em/Mstruc_em)
        return 0

    # fluxes

    def D_fruc_em(self, x_em, x, fruc_em, K_Dfructan, Vmax_Dfructan, C_em, delta_t):
        '''Flow from Cfruc_em to Csuc_em
        '''
        fc = (fruc_em / (fruc_em + 1))**2
        D_fruc_em = ((K_Dfructan * Vmax_Dfructan)/(max(0,C_em) + K_Dfructan)) * delta_t * fc

        return D_fruc_em

    def D_prot_em(self, Mstruc_em, delta_Dproteins, Prot_em, delta_t):
        '''Flow from Proteines_em to Naa_em
        '''
        if 0<Mstruc_em:
            return max(0,(delta_Dproteins * Prot_em)) * delta_t
        return 0

    def D_storage(self, x_em, x, delta_Dstorage, Amidon, delta_t):
        '''Flow from Camidon_em to Csuc_em
        '''
        if x_em<x:
            return delta_Dstorage *max(0,Amidon)*delta_t
        return 0

    def Export_Csuc_em(self, x, xend, C_em, Cpool_croi, q, Mstruc_em, conductance, delta_t):
        '''Flow from Csuc_em to Out_Csuc_em
        '''
        if x<xend:
            return (C_em - Cpool_croi)*(q*Mstruc_em)**(2/3)*conductance*delta_t
        return 0

    def Export_Csuc_pool_croi(self, In_Csuc_em):
        '''Flow from In_Csuc_em to Csuc_em
        '''
        return In_Csuc_em

    def Export_Mstruc_croi(self, In_Mstruc_em):
        '''Flow from In_Mstruc_em to Mstruc_em
        '''
        return In_Mstruc_em

    def Export_Naa_em(self, Mstruc_em, x, xend, N_em, Npool_croi, q, Mstruc_croi, conductance, delta_t):
        '''Flow from Naa_em to Out_Naa_em
        '''
        if 0<Mstruc_em and x<xend:
            return (N_em-Npool_croi)*(q*Mstruc_croi)**(2/3)*conductance *delta_t
        return 0

    def Export_Naa_pool_croi(self, In_Naa_em):
        '''Flow from In_Naa_em to Naa_em
        '''
        return In_Naa_em

    def Loading_Csuc_em(self, xend, x, C_em, Cphlo, q, Mstruc_em, conductance, delta_t):
        '''Flow from Csuc_em to Out_loading_Csuc_em
        '''
        if xend<=x:
            return (C_em - Cphlo)*(q*Mstruc_em)**(2/3)*conductance*delta_t
        return 0

    def Loading_Naa_em(self, xend, x, N_em, Nphlo, q, Mstruc_croi, conductance, delta_t):
        '''Flow from Naa_em to Out_loading_Naa_em
        '''
        if xend<=x:
            return (N_em-Nphlo)*(q*Mstruc_croi)**(2/3)*conductance *delta_t
        return 0

    def S_amino_acidsC_limbe(self, S_amino_acids, conv_aa_C, conv_aa_N):
        '''Flow from Ctri_em to Naa_em
        '''
        return S_amino_acids * conv_aa_C/ conv_aa_N

    def S_fruc_em(self, x_em, x, C_em, n_Sfructans, Vmax_Sfructans, K_Sfructans, delta_t, Regul_Sfructanes_limbe):
        '''Flow from Csuc_em to Cfruc_em
        '''
        if x_em<x:
            return (((max(0, C_em)) ** n_Sfructans) * Vmax_Sfructans)/ ( (max(0, C_em)) ** n_Sfructans + K_Sfructans**n_Sfructans) * delta_t *Regul_Sfructanes_limbe
        return 0

    def S_prot_em(self, Mstruc_em, Vmax_Sprotein, N_em, K_Sprotein, delta_t):
        '''Flow from Naa_em to Proteines_em
        '''
        if 0<Mstruc_em:
            return (Vmax_Sprotein*max(0,N_em))/(K_Sprotein + max(0, N_em)) * delta_t
        return 0

    def S_storage(self, x_em, x, Tri, Vmax_storage, K_storage, delta_t):
        '''Flow from Ctri_em to Camidon_em
        '''
        if x_em<x:
            return (max(0, Tri) * Vmax_storage) / (max(0,Tri)+K_storage)*delta_t
        return 0

    def S_sucrose(self, x_em, x, Tri, Vmax_sucrose, K_sucrose, delta_t):
        '''Flow from Ctri_em to Csuc_em
        '''
        if x_em<x:
            return (max(0,Tri) * Vmax_sucrose / (max(0,Tri) +K_sucrose))*delta_t
        return 0


    # compartments

    def calculate_Camidon_em_derivative(self, Mstruc_em, S_storage, D_storage):
        '''µmol de C
        '''
        if 0<Mstruc_em:
            return (+S_storage-D_storage) * Mstruc_em
        return 0

    def calculate_Cfruc_em_derivative(self, S_fruc_em, D_fruc_em, Mstruc_em):
        return (+S_fruc_em-D_fruc_em)*Mstruc_em

    def calculate_Csuc_em_derivative(self, x_em, x, Export_Csuc_pool_croi, Export_Csuc_em, S_sucrose, D_storage, Mstruc_em, Loading_Csuc_em, Respi, D_fruc_em, S_fruc_em):
        '''µmol de C
        '''
        if x_em<x:
            return +Export_Csuc_pool_croi - Export_Csuc_em  +(S_sucrose+D_storage )*Mstruc_em -Loading_Csuc_em - Respi+ Mstruc_em*(D_fruc_em-S_fruc_em)
        return 0

    def calculate_Ctri_em_derivative(self, x_em, x, PhotoS, S_sucrose, S_storage, S_amino_acidsC_limbe, Mstruc_em):
        '''µmol de C
        '''
        if x_em<x:
            return PhotoS+(-S_sucrose-S_storage -S_amino_acidsC_limbe ) * Mstruc_em
        return 0

    def calculate_Mstruc_em_derivative(self, Export_Mstruc_croi):
        '''g (uniquement C et N)
        '''
        return +Export_Mstruc_croi

    def calculate_Naa_em_derivative(self, Mstruc_em, Export_Naa_pool_croi, Export_Naa_em, Loading_Naa_em, S_prot_em, D_prot_em, S_amino_acids):
        '''µmol de N
        '''
        if 0<Mstruc_em:
            return Export_Naa_pool_croi - Export_Naa_em -Loading_Naa_em + (-S_prot_em +D_prot_em +S_amino_acids ) * Mstruc_em
        return 0

    def calculate_Proteines_em_derivative(self, S_prot_em, D_prot_em, Mstruc_em):
        '''µmol de N
        '''
        return (S_prot_em-D_prot_em) * Mstruc_em


class phloeme:

    def __init__(self, Csuc_phlo, Naa_phlo):

        self.Csuc_phlo = Csuc_phlo
        self.Naa_phlo = Naa_phlo

# variables

    def Cphlo(self, Csuc_phlo, Mass_plant):
        '''µmol/g
        '''
        return Csuc_phlo/Mass_plant

    def Nphlo(self, Naa_phlo, Mass_plant):
        '''µmol/g
        '''
        return Naa_phlo/Mass_plant

# fluxes

    def Loading_Csuc_em(self, In_loading_Csuc_em):
        '''Flow from In_loading_Csuc_em to Csuc_phlo
        '''
        return In_loading_Csuc_em

    def Loading_Naa_em(self, In_loading_Naa_em):
        '''Flow from In_loading_Naa_em to Naa_phlo
        '''
        return In_loading_Naa_em

    def Unloading_Csuc_phlo(self, x, Cphlo, Cpool_croi, q, Mstruc_croi, conductance, delta_t):
        '''Flow from Csuc_phlo to Out_Csuc_phlo
        '''
        if 0<=x:
            return (Cphlo-Cpool_croi)  *(q*Mstruc_croi)**(2/3)*conductance *delta_t
        return 0

    def Unloading_Naa_phlo(self, Nphlo, Npool_croi, q, Mstruc_croi, conductance, delta_t):
        '''Flow from Naa_phlo to Out_Naa_phlo
        '''
        return (Nphlo-Npool_croi) *(q*Mstruc_croi)**(2/3)*conductance *delta_t

    # compartments

    def calculate_Csuc_phlo_derivative(self, x_em, x, t):
        '''µmol
        '''
        if x_em<x:
            return 21*math.sin((t+18)*(2*3.14)/24)-0.5
        return 21*math.sin((t+18)*(2*3.14)/24)

    def calculate_Naa_phlo_derivative(self, t_em, t):
        '''µmol
        '''
        if t_em<t:
            return 10*math.sin((t+18)*(2*3.14)/24)-0.5
        return 10*math.sin((t+18)*(2*3.14)/24)


class Zone_cachee:

    def __init__(self, C_respi_croi, Camidon_croi, Csuc_pool_croi, Ctri_croi, Fruc_pool_croi, Fruc_pool_E, Mstruc_croi, Mstruc_E, Naa_pool_croi, Prot_pool_croi, Multi_dix_gaine=None):
        self.C_respi_croi = C_respi_croi
        self.Camidon_croi = Camidon_croi
        self.Csuc_pool_croi = Csuc_pool_croi
        self.Ctri_croi = Ctri_croi
        self.Fruc_pool_croi = Fruc_pool_croi
        self.Fruc_pool_E = Fruc_pool_E
        self.Mstruc_croi = Mstruc_croi
        self.Mstruc_E = Mstruc_E
        self.Naa_pool_croi = Naa_pool_croi
        self.Prot_pool_croi = Prot_pool_croi
        self.Multi_dix_gaine = Multi_dix_gaine

    # variables

    def Amidon_croi(self, Mstruc_em, Camidon_croi):
        if 0<Mstruc_em:
            return Camidon_croi/Mstruc_em
        return 0

    def Cpool_croi(self, Csuc_pool_croi, Mstruc_croi):
        '''µmol/gSDM
        '''
        return Csuc_pool_croi/Mstruc_croi

    def Fruc_croi(self, Fruc_pool_croi, Mstruc_croi):
        '''µmol/g
        '''
        return Fruc_pool_croi/Mstruc_croi

    def Fruc_E(self, Fruc_pool_E, Mstruc_E):
        '''µmol/g
        '''
        return Fruc_pool_E/Mstruc_E

    def Metabo_tot_croi(self, Prot_pool_croi, Naa_pool_croi, MN, conv_aa_C, conv_aa_N, MC, Fruc_pool_croi, Csuc_pool_croi, Camidon_croi, Ctri_croi):
        return (Prot_pool_croi + Naa_pool_croi )*(MN + (conv_aa_C/conv_aa_N) * MC) + (Fruc_pool_croi + Csuc_pool_croi + Camidon_croi +Ctri_croi) * MC

    def Npool_croi(self, Naa_pool_croi, Mstruc_croi):
        '''µmol/gSDM
        '''
        return Naa_pool_croi/Mstruc_croi

    def PhotoS_gaine(self, Multi_dix, S_gaine_em, delta_t):
        '''µmol de C/h
        '''
        return max(0,Multi_dix)*(1e-6)* S_gaine_em * delta_t*0.1

    def Prot_croi(self, Prot_pool_croi, Mstruc_croi):
        '''µmol/g
        '''
        return Prot_pool_croi/Mstruc_croi

    def Regul_Sfructanes_gaine(self, Vmax_Regul_Sfructans, K_Regul_Sfructans, n_Regul_Sfructans, Unloading_Csuc_phlo):
        return (Vmax_Regul_Sfructans*K_Regul_Sfructans**n_Regul_Sfructans)/(K_Regul_Sfructans**n_Regul_Sfructans + (max(0,(-Unloading_Csuc_phlo))**n_Regul_Sfructans))

    def Respi(self, x_em, x, ligulation, Ltot_max, L, Multi_dix, S_gaine_em, delta_t):
        if x_em<x and ligulation * Ltot_max<L:
            return max(0,Multi_dix) *(1e-6)* S_gaine_em * delta_t
        return 0

    def test(self, S_storage_croi, S_sucrose_croi, Mstruc_croi):
        return (S_storage_croi + S_sucrose_croi)*Mstruc_croi

    def Tri_croi(self, Mstruc_croi, Ctri_croi):
        if 0<Mstruc_croi:
            return Ctri_croi/Mstruc_croi
        return 0

    # fluxes

    def D_Fruc_pool_croi(self, Fruc_pool_croi, Mstruc_croi, K_Dfructan, Vmax_Dfructan, Cpool_croi, delta_t):
        '''Flow from Fruc_pool_croi to Csuc_pool_croi
        '''
        fc = (Fruc_pool_croi / (Fruc_pool_croi + 1))**2
        D_Fruc_pool_croi_pot = ((K_Dfructan * Vmax_Dfructan)/(Cpool_croi + K_Dfructan)) * delta_t * fc
        return D_Fruc_pool_croi_pot

    def D_Prot_pool_croi(self, delta_Dproteins, Prot_croi, delta_t):
        '''Flow from Prot_pool_croi to Naa_pool_croi
        '''
        return max(0,(delta_Dproteins * Prot_croi)) * delta_t

    def D_storage_croi(self, x_em, x, ligulation, Ltot_max, L, delta_Dstorage, Amidon_croi, delta_t):
        '''Flow from Camidon_croi to Csuc_pool_croi
        '''
        if x_em<=x and  ligulation*Ltot_max<=L:
            return delta_Dstorage *max(0,Amidon_croi)*delta_t
        return 0

    def Export_Csuc_em(self, In_Csuc_pool_croi_em):
        '''Flow from In_Csuc_pool_croi_em to Csuc_pool_croi
        '''
        return In_Csuc_pool_croi_em

    def Export_Csuc_pool_croi(self, Cpool_croi, Export_Mstruc_croi):
        '''Flow from Csuc_pool_croi to Out_Csuc_pool_croi
        '''
        return max(0,Cpool_croi) * Export_Mstruc_croi

    def Export_Mstruc_croi(self, L, ligulation, Ltot_max, correction_m, deltaS_photoS, SSLW):
        '''Flow from Mstruc_croi to Out_Mstruc_croi
        '''
        if L<= ligulation*Ltot_max:
            Export_Mstruc_croi =  correction_m*deltaS_photoS * SSLW
        else:
            Export_Mstruc_croi = 0
        return Export_Mstruc_croi

    def Export_Naa_em(self, In1):
        '''Flow from In1 to Naa_pool_croi
        '''
        return In1

    def Export_Naa_pool_croi(self, Export_Mstruc_croi, Npool_croi):
        '''Flow from Naa_pool_croi to Out_Naa_pool_croi
        '''
        return Export_Mstruc_croi * Npool_croi

    def Respiration(self, DeltaMstruc, correction_m, taux_respi, conv_suc_C):
        '''Flow from Csuc_pool_croi to C_respi_croi
        '''
        return (DeltaMstruc / correction_m)  * taux_respi *conv_suc_C

    def S_Fruc_pool_croi(self, Cpool_croi, n_Sfructans, Vmax_Sfructans, K_Sfructans, delta_t, Regul_Sfructanes_gaine):
        '''Flow from Csuc_pool_croi to Fruc_pool_croi
        '''
        return (((max(0, Cpool_croi)) ** n_Sfructans) * Vmax_Sfructans)/ ( (max(0, Cpool_croi)) ** n_Sfructans + K_Sfructans**n_Sfructans) * delta_t *Regul_Sfructanes_gaine

    def S_Mstruc_croi_Csuc(self, DeltaMstruc, taux_suc, conv_suc_C):
        '''Flow from Csuc_pool_croi to Mstruc_croi
        '''
        return DeltaMstruc * taux_suc *conv_suc_C

    def S_Mstruc_croi_Naa(self, DeltaMstruc, taux_aa, conv_aa_N):
        '''Flow from Naa_pool_croi to Mstruc_croi
        '''
        return DeltaMstruc * taux_aa *conv_aa_N

    def S_Prot_pool_croi(self, Vmax_Sprotein, Npool_croi, K_Sprotein,delta_t):
        '''Flow from Naa_pool_croi to Prot_pool_croi
        '''
        return (Vmax_Sprotein*max(0,Npool_croi))/(K_Sprotein + max(0, Npool_croi)) * delta_t

    def S_storage_croi(self, x_em, x, ligulation, Ltot_max, L, Tri_croi, Vmax_storage, K_storage, delta_t):
        '''Flow from Ctri_croi to Camidon_croi
        '''
        if x_em<=x and ligulation*Ltot_max<=L:
            return (max(0, Tri_croi) * Vmax_storage) / (max(0,Tri_croi)+K_storage)*delta_t
        return 0

    def S_sucrose_croi(self, x_em, x, ligulation, Ltot_max, L, Tri_croi, Vmax_sucrose, K_sucrose, delta_t):
        '''Flow from Ctri_croi to Csuc_pool_croi
        '''
        if x_em<=x and ligulation*Ltot_max<=L:
            return (max(0,Tri_croi) * Vmax_sucrose / (max(0,Tri_croi) +K_sucrose))*delta_t
        return 0

    def Unloading_Csuc_phlo(self, In_Csuc_pool_croi):
        '''Flow from In_Csuc_pool_croi to Csuc_pool_croi
        '''
        return In_Csuc_pool_croi

    def Unloading_Naa_phlo(self, In_Naa_pool_croi):
        '''Flow from In_Naa_pool_croi to Naa_pool_croi
        '''
        return In_Naa_pool_croi


    # compartments

    def calculate_C_respi_croi_derivative(self, Respiration, Msuc):
        '''g
        '''
        return +Respiration * Msuc

    def calculate_Camidon_croi_derivative(self, x_em, x, ligulation, Ltot_max, L, S_storage_croi, D_storage_croi, Mstruc_croi):
        '''µmol de C
        '''
        if x_em<=x and ligulation*Ltot_max<=L:
            return (+S_storage_croi-D_storage_croi) * Mstruc_croi
        return 0

    def calculate_Csuc_pool_croi_derivative(self, Unloading_Csuc_phlo, S_Mstruc_croi_Csuc, Respiration, S_Fruc_pool_croi, Mstruc_croi, D_Fruc_pool_croi, Export_Csuc_pool_croi, Export_Csuc_em, S_sucrose_croi, D_storage_croi, Respi):
        '''µmol
        '''
        Csuc_pool_croi_derivative = Unloading_Csuc_phlo  - S_Mstruc_croi_Csuc - Respiration - S_Fruc_pool_croi*Mstruc_croi + D_Fruc_pool_croi*Mstruc_croi -Export_Csuc_pool_croi+ Export_Csuc_em +(S_sucrose_croi+D_storage_croi) * Mstruc_croi -Respi
        return Csuc_pool_croi_derivative

    def calculate_Ctri_croi_derivative(self, x_em, x, ligulation, Ltot_max, L, PhotoS_gaine, S_sucrose_croi, S_storage_croi, Mstruc_croi):
        '''µmol de C
        '''
        if x_em<=x and  ligulation*Ltot_max<=L:
            return PhotoS_gaine+(-S_sucrose_croi -S_storage_croi  ) * Mstruc_croi
        return 0

    def calculate_Fruc_pool_croi_derivative(self, S_Fruc_pool_croi, Mstruc_croi, D_Fruc_pool_croi):
        '''µmol
        '''
        return S_Fruc_pool_croi * Mstruc_croi-D_Fruc_pool_croi*Mstruc_croi

    def calculate_Fruc_pool_E_derivative(self, x, xE, S_Fruc_pool_croi, Mstruc_croi, D_Fruc_pool_croi):
        '''µmol
        '''
        if x<=xE:
            return S_Fruc_pool_croi * Mstruc_croi-D_Fruc_pool_croi*Mstruc_croi
        return 0

    def calculate_Mstruc_croi_derivative(self, x_em, x, S_Mstruc_croi_Naa, MN, S_Mstruc_croi_Csuc, MC, conv_aa_C, conv_aa_N, Export_Mstruc_croi):
        '''g uniquement C et N
        '''
        if x_em<=x:
            return (+S_Mstruc_croi_Naa * MN +S_Mstruc_croi_Csuc * MC + S_Mstruc_croi_Naa * MC *(conv_aa_C/conv_aa_N))-Export_Mstruc_croi
        if x<x_em:
            return (+S_Mstruc_croi_Naa * MN +S_Mstruc_croi_Csuc * MC + S_Mstruc_croi_Naa * MC *(conv_aa_C/conv_aa_N))

    def calculate_Mstruc_E_derivative(self, x, xE, S_Mstruc_croi_Naa, MN, S_Mstruc_croi_Csuc, MC, conv_aa_C, conv_aa_N):
        if x<=xE:
            return (+S_Mstruc_croi_Naa * MN +S_Mstruc_croi_Csuc * MC + S_Mstruc_croi_Naa * MC *conv_aa_C/conv_aa_N)
        return 0

    def calculate_Naa_pool_croi_derivative(self, Unloading_Naa_phlo, S_Mstruc_croi_Naa, S_Prot_pool_croi, Mstruc_croi, D_Prot_pool_croi, Export_Naa_pool_croi, Export_Naa_em):
        '''µmol
        '''
        return Unloading_Naa_phlo- S_Mstruc_croi_Naa- S_Prot_pool_croi * Mstruc_croi+D_Prot_pool_croi*Mstruc_croi - Export_Naa_pool_croi + Export_Naa_em

    def calculate_Prot_pool_croi_derivative(self, S_Prot_pool_croi, Mstruc_croi, D_Prot_pool_croi):
        '''µmol
        '''
        return +S_Prot_pool_croi * Mstruc_croi-D_Prot_pool_croi *Mstruc_croi
