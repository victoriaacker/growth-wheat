# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    growthwheat.simulation
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`growthwheat.simulation` is the front-end to run the CN-Wheat model.

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

import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.integrate import ode

import model, parameters # TODO: parameters temporaire
import logging

class SimulationError(Exception): pass
class SimulationRunError(SimulationError): pass

class Simulation(object):
    """
    The Simulation class permits to initialize and run the model.

    Use :meth:`run` to run the model.

    :Parameters:

        - `population` (:class:`list`) - the :class:`population <cnwheat.model.Population` in the system.

    """

    #: the name of the compartments attributes in the model.
    MODEL_COMPARTMENTS_NAMES = {model.Limbe_emerge: ['Camidon_em', 'Cfruc_em', 'Csuc_em', 'Ctri_em', 'Naa_em', 'Proteines_em'],
                                model.phloeme: ['Csuc_phlo', 'Naa_phlo'],
                                model.Zone_cachee: ['C_respi_croi', 'C_respi_croi', 'Camidon_croi', 'Csuc_pool_croi', 'Ctri_croi', 'Fruc_pool_croi', 'Fruc_pool_E', 'Mstruc_E', 'Naa_pool_croi', 'Prot_pool_croi'],
                                model.croissance: ['L', 'Le', 'Lem', 'M_em', 'Mstruc_tot', 'S_photoS', 'S_tot', 'W_photoS_bis', 'Wbis', 'x']}


    T_INDEX = 't'

    LIMBE_EMERGE = [T_INDEX, 'Camidon_em', 'Cfruc_em', 'Csuc_em', 'Ctri_em', 'Mstruc_em', 'Naa_em', 'Proteines_em']
    PHLOEME = [T_INDEX, 'Csuc_phlo', 'Naa_phlo']
    ZONE_CACHEE =  [T_INDEX, 'C_respi_croi', 'Camidon_croi', 'Csuc_pool_croi', 'Ctri_croi', 'Fruc_pool_croi', 'Fruc_pool_E', 'Mstruc_croi', 'Mstruc_E', 'Naa_pool_croi', 'Prot_pool_croi', 'D_Fruc_pool_croi']
    CROISSANCE =  [T_INDEX, 'L', 'Le', 'Lem', 'M_em', 'Mstruc_tot', 'S_photoS', 'S_tot', 'W_photoS_bis', 'Wbis', 'x']


    LOGGERS_NAMES = {'compartments': {model.Limbe_emerge: 'growthwheat.compartments.Limbe_emerge',
                                      model.phloeme: 'growthwheat.compartments.phloeme',
                                      model.Zone_cachee: 'growthwheat.compartments.Zone_cachee',
                                      model.croissance: 'growthwheat.compartments.croissance'},
                     'derivatives': {model.Limbe_emerge: 'growthwheat.derivatives.Limbe_emerge',
                                     model.phloeme: 'growthwheat.derivatives.phloeme',
                                     model.Zone_cachee: 'growthwheat.derivatives.Zone_cachee',
                                     model.croissance: 'growthwheat.derivatives.croissance'}}


    def __init__(self):

        self.population = model.Population() #: the population to simulate on

        self.initial_conditions = [] #: the initial conditions of the compartments in the population
        self.initial_conditions_mapping = {} #: dictionary to map the compartments to their indexes in :attr:`initial_conditions`

        self._time_grid = np.array([]) #: the time grid of the simulation
        self._solver_output = np.array([]) #: the value of the compartments for each time step, with the initial conditions in the first row

    def initialize(self, population):
        """
        Initialize:

            * :attr:`population`,
            * :attr:`initial_conditions_mapping`,
            * and :attr:`initial_conditions`

        from `population`.

        :Parameters:

            - `population` (:class:`model.Population`) - a population of plants.

        """

        logger = logging.getLogger(__name__)

        logger.info('Initialization of the simulation...')


        # clean the attributes of the simulation
        del self.population.organs[:]
        del self.initial_conditions[:]
        self.initial_conditions_mapping.clear()

        self.population.organs.extend(population)

        # initialize initial conditions
        def _init_initial_conditions(model_object, i):
            class_ = model_object.__class__
            if issubclass(class_, model.Limbe_emerge):
                class_ = model.Limbe_emerge
            elif issubclass(class_, model.phloeme):
                class_ = model.phloeme
            elif issubclass(class_, model.Zone_cachee):
                class_ = model.Zone_cachee
            elif issubclass(class_, model.croissance):
                class_ = model.croissance
            compartments_names = Simulation.MODEL_COMPARTMENTS_NAMES[class_]
            self.initial_conditions_mapping[model_object] = {}
            for compartment_name in compartments_names:
                if hasattr(model_object, compartment_name):
                    self.initial_conditions_mapping[model_object][compartment_name] = i
                    self.initial_conditions.append(0)
                    i += 1
            return i

        i = 0

        for organ in self.population.organs:
            if organ is None:
                continue
            i = _init_initial_conditions(organ, i)

        logger.info('Initialization of the simulation DONE')

    def run(self, start_time, stop_time, number_of_output_steps, odeint_mxstep=0, show_progressbar=False):
        """
        Compute CN exchanges in :attr:`population` from `start_time` to `stop_time`, for `number_of_output_steps` steps.

        :Parameters:

            - `start_time` (:class:`int`) - The starting of the time grid.

            - `stop_time` (:class:`int`) - The end of the time grid.

            - `number_of_output_steps` (:class:`int`) - Number of time points for which to compute the CN exchanges in :attr:`population`.

            - `odeint_mxstep` (:class:`int`) - Maximum number of (internally defined) steps allowed for each integration point in time grid.
              `odeint_mxstep` is passed to :func:`scipy.integrate.odeint` as `mxstep`. If `odeint_mxstep` = 0 (the default), then `mxstep` is determined by the solver.
              Normally, the `mxstep` determined by the solver permits to solve the current model. User can try to increase this value if a more complex model is defined
              and if the integration failed. However, take care that the origin of an integration failure could be a discontinuity in the RHS function used
              by :func:`scipy.integrate.odeint`, and that this discontinuity could be due to a bug in your model. To summary: if the integration failed, first
              check the logs.

            - `show_progressbar` (:class:`bool`) - True: show the progress bar ; False: do not show the progress bar.

        :Returns:
            Dictionary containing output information from the solver.
            This is the dictionary returned by :func:`scipy.integrate.odeint` as second output.
            See the documentation of :func:`scipy.integrate.odeint` for more information.

        :Returns Type:
            :class:`dict`

        """
        logger = logging.getLogger(__name__)
        logger.info('Run of Growth-Wheat from {} to {}...'.format(start_time, stop_time))

        t = np.linspace(start_time, stop_time, number_of_output_steps)

        compartments_logger = logging.getLogger('growthwheat.compartments')
        derivatives_logger = logging.getLogger('growthwheat.derivatives')
        if compartments_logger.isEnabledFor(logging.DEBUG) or derivatives_logger.isEnabledFor(logging.DEBUG):
            sep = ','
            if compartments_logger.isEnabledFor(logging.DEBUG):
                Limbe_emerge_compartments_logger = logging.getLogger('growthwheat.compartments.Limbe_emerge')
                Limbe_emerge_compartments_logger.debug(sep.join([Simulation.T_INDEX] + Simulation.MODEL_COMPARTMENTS_NAMES[model.Limbe_emerge]))
                Zone_cachee_compartments_logger = logging.getLogger('growthwheat.compartments.Zone_cachee')
                Zone_cachee_compartments_logger.debug(sep.join([Simulation.T_INDEX] + Simulation.MODEL_COMPARTMENTS_NAMES[model.Zone_cachee]))
                croissance_compartments_logger = logging.getLogger('growthwheat.compartments.croissance')
                croissance_compartments_logger.debug(sep.join([Simulation.T_INDEX] + Simulation.MODEL_COMPARTMENTS_NAMES[model.croissance]))
                phloeme_compartments_logger = logging.getLogger('growthwheat.compartments.phloeme')
                phloeme_compartments_logger.debug(sep.join([Simulation.T_INDEX] + Simulation.MODEL_COMPARTMENTS_NAMES[model.phloeme]))
            if derivatives_logger.isEnabledFor(logging.DEBUG):
                Limbe_emerge_derivatives_logger = logging.getLogger('growthwheat.derivatives.Limbe_emerge')
                Limbe_emerge_derivatives_logger.debug(sep.join([Simulation.T_INDEX] + Simulation.MODEL_COMPARTMENTS_NAMES[model.Limbe_emerge]))
                Zone_cachee_derivatives_logger = logging.getLogger('growthwheat.derivatives.Zone_cachee')
                Zone_cachee_derivatives_logger.debug(sep.join([Simulation.T_INDEX] + Simulation.MODEL_COMPARTMENTS_NAMES[model.Zone_cachee]))
                croissance_derivatives_logger = logging.getLogger('growthwheat.derivatives.croissance')
                croissance_derivatives_logger.debug(sep.join([Simulation.T_INDEX] + Simulation.MODEL_COMPARTMENTS_NAMES[model.croissance]))
                phloeme_derivatives_logger = logging.getLogger('growthwheat.derivatives.phloeme')
                phloeme_derivatives_logger.debug(sep.join([Simulation.T_INDEX] + Simulation.MODEL_COMPARTMENTS_NAMES[model.phloeme]))

        self._update_initial_conditions()

        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(
                """Run the solver with:
                    - time grid = %s,
                    - odeint mxstep = %s""",
                t, odeint_mxstep)

        soln, infodict = odeint(self._calculate_all_derivatives, self.initial_conditions, t, full_output=True, mxstep=odeint_mxstep)

##        r = ode(self._calculate_all_derivatives).set_integrator('dop853')
##        r.set_initial_value(self.initial_conditions, start_time)
##        while r.successful() and r.t < stop_time:
##            soln = r.integrate(r.t + 1)






        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(
                """Run of the solver DONE: infodict = %s""",
                infodict)

        # update self._time_grid
        self._time_grid.resize(t.shape)
        np.copyto(self._time_grid, t)

        # update self._solver_output
        self._solver_output.resize(soln.shape)
        np.copyto(self._solver_output, soln)

        if not set(infodict['mused']).issubset([1,2]):
            message = "Integration failed. See the logs of lsoda or try to increase the value of 'mxstep'."
            raise SimulationRunError(message)

        last_compartments_values = self._solver_output[-1]
        self._update_population(last_compartments_values)

        logger.info('Run of Growth-Wheat from {} to {} DONE'.format(start_time, stop_time))

        #return infodict


    def _update_initial_conditions(self):
        """Update the compartments values in :attr:`initial_conditions` from the compartments values of :attr:`population`.
        """
        for model_object, compartments in self.initial_conditions_mapping.iteritems():
            for compartment_name, compartment_index in compartments.iteritems():
                self.initial_conditions[compartment_index] = getattr(model_object, compartment_name)

    def _log_compartments(self, t, y, loggers_names):
        """Log the values in `y` to the loggers in `loggers_names`.
        """
        def update_rows(model_object, indexes, rows, i):
            row = []
            class_ = model_object.__class__
            if issubclass(class_, model.Limbe_emerge):
                class_ = model.Limbe_emerge
            elif issubclass(class_, model.phloeme):
                class_ = model.phloeme
            elif issubclass(class_, model.Zone_cachee):
                class_ = model.Zone_cachee
            elif issubclass(class_, model.croissance):
                class_ = model.croissance
            compartments_names = Simulation.MODEL_COMPARTMENTS_NAMES[class_]
            for compartment_name in compartments_names:
                if hasattr(model_object, compartment_name):
                    row.append(str(y[i]))
                    i += 1
                else:
                    row.append('NA')
            rows.append([str(index) for index in indexes] + row)
            return i

        i = 0
        all_rows = dict([(class_, []) for class_ in loggers_names])

        for organ in self.population.organs:
            class_ = organ.__class__
            if issubclass(class_, model.Limbe_emerge):
                class_ = model.Limbe_emerge
            elif issubclass(class_, model.phloeme):
                class_ = model.phloeme
            elif issubclass(class_, model.Zone_cachee):
                class_ = model.Zone_cachee
            elif issubclass(class_, model.croissance):
                class_ = model.croissance

            i = update_rows(organ, [t], all_rows[class_], i)

        row_sep = '\n'
        column_sep = ','
        for class_, logger_name in loggers_names.iteritems():
            compartments_logger = logging.getLogger(logger_name)
            formatted_initial_conditions = row_sep.join([column_sep.join(row) for row in all_rows[class_]])
            compartments_logger.debug(formatted_initial_conditions)

    def _calculate_all_derivatives(self, y, t):
    #def _calculate_all_derivatives(self, y, t):
        """Compute the derivative of `y` at `t`.

        :meth:`_calculate_all_derivatives` is passed as **func** argument to
        :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`.
        :meth:`_calculate_all_derivatives` is called automatically by
        :func:`scipy.integrate.odeint <scipy.integrate.odeint>`.

        First call to :meth:`_calculate_all_derivatives` uses `y` = **y0** and
        `t` = **t** [0], where **y0** and **t** are arguments passed to :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`.

        Following calls to :meth:`_calculate_all_derivatives` use `t` in [min( **t** ), max( **t** ) + 1] where
        **t** is an argument passed to :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`. `y` is
        computed automatically by the solver.

        :Parameters:

            - `y` (:class:`list`) - The current y values. `y` is automatically set by
              :func:`scipy.integrate.odeint`. User does not have control over `y`.
              At first call to :meth:`_calculate_all_derivatives` by :func:`scipy.integrate.odeint`, `y` = **y0**
              where **y0** is one of the arguments passed to :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`.
              For each following call to :meth:`_calculate_all_derivatives`, `y` is
              computed automatically by the solver.

            - `t` (:class:`float`) - The current t at which we want to compute the derivatives.
              `t` is automatically set by :func:`scipy.integrate.odeint`.
              User does not have control over `t`.
              At first call to :meth:`_calculate_all_derivatives` :func:`scipy.integrate.odeint`,
              `t` = **t** [0], where **t** is one of the arguments passed to :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`.
              For each following call to :meth:`_calculate_all_derivatives`, `t` belongs
              to the interval [min( **t** ), max( **t** ) + 1], where **t** is an
              argument passed to :func:`odeint(func, y0, t, args=(),...) <scipy.integrate.odeint>`.

        :Returns:
            The derivatives of `y` at `t`.

        :Returns Type:
            :class:`list`


        """

        logger = logging.getLogger(__name__)

        logger.debug('t = {}'.format(t))

        compartments_logger = logging.getLogger('growthwheat.compartments')
        if compartments_logger.isEnabledFor(logging.DEBUG):
            self._log_compartments(t, y, Simulation.LOGGERS_NAMES['compartments'])

        # check that the solver is not crashed
        y_isnan = np.isnan(y)
        if y_isnan.any():
            message = 'The solver did not manage to compute a compartment. See the logs.'
            logger.exception(message)
            raise SimulationRunError(message)

        for organ in self.population.organs:
            if issubclass(organ.__class__, model.Limbe_emerge):
                Limbe_emerge = organ
            elif issubclass(organ.__class__, model.phloeme):
                phloeme = organ
            elif issubclass(organ.__class__, model.Zone_cachee):
                Zone_cachee = organ
            elif issubclass(organ.__class__, model.croissance):
                croissance = organ
        Main = model.Main()

        y_derivatives = np.zeros_like(y)
        delta_t = parameters.Main.delta_t

        # Initial values

        ## croissance
        croissance.L = y[self.initial_conditions_mapping[croissance]['L']]
        croissance.Le = y[self.initial_conditions_mapping[croissance]['Le']]
        croissance.Lem = y[self.initial_conditions_mapping[croissance]['Lem']]
        croissance.M_em = y[self.initial_conditions_mapping[croissance]['M_em']]
        croissance.Mstruc_tot = y[self.initial_conditions_mapping[croissance]['Mstruc_tot']]
        croissance.S_photoS = y[self.initial_conditions_mapping[croissance]['S_photoS']]
        croissance.S_tot = y[self.initial_conditions_mapping[croissance]['S_tot']]
        croissance.W_photoS_bis = y[self.initial_conditions_mapping[croissance]['W_photoS_bis']]
        croissance.Wbis = y[self.initial_conditions_mapping[croissance]['Wbis']]
        croissance.x = y[self.initial_conditions_mapping[croissance]['x']]
        ## Limbe emerge
        Limbe_emerge.Camidon_em = y[self.initial_conditions_mapping[Limbe_emerge]['Camidon_em']]
        Limbe_emerge.Cfruc_em = y[self.initial_conditions_mapping[Limbe_emerge]['Cfruc_em']]
        Limbe_emerge.Csuc_em = y[self.initial_conditions_mapping[Limbe_emerge]['Csuc_em']]
        Limbe_emerge.Ctri_em = y[self.initial_conditions_mapping[Limbe_emerge]['Ctri_em']]
        #Limbe_emerge.Mstruc_em = y[self.initial_conditions_mapping[Limbe_emerge]['Mstruc_em']]
        Limbe_emerge.Mstruc_em = Limbe_emerge.Mstruc_em_df(t)
        Limbe_emerge.Naa_em = y[self.initial_conditions_mapping[Limbe_emerge]['Naa_em']]
        Limbe_emerge.Proteines_em = y[self.initial_conditions_mapping[Limbe_emerge]['Proteines_em']]

        ## Phloeme
        phloeme.Csuc_phlo = y[self.initial_conditions_mapping[phloeme]['Csuc_phlo']]
        phloeme.Naa_phlo = y[self.initial_conditions_mapping[phloeme]['Naa_phlo']]
        ## Zone cachee
        Zone_cachee.C_respi_croi = y[self.initial_conditions_mapping[Zone_cachee]['C_respi_croi']]
        Zone_cachee.Camidon_croi = y[self.initial_conditions_mapping[Zone_cachee]['Camidon_croi']]
        Zone_cachee.Csuc_pool_croi = y[self.initial_conditions_mapping[Zone_cachee]['Csuc_pool_croi']]
        Zone_cachee.Ctri_croi = y[self.initial_conditions_mapping[Zone_cachee]['Ctri_croi']]
        Zone_cachee.Fruc_pool_croi = y[self.initial_conditions_mapping[Zone_cachee]['Fruc_pool_croi']]
        Zone_cachee.Fruc_pool_E = y[self.initial_conditions_mapping[Zone_cachee]['Fruc_pool_E']]
        #Zone_cachee.Mstruc_croi = y[self.initial_conditions_mapping[Zone_cachee]['Mstruc_croi']]
        Zone_cachee.Mstruc_croi = Zone_cachee.Mstruc_croi_df(t)
        Zone_cachee.Mstruc_E = y[self.initial_conditions_mapping[Zone_cachee]['Mstruc_E']]
        Zone_cachee.Naa_pool_croi = y[self.initial_conditions_mapping[Zone_cachee]['Naa_pool_croi']]
        Zone_cachee.Prot_pool_croi = y[self.initial_conditions_mapping[Zone_cachee]['Prot_pool_croi']]

        # Fruc_E
        Fruc_E = Zone_cachee.Fruc_E(Zone_cachee.Fruc_pool_E, Zone_cachee.Mstruc_E)
        # Cpool_croi & Npool_croi
        Cpool_croi = Zone_cachee.Cpool_croi(Zone_cachee.Csuc_pool_croi, Zone_cachee.Mstruc_croi)
        Npool_croi = Zone_cachee.Npool_croi(Zone_cachee.Naa_pool_croi, Zone_cachee.Mstruc_croi)

        # croissance
        xE = croissance.xE()
        xb = croissance.xb()
        x_em = croissance.x_em(xb)
        xmax = croissance.xmax(xb)
        xend = croissance.xend(xb)
        x = croissance.x

        ## Times
        t_em = croissance.t_em(parameters.Main.tE)
        tb = croissance.tb(parameters.Main.tE)
        tend = croissance.tend(t_em)
        tmax = croissance.tmax(t_em)

        ## Ltot_max
        Ltot_max = croissance.calculate_Ltot_max(xE, x, croissance.Le, parameters.croissance.Y0)

        ## ce sera la longueur maximale du limbe, et aussi la longueur de la partie émergée
        ligulation = croissance.calculate_ligulation()
        Lem_bis = croissance.calculate_Lem_bis(croissance.Le, Ltot_max, xend, x_em, xmax, xE)
        L_photoS_max = croissance.calculate_L_photoS_max(xE, x, ligulation, Ltot_max, Lem_bis)

        ## Position de Wmax le long du limbe
        Swmax = croissance.Swmax(parameters.Main.N_feuille)
        Lwmax = croissance.calculate_Lwmax(xE, x, Swmax, L_photoS_max)

        ## Maximal leaf width
        Wmax = croissance.Wmax(xE, x, L_photoS_max, parameters.croissance.EC_wmax, parameters.Main.Ksslw, Fruc_E)

        ## Width sheath
        Wlig = croissance.Wlig(xE, x, parameters.croissance.Klig, Wmax)

        ## Width
        W = croissance.W(croissance.L, Lwmax, xE, x, Wmax, parameters.croissance.c1, L_photoS_max, Wlig, parameters.croissance.c2, Ltot_max)

        ## deltaL
        deltaL = croissance.calculate_deltaL(Cpool_croi, xE, x, xend, Ltot_max, croissance.Le, xmax, xb, croissance.L, parameters.Main.Kc, Npool_croi, parameters.Main.Kn, parameters.Main.RERmax, delta_t)

        ## delta width
        deltaW = croissance.calculate_deltaW(croissance.L, Lwmax, Wmax, parameters.croissance.c1, deltaL, L_photoS_max, parameters.croissance.c2, Wlig)

        ## delta surface
        deltaS = croissance.calculate_deltaS(xE, x, xend, deltaL, W, deltaW)

        ## RER
        RER = croissance.calculate_RER(deltaL, croissance.L)

        ## Masse surfacique
        SSLW = croissance.SSLW(xE, x, parameters.Main.min_SSLW, parameters.Main.max_SSLW, Fruc_E, parameters.Main.Ksslw)

        ## DeltaMstruc
        DeltaMstruc = croissance.calculate_DeltaMstruc(x, xE, parameters.Main.alpha_mass_growth, parameters.Main.beta_mass_growth, croissance.L, deltaL, parameters.croissance.ratio_Mstruc, parameters.Main.correction_m, SSLW, deltaS)

        ## RER_m
        RER_m = croissance.calculate_RER_m(DeltaMstruc, croissance.Mstruc_tot)

        ## delta croissance.x (temps C)
        delta_x = croissance.calculate_delta_x(t, Zone_cachee.Csuc_pool_croi, parameters.Main.tE)

        ## longueur de la partie emergee de la feuille
        L_photoS = croissance.calculate_L_photoS(x, x_em, croissance.L, Ltot_max, ligulation, croissance.Lem)

        ## fonction de largeur de la partie de la feuille émergée en fonction de la longueur de la partie émergée
        W_photoS = croissance.W_photoS(L_photoS, Lwmax, Wmax, parameters.croissance.c1, L_photoS_max, Wlig, parameters.croissance.c2)

        ## accroissement de la largeur de la partie émergée de la feuille
        deltaW_photoS  = croissance.calculate_deltaW_photoS(x_em, x, L_photoS, Lwmax, Wmax, parameters.croissance.c1, deltaL, L_photoS_max, parameters.croissance.c2, Wlig)

        ## deltaS_photoS
        deltaS_photoS = croissance.calculate_deltaS_photoS(x_em, x, croissance.L, ligulation, Ltot_max, deltaL, W_photoS, deltaW_photoS)

        ## Longueur de la zone cachée
        L_zone_cachee = croissance.calculate_L_zone_cachee(croissance.L, L_photoS)

        ## pour contrôler que la longueur du limbe + la longueur de la gaine font bien la longueur totale de la feuille
        Lbis = croissance.calculate_Lbis(L_zone_cachee, L_photoS)

        ## Longueur de la partie émergée de la gaine
        L_gaine_em = croissance.calculate_L_gaine_em(x_em, x, L_zone_cachee, croissance.Lem)

        ## Surface de la partie émergée de la gaine.
        S_gaine_em = croissance.calculate_S_gaine_em(x_em, x, L_gaine_em, Wlig)

        # Zone_cachee
        ## Unloading of sucrose from phloem
        Cphlo = phloeme.Cphlo(phloeme.Csuc_phlo, parameters.Main.Mass_plant)
        Zone_cachee.Unloaded_Csuc_phlo = phloeme.Unloading_Csuc_phlo(x, Cphlo, Cpool_croi, parameters.Main.q, Zone_cachee.Mstruc_croi, parameters.Main.conductance, delta_t)

        ## Unloading of AA from phloem
        Nphlo = phloeme.Nphlo(phloeme.Naa_phlo, parameters.Main.Mass_plant)
        Zone_cachee.Unloaded_Naa_phlo = phloeme.Unloading_Naa_phlo(Nphlo, Npool_croi, parameters.Main.q, Zone_cachee.Mstruc_croi, parameters.Main.conductance, delta_t)

        ## Fructan synthesis
        Regul_Sfructanes_gaine = Zone_cachee.Regul_Sfructanes_gaine(parameters.Main.Vmax_Regul_Sfructans, parameters.Main.K_Regul_Sfructans, parameters.Main.n_Regul_Sfructans, Zone_cachee.Unloaded_Csuc_phlo)
        S_Fruc_pool_croi = Zone_cachee.S_Fruc_pool_croi(Cpool_croi, parameters.Main.n_Sfructans, parameters.Main.Vmax_Sfructans, parameters.Main.K_Sfructans, delta_t, Regul_Sfructanes_gaine)

        ## Fructan degradation
        Fruc_croi = Zone_cachee.Fruc_croi(Zone_cachee.Fruc_pool_croi, Zone_cachee.Mstruc_croi)
        D_Fruc_pool_croi = Zone_cachee.D_Fruc_pool_croi(Fruc_croi, Zone_cachee.Mstruc_croi, parameters.Main.K_Dfructan, parameters.Main.Vmax_Dfructan, Cpool_croi, delta_t)

        ## Respi (maintenance)
        Respi = Zone_cachee.Respi(x_em, x, ligulation, Ltot_max, croissance.L, Zone_cachee.Multi_dix_gaine, S_gaine_em, delta_t)

        ## Synthesis Mstruct in C
        taux_suc = Main.taux_suc()
        conv_suc_C = parameters.Main.conv_suc_C
        S_Mstruc_croi_Csuc = Zone_cachee.S_Mstruc_croi_Csuc(DeltaMstruc, taux_suc, conv_suc_C)

         ## Photosynthesis sheath
        PhotoS_gaine = Zone_cachee.PhotoS_gaine(Zone_cachee.Multi_dix_gaine, S_gaine_em, delta_t)

        ## Degradation starch
        Amidon_croi = Zone_cachee.Amidon_croi(Zone_cachee.Mstruc_croi, Zone_cachee.Camidon_croi)
        D_storage_croi = Zone_cachee.D_storage_croi(x_em, x, ligulation, Ltot_max, croissance.L, parameters.Main.delta_Dstorage, Amidon_croi, delta_t)

        ## Synthesis starch
        Tri_croi = Zone_cachee.Tri_croi(Zone_cachee.Mstruc_croi, Zone_cachee.Ctri_croi)
        S_storage_croi = Zone_cachee.S_storage_croi(x_em, x, ligulation, Ltot_max, croissance.L, Tri_croi, parameters.Main.Vmax_storage, parameters.Main.K_storage, delta_t)

        ## Synthesis sucrose
        S_sucrose_croi = Zone_cachee.S_sucrose_croi(x_em, x, ligulation, Ltot_max, croissance.L, Tri_croi, parameters.Main.Vmax_sucrose, parameters.Main.K_sucrose, delta_t)

        ## Respiration croissance
        Respiration = Zone_cachee.Respiration(DeltaMstruc, parameters.Main.correction_m, parameters.Main.taux_respi, conv_suc_C)

        ## Synthesis proteins
        S_Prot_pool_croi = Zone_cachee.S_Prot_pool_croi(parameters.Main.Vmax_Sprotein, Npool_croi, parameters.Main.K_Sprotein, delta_t)

        ## Degradation proteins
        Prot_croi = Zone_cachee.Prot_croi(Zone_cachee.Prot_pool_croi, Zone_cachee.Mstruc_croi)
        D_Prot_pool_croi = Zone_cachee.D_Prot_pool_croi(parameters.Main.delta_Dproteins, Prot_croi, delta_t)

        ## Synthesis Mstruct in N
        taux_aa = Main.taux_aa()
        conv_aa_N = parameters.Main.conv_aa_N
        S_Mstruc_croi_Naa = Zone_cachee.S_Mstruc_croi_Naa(DeltaMstruc, taux_aa, conv_aa_N)

        if Limbe_emerge.has_emerged:
            ## Export of Mstruct towards emerged lamina
            Export_Mstruc_croi = Zone_cachee.Export_Mstruc_croi(croissance.L, ligulation, Ltot_max, parameters.Main.correction_m, deltaS_photoS, SSLW)

            ## Export of AA towards emerged lamina
            Export_Naa_pool_croi = Zone_cachee.Export_Naa_pool_croi(Export_Mstruc_croi, Npool_croi)

            ## Export of sucrose towards emerged lamina
            Export_Csuc_pool_croi = Zone_cachee.Export_Csuc_pool_croi(Cpool_croi, Export_Mstruc_croi)

            # limbe emerge
            ## Export of sucrose towards zone cachee
            C_em = Limbe_emerge.C_em(Limbe_emerge.Mstruc_em, Limbe_emerge.Csuc_em)
            Export_Csuc_em = Limbe_emerge.Export_Csuc_em(x, xend, C_em, Cpool_croi, parameters.Main.q, Limbe_emerge.Mstruc_em, parameters.Main.conductance, delta_t)

            ## Export of AA towards zone cachee
            N_em = Limbe_emerge.N_em(Limbe_emerge.Mstruc_em, Limbe_emerge.Naa_em)
            Export_Naa_em = Limbe_emerge.Export_Naa_em(Limbe_emerge.Mstruc_em, x, xend, N_em, Npool_croi, parameters.Main.q, Zone_cachee.Mstruc_croi, parameters.Main.conductance, delta_t)

            ## Loading sucrose in phloem
            Limbe_emerge.Loaded_Csuc_em = Limbe_emerge.Loading_Csuc_em(xend, x, C_em, Cphlo, parameters.Main.q, Limbe_emerge.Mstruc_em, parameters.Main.conductance, delta_t)

            ## Loading AA in phloem
            Limbe_emerge.Loaded_Naa_em = Limbe_emerge.Loading_Naa_em(xend, x, N_em, Nphlo, parameters.Main.q, Zone_cachee.Mstruc_croi, parameters.Main.conductance, delta_t)

            ## Synthesis proteins
            S_prot_em = Limbe_emerge.S_prot_em(Limbe_emerge.Mstruc_em, parameters.Main.Vmax_Sprotein, N_em, parameters.Main.K_Sprotein, delta_t)

            ## Degradation proteins
            Prot_em = Limbe_emerge.Prot_em(Limbe_emerge.Mstruc_em, Limbe_emerge.Proteines_em)
            D_prot_em = Limbe_emerge.D_prot_em(Limbe_emerge.Mstruc_em, parameters.Main.delta_Dproteins, Prot_em, delta_t)

            ## Synthesis AA
            Tri = Limbe_emerge.Tri(Limbe_emerge.Mstruc_em, Limbe_emerge.Ctri_em)
            S_amino_acids = Limbe_emerge.S_amino_acids(parameters.Main.Vmax_aa, Tri, parameters.Main.K_Saa_trioses, delta_t)
            S_amino_acidsC_limbe = Limbe_emerge.S_amino_acidsC_limbe(S_amino_acids, parameters.Main.conv_aa_C, parameters.Main.conv_aa_N)

            ## Photosynthesis emerged lamina
            PhotoS = Limbe_emerge.PhotoS(x_em, x, Limbe_emerge.Photosynthesis_limbe, croissance.S_photoS, delta_t)

            ## Respiration
            Respi = Limbe_emerge.Respi(x_em, x, Limbe_emerge.Multi_dix_limbe, croissance.S_photoS, delta_t)

            ## Synthesis sucrose
            S_sucrose = Limbe_emerge.S_sucrose(x_em, x, Tri, parameters.Main.Vmax_sucrose, parameters.Main.K_sucrose, delta_t)

            ## Synthesis starch
            S_storage = Limbe_emerge.S_storage(x_em, x, Tri, parameters.Main.Vmax_storage, parameters.Main.K_storage, delta_t)

            ## Degradation starch
            Amidon = Limbe_emerge.Amidon(Limbe_emerge.Mstruc_em, Limbe_emerge.Camidon_em)
            D_storage = Limbe_emerge.D_storage(x_em, x, parameters.Main.delta_Dstorage, Amidon, delta_t)

            ## Synthesis fructans
            Regul_Sfructanes_limbe = Limbe_emerge.Regul_Sfructanes_limbe(parameters.Main.Vmax_Regul_Sfructans, parameters.Main.K_Regul_Sfructans, parameters.Main.n_Regul_Sfructans, Limbe_emerge.Loaded_Csuc_em, Export_Csuc_em)
            S_fruc_em = Limbe_emerge.S_fruc_em(x_em, x, C_em, parameters.Main.n_Sfructans, parameters.Main.Vmax_Sfructans, parameters.Main.K_Sfructans, delta_t, Regul_Sfructanes_limbe)

            ## Degradation fructans
            fruc_em = Limbe_emerge.Fruc_em(Limbe_emerge.Mstruc_em, Limbe_emerge.Cfruc_em)
            D_fruc_em = Limbe_emerge.D_fruc_em(x_em, x, fruc_em, parameters.Main.K_Dfructan, parameters.Main.Vmax_Dfructan, C_em, delta_t)

        else:
            # Zone_cachee
            ## Export of Mstruct towards emerged lamina
            Export_Mstruc_croi = 0

            ## Export of AA towards emerged lamina
            Export_Naa_pool_croi = 0

            ## Export of sucrose towards emerged lamina
            Export_Csuc_pool_croi = 0

            # limbe emerge
            ## Export of sucrose towards zone cachee
            Export_Csuc_em = 0

            ## Export of AA towards zone cachee
            Export_Naa_em = 0

            ## Loading sucrose in phloem
            Limbe_emerge.Loaded_Csuc_em = 0

            ## Loading AA in phloem
            Limbe_emerge.Loaded_Naa_em = 0


        # compute the derivative
        ## croissance
        y_derivatives[self.initial_conditions_mapping[croissance]['L']] = croissance.calculate_L_derivative(deltaL)
        y_derivatives[self.initial_conditions_mapping[croissance]['Le']] = croissance.calculate_Le_derivative(x, xE, deltaL)
        y_derivatives[self.initial_conditions_mapping[croissance]['Lem']] = croissance.calculate_Lem_derivative(x, x_em, deltaL)
        y_derivatives[self.initial_conditions_mapping[croissance]['M_em']] = croissance.calculate_M_em_derivative(x, x_em, DeltaMstruc)
        y_derivatives[self.initial_conditions_mapping[croissance]['Mstruc_tot']] = croissance.calculate_Mstruc_tot_derivative(DeltaMstruc)
        y_derivatives[self.initial_conditions_mapping[croissance]['S_photoS']] = croissance.calculate_S_photoS_derivative(deltaS_photoS)
        y_derivatives[self.initial_conditions_mapping[croissance]['S_tot']] = croissance.calculate_S_tot_derivative(deltaS)
        y_derivatives[self.initial_conditions_mapping[croissance]['W_photoS_bis']] = croissance.calculate_W_photoS_bis_derivative(deltaW_photoS)
        y_derivatives[self.initial_conditions_mapping[croissance]['Wbis']] = croissance.calculate_Wbis_derivative(deltaW)
        y_derivatives[self.initial_conditions_mapping[croissance]['x']] = croissance.calculate_x_derivative(delta_x)

        ## zone cachee
        y_derivatives[self.initial_conditions_mapping[Zone_cachee]['Ctri_croi']] = Zone_cachee.calculate_Ctri_croi_derivative(x_em, x, ligulation, Ltot_max, croissance.L, PhotoS_gaine, S_sucrose_croi, S_storage_croi, Zone_cachee.Mstruc_croi)
        y_derivatives[self.initial_conditions_mapping[Zone_cachee]['Camidon_croi']] = Zone_cachee.calculate_Camidon_croi_derivative(x_em, x, ligulation, Ltot_max, croissance.L, S_storage_croi, D_storage_croi, Zone_cachee.Mstruc_croi)
        y_derivatives[self.initial_conditions_mapping[Zone_cachee]['Fruc_pool_croi']] = Zone_cachee.calculate_Fruc_pool_croi_derivative(S_Fruc_pool_croi, Zone_cachee.Mstruc_croi, D_Fruc_pool_croi)
        y_derivatives[self.initial_conditions_mapping[Zone_cachee]['Fruc_pool_E']] = Zone_cachee.calculate_Fruc_pool_E_derivative(x, xE, S_Fruc_pool_croi, Zone_cachee.Mstruc_croi, D_Fruc_pool_croi)
        y_derivatives[self.initial_conditions_mapping[Zone_cachee]['Mstruc_E']] = Zone_cachee.calculate_Mstruc_E_derivative(x, xE, S_Mstruc_croi_Naa, parameters.Main.MN, S_Mstruc_croi_Csuc, parameters.Main.MC, parameters.Main.conv_aa_C, parameters.Main.conv_aa_N)
        y_derivatives[self.initial_conditions_mapping[Zone_cachee]['Csuc_pool_croi']] = Zone_cachee.calculate_Csuc_pool_croi_derivative(Zone_cachee.Unloaded_Csuc_phlo, S_Mstruc_croi_Csuc, Respiration, S_Fruc_pool_croi, Zone_cachee.Mstruc_croi, D_Fruc_pool_croi, Export_Csuc_pool_croi, Export_Csuc_em, S_sucrose_croi, D_storage_croi, Respi)
        y_derivatives[self.initial_conditions_mapping[Zone_cachee]['C_respi_croi']] = Zone_cachee.calculate_C_respi_croi_derivative(Respiration, parameters.Main.Msuc)
        #y_derivatives[self.initial_conditions_mapping[Zone_cachee]['Mstruc_croi']] = Zone_cachee.calculate_Mstruc_croi_derivative(x_em, x, S_Mstruc_croi_Naa, parameters.Main.MN, S_Mstruc_croi_Csuc, parameters.Main.MC, parameters.Main.conv_aa_C, parameters.Main.conv_aa_N, Export_Mstruc_croi)
        y_derivatives[self.initial_conditions_mapping[Zone_cachee]['Naa_pool_croi']] = Zone_cachee.calculate_Naa_pool_croi_derivative(Zone_cachee.Unloaded_Naa_phlo, S_Mstruc_croi_Naa, S_Prot_pool_croi, Zone_cachee.Mstruc_croi, D_Prot_pool_croi, Export_Naa_pool_croi, Export_Naa_em)
        y_derivatives[self.initial_conditions_mapping[Zone_cachee]['Prot_pool_croi']] = Zone_cachee.calculate_Prot_pool_croi_derivative(S_Prot_pool_croi, Zone_cachee.Mstruc_croi, D_Prot_pool_croi)

        ## Limbe emerge
        if Limbe_emerge.has_emerged:
            y_derivatives[self.initial_conditions_mapping[Limbe_emerge]['Ctri_em']] = Limbe_emerge.calculate_Ctri_em_derivative(x_em, x, PhotoS, S_sucrose, S_storage, S_amino_acidsC_limbe, Limbe_emerge.Mstruc_em)
            y_derivatives[self.initial_conditions_mapping[Limbe_emerge]['Camidon_em']] = Limbe_emerge.calculate_Camidon_em_derivative(Limbe_emerge.Mstruc_em, S_storage, D_storage)
            y_derivatives[self.initial_conditions_mapping[Limbe_emerge]['Cfruc_em']] = Limbe_emerge.calculate_Cfruc_em_derivative(S_fruc_em, D_fruc_em, Limbe_emerge.Mstruc_em)
            y_derivatives[self.initial_conditions_mapping[Limbe_emerge]['Csuc_em']] = Limbe_emerge.calculate_Csuc_em_derivative(x_em, x, Export_Csuc_pool_croi, Export_Csuc_em, S_sucrose, D_storage, Limbe_emerge.Mstruc_em, Limbe_emerge.Loaded_Csuc_em, Respi, D_fruc_em, S_fruc_em)
            #y_derivatives[self.initial_conditions_mapping[Limbe_emerge]['Mstruc_em']] = Limbe_emerge.calculate_Mstruc_em_derivative(Export_Mstruc_croi)
            y_derivatives[self.initial_conditions_mapping[Limbe_emerge]['Naa_em']] = Limbe_emerge.calculate_Naa_em_derivative(Limbe_emerge.Mstruc_em, Export_Naa_pool_croi, Export_Naa_em, Limbe_emerge.Loaded_Naa_em, S_prot_em, D_prot_em, S_amino_acids)
            y_derivatives[self.initial_conditions_mapping[Limbe_emerge]['Proteines_em']] = Limbe_emerge.calculate_Proteines_em_derivative(S_prot_em, D_prot_em, Limbe_emerge.Mstruc_em)
        else:
            y_derivatives[self.initial_conditions_mapping[Limbe_emerge]['Ctri_em']] = 0
            y_derivatives[self.initial_conditions_mapping[Limbe_emerge]['Camidon_em']] = 0
            y_derivatives[self.initial_conditions_mapping[Limbe_emerge]['Cfruc_em']] = 0
            y_derivatives[self.initial_conditions_mapping[Limbe_emerge]['Csuc_em']] = 0
            #y_derivatives[self.initial_conditions_mapping[Limbe_emerge]['Mstruc_em']] = Limbe_emerge.calculate_Mstruc_em_derivative(Export_Mstruc_croi)
            y_derivatives[self.initial_conditions_mapping[Limbe_emerge]['Naa_em']] = 0
            y_derivatives[self.initial_conditions_mapping[Limbe_emerge]['Proteines_em']] = 0

        ## phloeme (FLUX AVEC ORGANES IGNORES POUR PHLOEME)
        y_derivatives[self.initial_conditions_mapping[phloeme]['Csuc_phlo']] = phloeme.calculate_Csuc_phlo_derivative(x_em, x, t)
        y_derivatives[self.initial_conditions_mapping[phloeme]['Naa_phlo']] = phloeme.calculate_Naa_phlo_derivative(t_em, t)

        derivatives_logger = logging.getLogger('growthwheat.derivatives')
        if derivatives_logger.isEnabledFor(logging.DEBUG):
            self._log_compartments(t, y_derivatives, Simulation.LOGGERS_NAMES['derivatives'])


        return y_derivatives

    def _update_population(self, compartments_values):
        """Update the state of :attr:`population` from the values in `compartments_values`.
        """
        logger = logging.getLogger(__name__)
        logger.debug('Updating the state of the population...')
        for model_object, compartments in self.initial_conditions_mapping.iteritems():
            for compartment_name, compartment_index in compartments.iteritems():
                setattr(model_object, compartment_name, compartments_values[compartment_index])
        logger.debug('Updating the state of the population DONE')

    def postprocessings(self):
        """
        Compute:

            * intermediate variables (see :attr:`Simulation:PLANTS_INTERMEDIATE_VARIABLES`, :attr:`Simulation:AXES_INTERMEDIATE_VARIABLES`, :attr:`Simulation:PHYTOMERS_INTERMEDIATE_VARIABLES`, :attr:`Simulation:ORGANS_INTERMEDIATE_VARIABLES` and :attr:`Simulation:ELEMENTS_INTERMEDIATE_VARIABLES`),
            * fluxes (see :attr:`Simulation:PLANTS_FLUXES`, :attr:`Simulation:AXES_FLUXES`, :attr:`Simulation:PHYTOMERS_FLUXES`, :attr:`Simulation:ORGANS_FLUXES` and :attr:`Simulation:ELEMENTS_FLUXES`),
            * and integrative variables (see :attr:`Simulation:PLANTS_INTEGRATIVE_VARIABLES`, :attr:`Simulation:AXES_INTEGRATIVE_VARIABLES`, :attr:`Simulation:PHYTOMERS_INTEGRATIVE_VARIABLES`, :attr:`Simulation:ORGANS_INTEGRATIVE_VARIABLES` and :attr:`Simulation:ELEMENTS_INTEGRATIVE_VARIABLES`),

        from :attr:`_solver_output` and format them to :class:`dataframes <pandas.DataFrame>`.

        :Returns:
            :class:`dataframes <pandas.DataFrame>` of post-processing outputs at each scale:

                * plant (see :attr:`Simulation:PLANTS_ALL_VARIABLES`)
                * axis (see :attr:`Simulation:AXES_ALL_VARIABLES`)
                * phytomer (see :attr:`Simulation:PHYTOMERS_ALL_VARIABLES`)
                * organ (see :attr:`Simulation:ORGANS_ALL_VARIABLES`)
                * and element (see :attr:`Simulation:ELEMENTS_ALL_VARIABLES`)

        :Returns Type:
            :class:`tuple` of :class:`pandas.DataFrame`

        """
        logger = logging.getLogger(__name__)
        logger.debug('Formatting of outputs...')

        solver_output_transposed = self._solver_output.T

        croissance_df = pd.DataFrame(columns=Simulation.CROISSANCE)
        Zone_cachee_df = pd.DataFrame(columns=Simulation.ZONE_CACHEE)
        Limbe_emerge_df = pd.DataFrame(columns=Simulation.LIMBE_EMERGE)
        phloeme_df = pd.DataFrame(columns=Simulation.PHLOEME)

        for organ in self.population.organs:
            if issubclass(organ.__class__, model.Limbe_emerge):
                Limbe_emerge = organ
            elif issubclass(organ.__class__, model.phloeme):
                phloeme = organ
            elif issubclass(organ.__class__, model.Zone_cachee):
                Zone_cachee = organ
            elif issubclass(organ.__class__, model.croissance):
                croissance = organ

        croissance_df['t'] = self._time_grid
        Zone_cachee_df['t'] = self._time_grid
        Limbe_emerge_df['t'] = self._time_grid
        phloeme_df['t'] = self._time_grid


        # format phloem outputs
        phloeme_df['Csuc_phlo'] = solver_output_transposed[self.initial_conditions_mapping[phloeme]['Csuc_phlo']]
        phloeme_df['Naa_phlo'] = solver_output_transposed[self.initial_conditions_mapping[phloeme]['Naa_phlo']]

        # format croissance output
        croissance_df['L'] = solver_output_transposed[self.initial_conditions_mapping[croissance]['L']]
        croissance_df['Le'] = solver_output_transposed[self.initial_conditions_mapping[croissance]['Le']]
        croissance_df['Lem'] = solver_output_transposed[self.initial_conditions_mapping[croissance]['Lem']]
        croissance_df['M_em'] = solver_output_transposed[self.initial_conditions_mapping[croissance]['M_em']]
        croissance_df['Mstruc_tot'] = solver_output_transposed[self.initial_conditions_mapping[croissance]['Mstruc_tot']]
        croissance_df['S_photoS'] = solver_output_transposed[self.initial_conditions_mapping[croissance]['S_photoS']]
        croissance_df['S_tot'] = solver_output_transposed[self.initial_conditions_mapping[croissance]['S_tot']]
        croissance_df['W_photoS_bis'] = solver_output_transposed[self.initial_conditions_mapping[croissance]['W_photoS_bis']]
        croissance_df['Wbis'] = solver_output_transposed[self.initial_conditions_mapping[croissance]['Wbis']]
        croissance_df['x'] = solver_output_transposed[self.initial_conditions_mapping[croissance]['x']]

        # format Zone cachee output
        Zone_cachee_df['Ctri_croi'] = solver_output_transposed[self.initial_conditions_mapping[Zone_cachee]['Ctri_croi']]
        Zone_cachee_df['Camidon_croi'] = solver_output_transposed[self.initial_conditions_mapping[Zone_cachee]['Camidon_croi']]
        Zone_cachee_df['Fruc_pool_croi'] = solver_output_transposed[self.initial_conditions_mapping[Zone_cachee]['Fruc_pool_croi']]
        Zone_cachee_df['Fruc_pool_E'] = solver_output_transposed[self.initial_conditions_mapping[Zone_cachee]['Fruc_pool_E']]
        Zone_cachee_df['Mstruc_E'] = solver_output_transposed[self.initial_conditions_mapping[Zone_cachee]['Mstruc_E']]
        Zone_cachee_df['Csuc_pool_croi'] = solver_output_transposed[self.initial_conditions_mapping[Zone_cachee]['Csuc_pool_croi']]
        Zone_cachee_df['C_respi_croi'] = solver_output_transposed[self.initial_conditions_mapping[Zone_cachee]['C_respi_croi']]
        #Zone_cachee_df['Mstruc_croi'] = solver_output_transposed[self.initial_conditions_mapping[Zone_cachee]['Mstruc_croi']]
        Zone_cachee_df['Naa_pool_croi'] = solver_output_transposed[self.initial_conditions_mapping[Zone_cachee]['Naa_pool_croi']]
        Zone_cachee_df['Prot_pool_croi'] = solver_output_transposed[self.initial_conditions_mapping[Zone_cachee]['Prot_pool_croi']]

        Cpool_croi = map(Zone_cachee.Cpool_croi, Zone_cachee_df['Csuc_pool_croi'], Zone_cachee_df['Mstruc_croi'])
        Zone_cachee_df['D_Fruc_pool_croi'] = map(Zone_cachee.D_Fruc_pool_croi, Zone_cachee_df['Fruc_pool_croi'], Zone_cachee_df['Mstruc_croi'], [parameters.Main.K_Dfructan]* len(self._time_grid), [parameters.Main.Vmax_Dfructan]* len(self._time_grid), Cpool_croi, [3600]* len(self._time_grid))

        # format Limbe emerge output
        Limbe_emerge_df['Ctri_em'] = solver_output_transposed[self.initial_conditions_mapping[Limbe_emerge]['Ctri_em']]
        Limbe_emerge_df['Camidon_em'] = solver_output_transposed[self.initial_conditions_mapping[Limbe_emerge]['Camidon_em']]
        Limbe_emerge_df['Cfruc_em'] = solver_output_transposed[self.initial_conditions_mapping[Limbe_emerge]['Cfruc_em']]
        Limbe_emerge_df['Csuc_em'] = solver_output_transposed[self.initial_conditions_mapping[Limbe_emerge]['Csuc_em']]
        #Limbe_emerge_df['Mstruc_em'] = solver_output_transposed[self.initial_conditions_mapping[Limbe_emerge]['Mstruc_em']]
        Limbe_emerge_df['Naa_em'] = solver_output_transposed[self.initial_conditions_mapping[Limbe_emerge]['Naa_em']]
        Limbe_emerge_df['Proteines_em'] = solver_output_transposed[self.initial_conditions_mapping[Limbe_emerge]['Proteines_em']]

        # set the order of the columns
        croissance_df = croissance_df.reindex_axis(Simulation.CROISSANCE, axis=1, copy=False)
        Zone_cachee_df = Zone_cachee_df.reindex_axis(Simulation.ZONE_CACHEE, axis=1, copy=False)
        Limbe_emerge_df = Limbe_emerge_df.reindex_axis(Simulation.LIMBE_EMERGE, axis=1, copy=False)
        phloeme_df = phloeme_df.reindex_axis(Simulation.PHLOEME, axis=1, copy=False)

        # sort the rows by the columns
        croissance_df.sort_index(by=Simulation.CROISSANCE, inplace=True)
        Zone_cachee_df.sort_index(by=Simulation.ZONE_CACHEE, inplace=True)
        Limbe_emerge_df.sort_index(by=Simulation.LIMBE_EMERGE, inplace=True)
        phloeme_df.sort_index(by=Simulation.PHLOEME, inplace=True)

        # infer the right types of the columns
        croissance_df = croissance_df.convert_objects(copy=False)
        Zone_cachee_df = Zone_cachee_df.convert_objects(copy=False)
        Limbe_emerge_df = Limbe_emerge_df.convert_objects(copy=False)
        phloeme_df = phloeme_df.convert_objects(copy=False)

        logger.debug('Formatting of outputs DONE')

        return croissance_df, Zone_cachee_df, Limbe_emerge_df, phloeme_df

