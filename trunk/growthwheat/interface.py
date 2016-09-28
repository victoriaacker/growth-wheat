import converter, simulation


def initialize(g, inputs):
    """Update `g` from `inputs` in place.

    :Parameters:
        - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to update.

        - `inputs` (:class:`dict` of :class:`pandas.DataFrame`) - The inputs of the model. Required keys are 'hiddenzone_inputs', 'organ_inputs' and 'root_inputs'

    """
    all_inputs_dict = converter.from_dataframes(**inputs)
    converter.update_MTG(g, inputs=all_inputs_dict)

def run(g, dt, copy=True):
    """Run the model from `g` and `dt`.

    :Parameters:
        - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to run the model on.

        - `dt` (:class:`int`) - the delta t of the simulation (in seconds).

        - `copy` (:class:`bool`) - If `True.`, return a new instance of :class:`g <openalea.mtg.mtg.MTG>` (the default).
          If `False`, update `g` in place.

    :Returns:
        A copy of `g` if `copy` is True (the default), the same instance otherwise.

    :Returns Type:
        :class:`openalea.mtg.mtg.MTG`

    .. seealso:: see :func:`adelgrowthwheat.adelgeom.interface.add_metamer` for the type signature of the function. #TODO: update
    """
    # create a simulation object
    simulation_ = simulation.Simulation(delta_t=dt)
    # convert the MTG to simulation inputs format
    inputs = converter.from_MTG(g)
    # initialize the simulation from the inputs
    simulation_.initialize(inputs)
    # run the simulation
    simulation_.run()
    # update the MTG from the outputs of the simulation, adding new metamer(s) if needed
    converter.update_MTG(g, outputs=simulation_.outputs)
    # Conversion of results to dataframes
    growthwheat_hiddenzones_outputs, growthwheat_organs_outputs, growthwheat_roots_outputs = converter.to_dataframes(simulation_.outputs)
    # return the updated MTG
    return g, growthwheat_hiddenzones_outputs, growthwheat_organs_outputs, growthwheat_roots_outputs