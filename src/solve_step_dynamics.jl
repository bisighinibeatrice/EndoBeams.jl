function solve_step_dynamics!(conf::SimulationConfiguration, state::SimulationState, params::SimulationParams, inter::Union{Nothing, Interaction}, tⁿ⁺¹, Δt, solver) 
    
    # Update external loads at the current time step
    update_loads!(conf, state, tⁿ⁺¹)
    update_boundary_conditions!(conf, tⁿ⁺¹)

    # Predictor phase: Estimate the solution for the current step
    @timeit_debug "Predictor" predictor!(conf, state, params, Δt, solver)

    # Corrector phase: Iteratively refine the solution until convergence
    @timeit_debug "Corrector" k = corrector!(conf, state, params, inter, Δt, solver)

    return k
end
