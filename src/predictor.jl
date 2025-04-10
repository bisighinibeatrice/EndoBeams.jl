function predictor!(conf::SimulationConfiguration, state::SimulationState, params::SimulationParams, Δt, solver)

    # Unpack parameters
    @unpack tolerance_residual, tolerance_displacement, max_iterations, α, β, γ = params

    # Initialize local predictor (displacements, velocities, and accelerations)
    initialize_global_predictor!(conf, state)

    # Compute tangent matrix and residuals (needed for prediction step)
    compute_tangent_and_residuals_predictor!(state, Δt, α, β, γ)

    # Apply boundary conditions (fixing displacement or forces)
    apply_boundary_conditions!(conf, state)

    # Compute free tangent matrix and residual for the unconstrained degrees of freedom
    extract_free_dofs!(conf, state)

    # Solve the system for free degrees of freedom (computes displacement, velocity, and acceleration)
    solve_free_dofs!(conf, state, solver)

    # Update global displacement, velocity, and acceleration variables (predicted values)
    compute_global_predictor!(state, Δt, β, γ)

    # Update local variables with the global predicted values
    compute_local_predictor!(conf, state, γ, β, Δt)

end

