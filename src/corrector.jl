function corrector!(conf::SimulationConfiguration, state::SimulationState, params::SimulationParams, inter::Union{Nothing, Interaction}, Δt, solver)
    
    # Unpack parameters for the corrector loop
    @unpack tolerance_residual, tolerance_displacement, max_iterations, α, β, γ = params
    
    # Initialize global corrector variables (sets up solver matrices and forces)
    @timeit_debug "Initialize Global Corrector" initialize_global_corrector!(conf, state)
    
    # Initialize iteration variables
    k = 1                   # Current iteration number
    ΔD_norm = Inf           # Displacement increment norm
    res_norm = Inf          # Residual norm
        
    # Main loop for corrector iteration
    while (res_norm > tolerance_residual || ΔD_norm > tolerance_displacement) && k ≤ max_iterations

        # Assemble element contributions to the stiffness matrix and residuals
        @timeit_debug "Assemble" assemble!(conf, state, params)
        if isa(inter, Interaction) 
            @timeit_debug "Assemble Contact" assemble_contact!(conf, state, inter.master, inter.slave, inter.properties, inter, params)
        end  
        if !isnothing(conf.constraints)
            @timeit_debug "Assemble Constraints" assembly_constraints!(conf, state)
        end 

        # Compute tangent matrix and residuals for corrector
        @timeit_debug "Compute Tangent and Residuals" compute_tangent_and_residuals_corrector!(state, Δt, α, β, γ)

        # Apply cylindrical to carthesian coordinate system
        if conf.bcs.use_cylindrical_coords 
            @timeit_debug "Apply cylindrical to carthesian coordinate system" apply_cylindrical_coordinate_system!(conf, state) 
        end
        
        # Apply boundary conditions to restrict degrees of freedom
        @timeit_debug "Apply Boundary Conditions" apply_boundary_conditions!(conf, state)
        
        # Extract free dof tangent matrix and residual (for solving the free system)
        @timeit_debug "Extract Free DOFs" extract_free_dofs!(conf, state)

        # Revert to carthesian to cylindrical coordinate system
        if conf.bcs.use_cylindrical_coords 
            @timeit_debug "Revert to carthesian to cylindrical coordinate system" revert_to_carthesian_coordinate_system!(conf, state) 
        end
        
        # Solve for displacement increments at the free dofs
        @timeit_debug "Solve Free DOFs" solve_free_dofs!(conf, state, solver)
        
        # Update global displacement, velocity, and acceleration variables after solving
        @timeit_debug "Compute Global Corrector" compute_global_corrector!(state, Δt, β, γ)
        
        # Update local variables based on the computed global values
        @timeit_debug "Compute Local Corrector" compute_local_corrector!(conf, state, γ, β, Δt)

        # Compute the norms of the residual and displacement increments
        @timeit_debug "Compute Norms" res_norm, ΔD_norm = compute_norms_corrector(conf, state)
 
        # Print iteration information if verbose flag is set
        if params.verbose
            k == 1 && @printf "%4s\t%8s\t%8s\n" "iter" "‖res‖" "‖ΔD‖"
            @printf "%4d\t%1.2e\t%1.2e\n"  k res_norm ΔD_norm
        end
 
        # Increment the iteration counter
        k += 1
        
    end 

    # Release solver resources after convergence or max iterations
    release!(solver)
    
    return k
end
