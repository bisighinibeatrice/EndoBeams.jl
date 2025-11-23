function run_simulation!(conf::SimulationConfiguration, params::SimulationParams, inter::Union{Nothing, Interaction}=nothing)

    # Unpack relevant parameters for simulation control
    @unpack stop_on_energy_threshold, verbose, record_timings, output_dir, min_timestep, initial_timestep, 
    simulation_end_time, output_timestep, accelerate_after_success_it, max_timestep, max_iterations, energy_threshold, 
    tcompt_max, stop_long_simulation = params

    # If timing information is required, initialize timing tools
    if record_timings
        TimerOutputs.enable_debug_timings(@__MODULE__)
        global timer = TimerOutputs.TimerOutput()  # Create a new TimerOutput instance        
        reset_timer!()
    end 

    # Setup variables for the solver and for the visualization 
    state, vtkdata = setup_state_simulation(conf, inter, output_dir)

    # Setup the linear solver 
    solver = init(:MKL)     

    # Initialize simulation time variables
    tⁿ = 0.0                     # Current simulation time
    write_t = 0.0               # Time of the last output write
    start = time()              # Start time for computational timing
    write_counter = 0         # Counter for output writes
    Δt = initial_timestep                 # Initial time step size
    successive_success_it = 0   # Consecutive successful iterations counter
    step = 1                    # Simulation step counter
    t⁺_old = 0 # parameters to change SDF

    # Initialize visualization variables and save initial configuration
    write_VTK(write_counter, step, tⁿ, conf, vtkdata)
    write_counter += 1     

    # Initialize simulation state variables
    initialize_state_simulation!(conf, state, params)

    # Temporal loop: run simulation until the end time is reached
    while simulation_end_time - tⁿ > 1.0e-10

        tⁿ⁺¹ = tⁿ + Δt    # Define current and next time points
        tcomp = time() - start  # Compute elapsed computation time

        # Display step information if verbose mode is enabled
        if verbose
            println("------------------------------------------------")
            @printf "%8s\t%8s\t%4s\t%4s\n" "tⁿ" "Δt" "step" "tcomp"
            @printf "%1.2e\t%1.2e\t%4d\t%1.2f\n" tⁿ⁺¹ Δt step tcomp
            println("------------------------------------------------")
        end 
        
        # Solve the system for the current time step
        n_it = solve_step_dynamics!(conf, state, params, inter, tⁿ⁺¹, Δt, solver)


        t⁺ = convert(Int, round(tⁿ⁺¹, RoundDown))
        if !isnothing(inter) && isa(inter.master, DiscreteSignedDistanceField)
            if inter.master.flag_load_from_file_iterative && t⁺ != t⁺_old && t⁺ <= 500
                println("Updating SDF $t⁺")
                surface_master = DiscreteSignedDistanceField("stent_3PB/inputSDF/sdf_$t⁺.vtk", false, true)
                                println("Updated SDF")
                inter = RigidInteraction(surface_master, inter.slave, inter.properties)
            end 
        end 

        # Check convergence status
        if n_it > max_iterations
            # If the solver fails to converge within the maximum allowed iterations
            verbose && printstyled("Time step did not converge...\n"; color = :yellow)
            Δt = Δt / 2  # Halve the time step size
            verbose && printstyled((@sprintf "Δt decreased: Δt = %1.2e\n" Δt); color = :yellow)

            # If the time step is too small, abort the simulation
            if Δt < min_timestep 
                printstyled((@sprintf "Δt < %1.2e: aborting simulation...\n" min_timestep); color = :red)
                error("Simulation aborted.")
                break
            end

            # Revert to the previous state and reattempt convergence
            update_not_converged!(conf, state)
            successive_success_it = 0
            continue
        else
            # Solver converged
            verbose && printstyled("Time step converged!\n"; color = :green) 
            successive_success_it += 1
        end

        # Increase the time step if multiple consecutive successful steps occur
        if (successive_success_it ≥ accelerate_after_success_it) && n_it <= 8
            Δtold = Δt
            Δt = min(Δt * 1.5, max_timestep)  # Gradually increase the time step size
            verbose && Δt > Δtold && printstyled((@sprintf "Δt increased: Δt = %1.2e\n" Δt); color = :green)
            successive_success_it = 0
        end

        if  tⁿ⁺¹ > write_t + output_timestep  ||  tⁿ⁺¹ ≈ write_t + output_timestep
            write_VTK(write_counter, step, tⁿ, conf, vtkdata)
            write_counter += 1     
            write_t = tⁿ⁺¹
        end  

        # Check computational time to avoid excessively long simulations
        if tcomp > tcompt_max && stop_long_simulation
            printstyled((@sprintf "tcomp > %1.2e: aborting simulation...\n" tcompt_max); color = :red)
            printstyled("Simulation aborted."; color = :red)
            break
        end
        
        # Update the system to the converged state for the next time step
        update_converged!(conf, state)

        # Advance to the next time step
        tⁿ = tⁿ⁺¹
        step += 1
        t⁺_old = t⁺

        # Terminate simulation early if the energy threshold is met
        if stop_on_energy_threshold && energyⁿ⁺¹.kinetic_energy < energy_threshold
            printstyled((@sprintf "Energy threshold reached! Ending simulation at tcomp = %1.2e.\n" tcomp); color = :green) 
            break
        end
        
    end 
    
    # Show sections time
    if record_timings
        print_timer()
    end 
    
end
