"Preallocates and initialises all the variables and starts the temporal loop"
function solver!(conf, params)

    nodes = conf.nodes
    
    # -------------------------------------------------------------------------------------------
    # INITIALISATION
    # -------------------------------------------------------------------------------------------
  
    @unpack stop_on_energy_threshold, verbose, record_timings, output_dir, min_Δt, ini_Δt, tᵉⁿᵈ, Δt_plot, accelerate_after_success_it, max_Δt, max_it, energy_threshold = params

    if record_timings
        TimerOutputs.enable_debug_timings(@__MODULE__)
        reset_timer!()
    end 

    @timeit_debug "Solver initialization" begin
        
        # initialization of time loop variables
        t = 0.
        write_t = t
        start = time()
        write_counter = 0
        Δt = ini_Δt
        Δt_plot = Δt_plot
        successive_success_it = 0
        step = 1
        
        # preallocation of the vectors and matrices used in the computaion
        solⁿ, solⁿ⁺¹, nodes_sol, matrices, energy, vtkdata = solver_initialisation(conf, output_dir)

        write_VTK(write_counter, 0, t, conf, energy, vtkdata)
        write_counter += 1
        
        # Linear solver
        
        solver = init(:MKL)
        
    end 
    
    # -------------------------------------------------------------------------------------------
    # TIME LOOP 
    # -------------------------------------------------------------------------------------------

    @timeit_debug "Solver time loop" begin
        
        while t<tᵉⁿᵈ
            
            @timeit_debug "Dynamic step" begin
                
                tⁿ = t
                tⁿ⁺¹ = tⁿ + Δt
                tcomp = time() - start
                
                if verbose 
                    println("------------------------------------------------")
                    @printf "%8s\t%8s\t%4s\t%4s\n" "t" "Δt" "step" "tcomp"
                    @printf "%1.2e\t%1.2e\t%4d\t%1.2f\n" tⁿ⁺¹ Δt step tcomp
                    println("------------------------------------------------")
                end 
                
                # solve system @tⁿ
                n_it = solve_step_dynamics!(conf, tⁿ⁺¹, Δt, solⁿ, solⁿ⁺¹, nodes_sol, matrices, energy, solver, params)

                # if not converged, halve the time step and re-solve the system until it converges
                if n_it > max_it

                    verbose && printstyled("Time step did not converge...\n"; color = :yellow)
                    
                    Δt = Δt/2
                    
                    verbose && printstyled((@sprintf "Δt decreased: Δt=%1.2e\n" Δt); color = :yellow)

                    if Δt<min_Δt
                        printstyled((@sprintf "Δt < %1.2e: aborting simulation...\n" min_Δt); color = :red)
                        break
                    end
                    
                    solⁿ⁺¹ = deepcopy(solⁿ)
                    update_nodes_not_converged!(nodes)
                    successive_success_it = 0

                    continue

                else

                    verbose && printstyled("Time step converged!\n"; color = :green) 
                    successive_success_it += 1

                end
                
                
                if successive_success_it ≥ accelerate_after_success_it
                    
                    Δtold = Δt
                    Δt = min(Δt*1.5, max_Δt)
                    verbose && Δt > Δtold && printstyled((@sprintf "Δt increased: Δt=%1.2e\n" Δt); color = :green)
                    successive_success_it = 0
                    
                end
                
            end
            
            @timeit_debug "Write results to file" begin
                
                # save VTK with frequency 1/Δt_plot: 
                if  tⁿ⁺¹ > write_t + Δt_plot  ||  tⁿ⁺¹ ≈ write_t + Δt_plot
                    write_VTK(write_counter, step, tⁿ⁺¹, conf, energy, vtkdata)
                    write_counter += 1
                    write_t = tⁿ⁺¹
                end
 
                
            end 
            
            @timeit_debug "Update variables" begin
                
                # update variables after the step
                for f in fieldnames(typeof(solⁿ⁺¹))
                    getfield(solⁿ, f) .= getfield(solⁿ⁺¹, f)
                end
                update_nodes_converged!(nodes, solⁿ, matrices)
                t = tⁿ⁺¹
                step += 1
                
            end 

            if stop_on_energy_threshold && energy.kinetic_energy < energy_threshold
                printstyled("Energy threshold reached! Ending simulation.\n"; color = :green) 
                break
            end

        end 

    end
    
    # show  sections time
    if record_timings
        print_timer()
    end 

    vtk_save(vtkdata.VTKcollection)
    
end

#  Solves the current step
function solve_step_dynamics!(conf, tⁿ⁺¹, Δt, solⁿ, solⁿ⁺¹, nodes_sol, matrices, energy, solver, params) 

    @unpack ext_forces, bcs, sdf = conf

    # -------------------------------------------------------------------------------------------
    # INITIALIZATION
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug "Time step initialization" begin
        
        
        # update external force value 
        for i in ext_forces.loaded_dofs
            solⁿ⁺¹.fᵉˣᵗ[i] = ext_forces.f(tⁿ⁺¹, i)     
        end

        for i in bcs.disp_dofs
            bcs.disp_vals[i] = bcs.u(tⁿ⁺¹, i)     
        end
        
    end
    
    # -------------------------------------------------------------------------------------------
    # PREDICTOR
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug "Predictor" begin
        
        # predict the solution @n+1
        predictor!(conf, matrices, energy, solⁿ⁺¹, solⁿ, Δt, nodes_sol, solver, params)
        
    end
    
    # -------------------------------------------------------------------------------------------
    # CORRECTOR
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug "Corrector" begin
        
        # corrector loop: output = number of iterations
        k = corrector!(conf, matrices, energy, solⁿ⁺¹, solⁿ, Δt, nodes_sol, solver, params)
        
    end
   
    return k
    
end 

# Predicts the solution at the current step
function predictor!(conf, matrices, energy, solⁿ⁺¹, solⁿ, Δt, nodes_sol, solver, params)
    # -------------------------------------------------------------------------------------------
    # INITIALISATION
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug "Initialization" begin

        @unpack nodes, beams, constraints, bcs = conf
        
        # dofs
        free_dofs = bcs.free_dofs

        # parameters for the numerical integration 
        @unpack α, β, γ = params
        
        # solutions @n
        for n in LazyRows(nodes)
            nodes_sol.D[n.idof_disp] .= n.u
            nodes_sol.Ḋ[n.idof_disp] .= n.u̇
            nodes_sol.D̈[n.idof_disp] .= n.ü
            nodes_sol.D[n.idof_rot] .= 0
            nodes_sol.Ḋ[n.idof_rot] .= n.ẇ
            nodes_sol.D̈[n.idof_rot] .= n.ẅ
        end
        
        nodes_sol.ΔD .= 0  

        
    end
    
    # -------------------------------------------------------------------------------------------
    # SOLVE
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug "Assemble element contributions" assemble!(conf, matrices, energy, params)
    
    @timeit_debug "Compute penalty constraints contributions" constraints!(matrices, nodes, constraints)        
    
    @timeit_debug "Compute tangent matrix and residuals" tangent_and_residuals_predictor!(nodes_sol, matrices, solⁿ, solⁿ⁺¹, Δt, α, β, γ)
    
    @timeit_debug "Impose BCs" apply_BCs!(nodes_sol, bcs)

    @timeit_debug "Compute free tangent matrix and residual" begin

        @views nodes_sol.r_free .= nodes_sol.r[free_dofs]              
        @views nodes_sol.Ktan_free.nzval .= nodes_sol.Ktan[matrices.sparsity_free]

    end 
    
    @timeit_debug "Linear solve" begin
        analyze!(solver, nodes_sol.Ktan_free, nodes_sol.r_free)
        solve!(solver, nodes_sol.ΔD_free, nodes_sol.Ktan_free, nodes_sol.r_free)
    end
    
    @timeit_debug "Update" begin 

        # fill whole dofs solution vector with free dofs solution
        nodes_sol.ΔD[free_dofs] .= nodes_sol.ΔD_free

        @inbounds for i in eachindex(nodes_sol.D)
            Ḋⁿ = nodes_sol.Ḋ[i]
            nodes_sol.D[i] += nodes_sol.ΔD[i]
            nodes_sol.Ḋ[i] += (γ/(β*Δt))*nodes_sol.ΔD[i] - (γ/β)*Ḋⁿ + (Δt*(2*β-γ)/(2*β))*nodes_sol.D̈[i]
            nodes_sol.D̈[i] += (1/(β*Δt^2))*(nodes_sol.ΔD[i] - Δt*Ḋⁿ - (Δt^2)/2*nodes_sol.D̈[i])
        end

        # update the local nodes 
        update_local_predictor!(nodes, nodes_sol, Δt, β, γ)
        
    end
    
end 

# Corrects the solution at the current step
function corrector!(conf, matrices, energy, solⁿ⁺¹, solⁿ, Δt, nodes_sol, solver, params) 
    
    # -------------------------------------------------------------------------------------------
    # INITIALISATION
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug "Initialization" begin
        
        # tolerances
        @unpack nodes, beams, constraints, bcs, disp_dofs = conf
        
        # dofs
        free_dofs = bcs.free_dofs

        # parameters for the numerical integration 
        @unpack tol_res, tol_ΔD, max_it, α, β, γ = params
        

        # update @n+1   
        for n in LazyRows(nodes)
            nodes_sol.D[n.idof_disp] .= n.u
            nodes_sol.Ḋ[n.idof_disp] .= n.u̇
            nodes_sol.D̈[n.idof_disp] .= n.ü
            nodes_sol.D[n.idof_rot] .= n.w
            nodes_sol.Ḋ[n.idof_rot] .= n.ẇ
            nodes_sol.D̈[n.idof_rot] .= n.ẅ
        end

        
        # current iteration number, initalised to 1
        k = 1

        # randomly initialised the variables for the while cycle 
        ΔD_norm = Inf
        res_norm = Inf
        
    end
    
    # -------------------------------------------------------------------------------------------
    # CORRECTOR LOOP
    # -------------------------------------------------------------------------------------------

    @timeit_debug "Newton solver" begin
        
        while ( res_norm>tol_res || ΔD_norm>tol_ΔD ) && k≤max_it        
            
            @timeit_debug "Assemble element contributions" assemble!(conf, matrices, energy, params)
            
            @timeit_debug "Compute constraints contributions" constraints!(matrices, nodes, constraints)             
                
            @timeit_debug "Compute tangent matrix and residual" tangent_and_residuals_corrector!(nodes_sol, matrices, solⁿ, solⁿ⁺¹, Δt, α, β, γ)

            @timeit_debug "Impose BCs" apply_BCs!(nodes_sol, bcs)
                
            @timeit_debug "Compute free tangent matrix and residual" begin
                @views nodes_sol.r_free .= nodes_sol.r[free_dofs]               
                @views nodes_sol.Ktan_free.nzval .= nodes_sol.Ktan[matrices.sparsity_free]
            end 
            
            @timeit_debug "Linear solve" solve!(solver, nodes_sol.ΔD_free, nodes_sol.Ktan_free, nodes_sol.r_free)

            
            @timeit_debug "Update global and local variables" begin
    
                # fill whole dofs solution vector with free dofs solution
                nodes_sol.ΔD[free_dofs] .= nodes_sol.ΔD_free   
                
                # update global displacement variables
                @inbounds for i in disp_dofs
                    nodes_sol.D[i] += nodes_sol.ΔD[i]
                    nodes_sol.Ḋ[i] += (γ/(β*Δt))*nodes_sol.ΔD[i]
                    nodes_sol.D̈[i] += (1/(β*Δt^2))*nodes_sol.ΔD[i]
                end

                # update the nodes displacement, velocity, acceleration and rotation
                update_local_corrector!(nodes, nodes_sol, Δt, β, γ)
                
            end 
            
            @timeit_debug "Compute norms" res_norm, ΔD_norm = compute_norms_corrector(k, nodes_sol, params.verbose)

            k += 1

        end
        
        
    end 
    
    release!(solver)
    
    return k
    
end 