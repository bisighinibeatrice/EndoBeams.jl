"Preallocates and initialises all the variables and starts the temporal loop"
function solver!(nodes, beams, conf, comp, sdf, constraints, params, T=Float64)
    
    # -------------------------------------------------------------------------------------------
    # INITIALISATION
    # -------------------------------------------------------------------------------------------
  
    @unpack scale, ENERGY_STOP, SHOW_COMP_TIME, SHOW_TIME_SECTIONS, output_dir = params

    if SHOW_TIME_SECTIONS
        TimerOutputs.enable_debug_timings(@__MODULE__)
        reset_timer!()
    end 

    @timeit_debug "Solver initialization" begin
        
        # initialization of time loop variables
        t = 0.
        start = time()
        write_counter = 0
        Δt = comp.Δt
        Δt_plot = comp.Δt_plot 
        fact_div = 1.   
        successive_convergence = 0
        step = 1
        
        # preallocation of the vectors and matrices used in the computaion
        solⁿ, solⁿ⁺¹, nodes_sol, matrices, energy, vtkdata = solver_initialisation(conf, nodes, beams, constraints, sdf, output_dir, T)

        write_VTK(write_counter, 0, t, nodes, beams, energy, conf, sdf, comp, vtkdata)
        write_counter += 1
        
        # Linear solver
        ps = MKLPardisoSolver()
        
    end 
    
    # -------------------------------------------------------------------------------------------
    # TIME LOOP 
    # -------------------------------------------------------------------------------------------

    STOP_SIMULATION = false

    @timeit_debug "Solver time loop" begin
        
        while (t<comp.tᵉⁿᵈ && fact_div < 1E20 && !STOP_SIMULATION)
            
            @timeit_debug "Dynamic step" begin
                
                tⁿ = t
                tⁿ⁺¹ = tⁿ + Δt
                tcomp = time() - start
                
                if SHOW_COMP_TIME 
                    println("------------------------------------------------")
                    @printf "%8s\t%8s\t%4s\t%4s\n" "t" "Δt" "step" "tcomp"
                    @printf "%1.2e\t%1.2e\t%4d\t%1.2f\n" tⁿ⁺¹ Δt step tcomp
                    println("------------------------------------------------")
                end 
                
                # solve system @tⁿ
                flag_conv = solve_step_dynamics!(nodes, beams, constraints, tⁿ⁺¹, Δt, solⁿ, solⁿ⁺¹, nodes_sol, matrices, energy, conf, comp, sdf, ps, SHOW_COMP_TIME)

                # if not converged, halve the time step and re-solve the system until it converges
                if !flag_conv

                    SHOW_COMP_TIME && printstyled("Time step did not converge...\n"; color = :yellow)
                    
                    fact_div = fact_div*scale
                    Δt = Δt/scale
                    
                    SHOW_COMP_TIME && printstyled("Δt decreased: Δt₀/$fact_div\n"; color = :yellow)

                    if Δt<1E-10
                        printstyled("Δt < 1E-10: aborting simulation...\n"; color = :red)
                        break
                    end
                    
                    solⁿ⁺¹ = deepcopy(solⁿ)
                    update_nodes_not_converged!(nodes)
                    successive_convergence = 0

                    continue

                else

                    SHOW_COMP_TIME && printstyled("Time step converged!\n"; color = :green) 
                    successive_convergence = successive_convergence + 1

                end
                
                # if the time step has been halved but the simulation has converged the last 8 times, double it
                if fact_div > 1 && successive_convergence > 8
                        
                        SHOW_COMP_TIME && printstyled("Δt increased: Δt₀/$fact_div\n"; color = :green)

                        fact_div = fact_div/scale
                        Δt = Δt*scale
                        
                        successive_convergence = 0
                    
                end
                
            end
            
            @timeit_debug "Write results to file" begin
                
                # save VTK with frequency 1/Δt_plot: 
                if  tⁿ⁺¹ > write_counter*Δt_plot  ||  tⁿ⁺¹ ≈ write_counter*Δt_plot
                    write_VTK(write_counter, step, tⁿ⁺¹, nodes, beams, energy, conf, sdf, comp, vtkdata)
                    write_counter += 1
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

            if ENERGY_STOP && energy.kinetic_energy < 1e-6
                printstyled("Energy threshold reached! Ending simulation.\n"; color = :green) 
                break
            end

        end 

    end
    
    # show  sections time
    if SHOW_TIME_SECTIONS
        print_timer()
    end 

    vtk_save(vtkdata.VTKcollection)
    
end

#  Solves the current step
function solve_step_dynamics!(nodes, beams, constraints, tⁿ⁺¹, Δt, solⁿ, solⁿ⁺¹, nodes_sol, matrices, energy, conf, comp, sdf, ps, SHOW_COMP_TIME) 

    # -------------------------------------------------------------------------------------------
    # INITIALIZATION
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug "Time step initialization" begin
        
        # maximum number of iterations
        max_it = comp.max_it
        
        # update external force value 
        for i in conf.ext_forces.loaded_dofs
            solⁿ⁺¹.fᵉˣᵗ[i] = conf.ext_forces.f(tⁿ⁺¹, i)     
        end

        for i in conf.bcs.disp_dofs
            conf.bcs.disp_vals[i] = conf.bcs.u(tⁿ⁺¹, i)     
        end
        
        
    end
    
    # -------------------------------------------------------------------------------------------
    # PREDICTOR
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug "Predictor" begin
        
        # predict the solution @n+1
        predictor!(nodes, beams, constraints, matrices, energy, solⁿ⁺¹, solⁿ, Δt, conf, comp, sdf, nodes_sol, ps) 
        
    end
    
    # -------------------------------------------------------------------------------------------
    # CORRECTOR
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug "Corrector" begin
        
        # corrector loop: output = number of iterations
        k = corrector_loop!(nodes, beams, constraints, matrices, energy, solⁿ⁺¹, solⁿ, Δt, conf, comp, sdf, nodes_sol, ps, SHOW_COMP_TIME)
        
    end
   
    return k ≤ max_it
    
end 

# Predicts the solution at the current step
function predictor!(nodes, beams, pencons, matrices, energy, solⁿ⁺¹, solⁿ, Δt, conf, comp, sdf, nodes_sol, ps)
    
    # -------------------------------------------------------------------------------------------
    # INITIALISATION
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug "Initialization" begin
        
        # dofs
        free_dofs = conf.bcs.free_dofs
        fixed_dofs = conf.bcs.fixed_dofs

        # parameters for the numerical integration 
        β = comp.β
        γ = comp.γ
        α = comp.α
        
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
    
    @timeit_debug "Assemble element contributions" assemble!(nodes, beams, matrices, energy, conf, sdf, comp)
    
    @timeit_debug "Compute penalty constraints contributions" constraints!(matrices, nodes, pencons)        
    
    @timeit_debug "Compute tangent matrix and residuals" tangent_and_residuals_predictor!(nodes_sol, matrices, solⁿ, solⁿ⁺¹, Δt, α, β, γ)
    
    @timeit_debug "Impose BCs" apply_BCs!(nodes_sol, conf.bcs)

    @timeit_debug "Compute free tangent matrix and residual" begin

        @views nodes_sol.r_free .= nodes_sol.r[free_dofs]              
        @views nodes_sol.Ktan_free.nzval .= nodes_sol.Ktan[matrices.sparsity_free]

    end 
    
    @timeit_debug "Linear solve" solve!(ps, nodes_sol.ΔD_free, nodes_sol.Ktan_free, nodes_sol.r_free)
    
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
function corrector_loop!(nodes, beams, constraints, matrices, energy, solⁿ⁺¹, solⁿ, Δt, conf, comp, sdf, nodes_sol, ps, SHOW_COMP_TIME) 
    
    # -------------------------------------------------------------------------------------------
    # INITIALISATION
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug "Initialization" begin
        
        # tollerances
        tol_res = comp.tol_res
        ΔD_tol = comp.tol_ΔD
        max_it = comp.max_it
        
        # dofs
        free_dofs = conf.bcs.free_dofs
        fixed_dofs = conf.bcs.fixed_dofs
        disp_dofs = conf.disp_dofs

        # update @n+1   
        for n in LazyRows(nodes)
            nodes_sol.D[n.idof_disp] .= n.u
            nodes_sol.Ḋ[n.idof_disp] .= n.u̇
            nodes_sol.D̈[n.idof_disp] .= n.ü
            nodes_sol.D[n.idof_rot] .= n.w
            nodes_sol.Ḋ[n.idof_rot] .= n.ẇ
            nodes_sol.D̈[n.idof_rot] .= n.ẅ
        end

        
        # constants for the numerical integration
        β = comp.β
        γ = comp.γ
        α = comp.α

        
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
        
        while ( res_norm>tol_res || ΔD_norm>ΔD_tol ) && k≤max_it        
            
            @timeit_debug "Assemble element contributions" assemble!(nodes, beams, matrices, energy, conf, sdf, comp)
            
            @timeit_debug "Compute constraints contributions" constraints!(matrices, nodes, constraints)             
                
            @timeit_debug "Compute tangent matrix and residual" tangent_and_residuals_corrector!(nodes_sol, matrices, solⁿ, solⁿ⁺¹, Δt, α, β, γ)

            @timeit_debug "Impose BCs" apply_BCs!(nodes_sol, conf.bcs)
                
            @timeit_debug "Compute free tangent matrix and residual" begin
                @views nodes_sol.r_free .= nodes_sol.r[free_dofs]               
                @views nodes_sol.Ktan_free.nzval .= nodes_sol.Ktan[matrices.sparsity_free]
            end 
            
            @timeit_debug "Linear solve" solve!(ps, nodes_sol.ΔD_free, nodes_sol.Ktan_free, nodes_sol.r_free)
            
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
            
            @timeit_debug "Compute norms" res_norm, ΔD_norm = compute_norms_corrector(k, solⁿ⁺¹, nodes_sol, matrices, SHOW_COMP_TIME)

            k += 1

        end
        
        
    end 
    
    # -------------------------------------------------------------------------------------------
    # OUTPUT
    # -------------------------------------------------------------------------------------------
    
    return k
    
end 