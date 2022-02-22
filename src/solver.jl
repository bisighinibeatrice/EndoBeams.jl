"Preallocates and initialises all the variables and starts the temporal loop"
function solver!(nodes, beams, conf, comp, sdf, constraints, params, T=Float64)
    
    # -------------------------------------------------------------------------------------------
    # INITIALISATION
    # -------------------------------------------------------------------------------------------
  
    @unpack scale, ENERGY_STOP, SHOW_COMP_TIME, SHOW_TIME_SECTIONS, SAVE_NODES_VTK, SAVE_ENERGY, SAVE_INTERPOLATION_VTK, SAVE_GP_VTK, thisDirOutputPath = params

    if SHOW_TIME_SECTIONS
        TimerOutputs.enable_debug_timings(@__MODULE__)
        reset_timer!()
    end 

    @timeit_debug "Solver initialization" begin

        # vectors to save energy evolution
        plot_Phi_energy  = T[]
        plot_K_energy = T[]
        plot_contact_energy = T[]
        
        # initialization of time loop variables
        t = 0.
        start = time()
        write_counter = 1
        Δt = comp.Δt
        Δt_plot = comp.Δt_plot 
        fact_div = 1.   
        successive_convergence = 0
        step = 1

        # preallocated vectors for beam interpolation
        int_pos = zeros(T, (beams.numberInterpolationPoints[1]+1)*length(beams), 3)
        int_conn = zeros(Int, beams.numberInterpolationPoints[1]+1, length(beams))
        
        # preallocation of the vectors and matrices used in the computaion
        solⁿ, solⁿ⁺¹, sol_GP, nodes_sol, matrices, energy, u_imp = solver_initialisation(conf, nodes, beams, int_pos, int_conn, comp, constraints, thisDirOutputPath, SAVE_INTERPOLATION_VTK, SAVE_GP_VTK, SAVE_NODES_VTK, T)
        
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
                flag_conv = solve_step_dynamics!(nodes, beams, constraints, tⁿ⁺¹, Δt, solⁿ, solⁿ⁺¹, sol_GP, nodes_sol, matrices, energy, u_imp, conf, comp, sdf, ps, SHOW_COMP_TIME, T)

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
                    save_VTK(write_counter, nodes, beams, sol_GP, int_pos, int_conn, thisDirOutputPath, SAVE_INTERPOLATION_VTK, SAVE_GP_VTK, SAVE_NODES_VTK)
                    write_counter += 1
                end

                # save the energy values at the current step
                push!(plot_Phi_energy, energy.strain_energy)
                push!(plot_K_energy, energy.kinetic_energy)
                push!(plot_contact_energy, energy.contact_energy)
                
                if SAVE_ENERGY
                    save_energy(plot_Phi_energy, plot_K_energy, plot_contact_energy, thisDirOutputPath)           
                end        
                
            end 
            
            @timeit_debug "Update variables" begin
                
                # update variables after the step
                solⁿ = deepcopy(solⁿ⁺¹)
                update_nodes_converged!(nodes, solⁿ, matrices)
                t = tⁿ⁺¹
                step += 1
                
            end 

            if ENERGY_STOP && plot_K_energy[end] < 1e-6
                printstyled("Energy threshold reached! Ending simulation.\n"; color = :green) 
                break
            end

        end 

    end
    
    # show  sections time
    if SHOW_TIME_SECTIONS
        print_timer()
    end 
    
end

#  Solves the current step
function solve_step_dynamics!(nodes, beams, constraints, tⁿ⁺¹, Δt, solⁿ, solⁿ⁺¹, sol_GP, nodes_sol, matrices, energy, u_imp, conf, comp, sdf, ps, SHOW_COMP_TIME, T=Float64) 

    # -------------------------------------------------------------------------------------------
    # INITIALIZATION
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug "Time step initialization" begin
        
        # maximum number of iterations
        max_it = comp.max_it
        
        # update external force value
        if conf.ext_forces.flag_crimping == false
            update_current_solution_external_force!(solⁿ⁺¹, tⁿ⁺¹, conf)
        else 
            update_current_solution_external_force_crimping!(solⁿ⁺¹, tⁿ⁺¹, conf, nodes)
        end 
        
        # update boundary conditions value
        if conf.bc.flag_disp_vector == false 
            update_current_boundary_conditions!(u_imp, tⁿ⁺¹, conf)
        else
            update_current_boundary_conditions_vector!(u_imp, tⁿ⁺¹, conf)
        end 
        
    end
    
    # -------------------------------------------------------------------------------------------
    # PREDICTOR
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug "Predictor" begin
        
        # predict the solution @n+1
        predictor!(nodes, beams, constraints, matrices, energy, solⁿ⁺¹, solⁿ, sol_GP, u_imp, Δt, conf, comp, sdf, nodes_sol, tⁿ⁺¹, ps, T) 
        
    end
    
    # -------------------------------------------------------------------------------------------
    # CORRECTOR
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug "Corrector" begin
        
        # corrector loop: output = number of iterations
        k = corrector_loop!(nodes, beams, constraints, matrices, energy, solⁿ⁺¹, solⁿ, sol_GP, u_imp, Δt, conf, comp, sdf, nodes_sol, tⁿ⁺¹, ps, SHOW_COMP_TIME, T)
        
    end
   
    return k ≤ max_it
    
end 

# Predicts the solution at the current step
function predictor!(nodes, beams, pencons, matrices, energy, solⁿ⁺¹, solⁿ, sol_GP, u_imp, Δt, conf, comp, sdf, nodes_sol, tⁿ⁺¹, ps, T=Float64)
    
    # -------------------------------------------------------------------------------------------
    # INITIALISATION
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug "Initialization" begin
        
        # dofs
        free_dofs = conf.bc.free_dofs
        fixed_dofs = conf.bc.fixed_dofs

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
    
    @timeit_debug "Assemble element contributions" assemble!(nodes, beams, matrices, energy, conf, sdf, comp, sol_GP)
    
    @timeit_debug "Compute penalty constrains contributions" constraints!(matrices, nodes, pencons)        
    
    @timeit_debug "Compute tangent matrix and residuals" tangent_and_residuals_predictor!(nodes_sol, matrices, solⁿ, solⁿ⁺¹, Δt, α, β, γ)
    
    @timeit_debug "Impose BCs" apply_BCs!(nodes_sol, u_imp, fixed_dofs)

    @timeit_debug "Compute free tangent matrix and residual" begin

        @views nodes_sol.r_free .= nodes_sol.r[free_dofs]              
        @views nodes_sol.Ktan_free.nzval .= nodes_sol.Ktan[nodes_sol.sparsity_map_free]

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
function corrector_loop!(nodes, beams, pncons, matrices, energy, solⁿ⁺¹, solⁿ, sol_GP, u_imp, Δt, conf, comp, sdf, nodes_sol, tⁿ⁺¹, ps, SHOW_COMP_TIME, T=Float64) 
    
    # -------------------------------------------------------------------------------------------
    # INITIALISATION
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug "Initialization" begin
        
        # tollerances
        res_tol = comp.res_tol
        ΔD_tol = comp.tol_ΔDk
        max_it = comp.max_it
        
        # dofs
        free_dofs = conf.bc.free_dofs
        fixed_dofs = conf.bc.fixed_dofs
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
        
        while ( res_norm>res_tol || ΔD_norm>ΔD_tol ) && k≤max_it        
            
            @timeit_debug "Assemble element contributions" assemble!(nodes, beams, matrices, energy, conf, sdf, comp, sol_GP)  
            
            @timeit_debug "Compute constraints contributions" constraints!(matrices, nodes, pncons)             
                
            @timeit_debug "Compute tangent matrix and residual" tangent_and_residuals_corrector!(nodes_sol, matrices, solⁿ, solⁿ⁺¹, Δt, α, β, γ)

            @timeit_debug "Impose BCs" apply_BCs!(nodes_sol, u_imp, fixed_dofs)
                
            @timeit_debug "Compute free tangent matrix and residual" begin
                @views nodes_sol.r_free .= nodes_sol.r[free_dofs]               
                @views nodes_sol.Ktan_free.nzval .= nodes_sol.Ktan[nodes_sol.sparsity_map_free]
            end 
            
            @timeit_debug "Linear solve" solve!(ps, nodes_sol.ΔD_free, nodes_sol.Ktan_free, nodes_sol.r_free)
            
            @timeit_debug "Update global and local variables" begin
    
                # fill whole dofs solution vector with free dofs solution
                nodes_sol.ΔD[free_dofs] .= nodes_sol.ΔD_free  
                
                # if necessary, return to cartesian coordinates        
                # if conf.bc.flag_cylindrical == 1   
                #     dirichlet_local_to_global!(nodes_sol,  nodes)            
                # end     
                
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