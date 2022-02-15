"Preallocates and initialises all the variables and starts the temporal loop"
function solver!(nodes, allbeams, conf, comp, sdf, pn_constrains, params, T=Float64)
    
    # -------------------------------------------------------------------------------------------
    # INITIALISATION
    # -------------------------------------------------------------------------------------------
  
    @unpack scale, ENERGY_STOP, SHOW_COMP_TIME, SHOW_TIME_SECTIONS, SAVE_NODES_VTK, SAVE_ENERGY, SAVE_INTERPOLATION_VTK, SAVE_GP_VTK, thisDirOutputPath = params

    if SHOW_TIME_SECTIONS
        TimerOutputs.enable_debug_timings(@__MODULE__)
    end 

    to = TimerOutput()

    @timeit_debug to "Solver initialization" begin

        # vectors to save energy evolution
        plot_Phi_energy  = T[]
        plot_K_energy = T[]
        plot_contact_energy = T[]
        
        # initialization of time loop variables
        t = zero(T)
        start = time()
        i = 1
        Δt = comp.Δt
        Δt_plot = comp.Δt_plot 
        fact_div = 1   
        k_count_OKit = 0
        step = 1

        # preallocated vectors for beam interpolation
        int_pos = zeros(T, (allbeams.numberInterpolationPoints[1]+1)*length(allbeams), 3)
        int_conn = zeros(Int, allbeams.numberInterpolationPoints[1]+1, length(allbeams))
        
        # preallocation of the vectors and matrices used in the computaion
        solⁿ, solⁿ⁺¹, sol_GP, nodes_sol, matrices, energy, u_imp = solver_initialisation(conf, nodes, allbeams, int_pos, int_conn, comp, pn_constrains, thisDirOutputPath, SAVE_INTERPOLATION_VTK, SAVE_GP_VTK, SAVE_NODES_VTK, T)
        
        # Linear solver
        ps = MKLPardisoSolver()
        
    end 
    
    # -------------------------------------------------------------------------------------------
    # TIME LOOP 
    # -------------------------------------------------------------------------------------------

    STOP_SIMULATION = false

    @timeit_debug to "Solver time loop" begin
        
        while (t<comp.tᵉⁿᵈ && fact_div < 1E20 && !STOP_SIMULATION)
            
            @timeit_debug to "Dynamic step" begin
                
                tⁿ = t
                tⁿ⁺¹ = tⁿ + Δt
                tcomp = time() - start
                
                if SHOW_COMP_TIME
                    println("-------------------------------------------------")
                    println("tphys = $tⁿ⁺¹, tcomp = $tcomp, ratio = $(tcomp/t), Δt = $Δt, step = $step ")
                end 
                
                # solve system @tⁿ
                flag_conv = solve_step_dynamics!(nodes, allbeams, pn_constrains, tⁿ⁺¹, Δt, solⁿ, solⁿ⁺¹, sol_GP, nodes_sol, matrices, energy, u_imp, conf, comp, sdf, to, ps, SHOW_COMP_TIME, T)

                # if not converged, halve the time step and re-solve the system until it converges
                while (flag_conv == 0 && fact_div < 1E20)
                    
                    fact_div = fact_div*scale
                    Δt = Δt/scale
                    
                    tⁿ⁺¹ = tⁿ + Δt
                    
                    if SHOW_COMP_TIME
                        println("-----TIME STEP DECREASED Δt = Δt0/$fact_div")
                        println("tphys = $tⁿ⁺¹, tcomp = $tcomp  Δt = $Δt, ratio = $(tcomp/t)")
                    end

                    if Δt<1E-10
                        @error "Δt < 1E-10: stopping simulation" 
                        STOP_SIMULATION = true
                        break
                    end
                    
                    solⁿ⁺¹ = deepcopy(solⁿ)
                    update_nodes_not_converged!(nodes)
                    
                    flag_conv = solve_step_dynamics!(nodes, allbeams, pn_constrains, tⁿ⁺¹, Δt, solⁿ, solⁿ⁺¹, sol_GP, nodes_sol, matrices, energy, u_imp, conf, comp, sdf, to, ps, SHOW_COMP_TIME), T
                    k_count_OKit = 0
                    
                end
                
                # if the time step has been halved but the simulation has converged the last 8 times, double it
                #if fact_div > 1
                    
                    k_count_OKit = k_count_OKit + 1
                
                    if k_count_OKit > 15
                        
                        if SHOW_COMP_TIME
                            println("-----TIME STEP INCREASED Δt=Δt0/$fact_div")
                        end 

                        fact_div = fact_div/scale
                        Δt = Δt*scale
                        
                        k_count_OKit = 0
                        
                    end
                    
                #end
                
            end
            
            @timeit_debug to "Write results to file" begin
                
                # save VTK with frequency 1/Δt_plot: 
                if  abs(tⁿ⁺¹-i*Δt_plot)<1e-9 || tⁿ⁺¹>i*Δt_plot
                    save_VTK(i, nodes, allbeams, sol_GP, int_pos, int_conn, thisDirOutputPath, SAVE_INTERPOLATION_VTK, SAVE_GP_VTK, SAVE_NODES_VTK)
                    i = i+1
                end

                # save the energy values at the current step
                push!(plot_Phi_energy, energy.strain_energy)
                push!(plot_K_energy, energy.kinetic_energy)
                push!(plot_contact_energy, energy.contact_energy)
                
                if SAVE_ENERGY
                    save_energy(plot_Phi_energy, plot_K_energy, plot_contact_energy, thisDirOutputPath)           
                end        
                
            end 
            
            @timeit_debug to "Update variables" begin
                
                # update variables after the step
                solⁿ = deepcopy(solⁿ⁺¹)
                update_nodes_converged!(nodes)
                t = tⁿ⁺¹
                step = step+1
                
            end 

            if ENERGY_STOP
                if plot_K_energy[end] < 1e-10
                    println("ENERGY THRESHOLD REACHED")
                    STOP_SIMULATION = true
                end 
            end

        end 

    end
    
    # show  sections time
    if SHOW_TIME_SECTIONS
        show(to)
    end 
    
end

#  Solves the current step
function solve_step_dynamics!(nodes, allbeams, pn_constrains, tⁿ⁺¹, Δt, solⁿ, solⁿ⁺¹, sol_GP, nodes_sol, matrices, energy, u_imp, conf, comp, sdf, to, ps, SHOW_COMP_TIME, T=Float64) 

    # -------------------------------------------------------------------------------------------
    # INITIALIZATION
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug to "Time step initialization" begin
        
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
    
    @timeit_debug to "Predictor" begin
        
        # predict the solution @n+1
        predictor!(nodes, allbeams, pn_constrains, matrices, energy, solⁿ⁺¹, solⁿ, sol_GP, u_imp, Δt, conf, comp, sdf, nodes_sol, tⁿ⁺¹, to, ps, T) 
        
    end
    
    # -------------------------------------------------------------------------------------------
    # CORRECTOR
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug to "Corrector" begin
        
        # corrector loop: output = number of iterations
        k = corrector_loop!(nodes, allbeams, pn_constrains, matrices, energy, solⁿ⁺¹, solⁿ, sol_GP, u_imp, Δt, conf, comp, sdf, nodes_sol, tⁿ⁺¹, to, ps, SHOW_COMP_TIME, T)
        
    end
    
    return k < max_it
    
end 

# Predicts the solution at the current step
function predictor!(nodes, allbeams, pencons, matrices, energy, solⁿ⁺¹, solⁿ, sol_GP, u_imp, Δt, conf, comp, sdf, nodes_sol, tⁿ⁺¹, to, ps, T=Float64)
    
    # -------------------------------------------------------------------------------------------
    # INITIALISATION
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug to "Initialization" begin
        
        # dofs
        free_dofs = conf.bc.free_dofs
        fixed_dofs = conf.bc.fixed_dofs

        # parameters for the numerical integration 
        β = comp.β
        γ = comp.γ
        α = comp.α
        
        # solutions @n
        update_global_predictor!(nodes_sol, nodes)
        
    end
    
    # -------------------------------------------------------------------------------------------
    # SOLVE
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug to "Compute matrices" begin
        
        # compute the matrices    
        compute_K_T!(nodes, allbeams, matrices, energy, conf, sdf, comp, sol_GP, to, T)
        
    end
    
    @timeit_debug to "Compute penalty constrains contributions" begin
        
        # impose multi freedom constrains
        compute_multifreedom_constraints!(matrices, nodes, pencons, tⁿ⁺¹)        
    end 
    
    @timeit_debug to "Compute tangent matrix" begin
        
        # eq142 in [2]: compute tangent matrix
        compute_Ktan_sparse!(nodes_sol, matrices, α, β, γ, Δt)

    end
    
    @timeit_debug to "Compute residual" begin

        mul!(nodes_sol.r, nodes_sol.Ktan, nodes_sol.ΔD)


    end
    
    @timeit_debug to "Impose BCs" begin
        
        # impose BCs   
        if conf.bc.flag_cylindrical == 1       
            dirichlet_global_to_local!(nodes_sol, nodes)     
        end  
        
        # impose single freedom constrain
        impose_BC_displacements!(nodes_sol, u_imp, fixed_dofs)
    end 


    @timeit_debug to "Compute free tangent matrix and residual" begin

        fill_r_free!(nodes_sol, free_dofs)                
        fill_Ktan_free!(nodes_sol)

    end 
    
    @timeit_debug to "Linear solve" begin 

        solve!(ps, nodes_sol.ΔD_free, nodes_sol.Ktan_free, nodes_sol.r_free)

    end
    
    @timeit_debug to "Update" begin 

        # fill whole dofs solution vector with free dofs solution
        fill_ΔD_free_dofs!(nodes_sol, free_dofs)
                        
        if conf.bc.flag_cylindrical == 1
            dirichlet_local_to_global!(nodes_sol, nodes)    
        end  
        
        # eq114-115 in [3]
        update_nodal_solutions_predictor!(nodes_sol, β, γ, Δt)
        
        # update the local nodes 
        update_local_predictor!(nodes, nodes_sol)
        
    end
    
end 

# Corrects the solution at the current step
function corrector_loop!(nodes, allbeams, pncons, matrices, energy, solⁿ⁺¹, solⁿ, sol_GP, u_imp, Δt, conf, comp, sdf, nodes_sol, tⁿ⁺¹, to, ps, SHOW_COMP_TIME, T=Float64) 
    
    # -------------------------------------------------------------------------------------------
    # INITIALISATION
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug to "Initialization" begin
        
        # tollerances
        tol_res = comp.tol_res
        tol_ΔD_k = comp.tol_ΔDk
        max_it = comp.max_it
        
        # dofs
        ndofs = conf.ndofs
        free_dofs = conf.bc.free_dofs
        fixed_dofs = conf.bc.fixed_dofs
        disp_dofs = conf.disp_dofs

        # update @n+1   
        update_global_corrector!(nodes_sol, nodes, disp_dofs)
        
        # constants for the numerical integration
        β = comp.β
        γ = comp.γ
        α = comp.α
        
        # current iteration number, initalised to 0
        k = 0
        
        # randomly initialised the variables for the while cycle 
        nodes_sol.ΔD .= 5 
        aux_tol = 1.0e10
        
    end
    
    # -------------------------------------------------------------------------------------------
    # CORRECTOR LOOP
    # -------------------------------------------------------------------------------------------
    
    @timeit_debug to "Newton solver" begin
        
        while ((aux_tol>tol_res) || (((norm(nodes_sol.ΔD[free_dofs])))>tol_ΔD_k)) && (k<max_it)         

            k += 1
            
            @timeit_debug to "Compute matrices" begin
                
                # compute the contributions for each beam and assemble into the Matrix structure
                compute_K_T!(nodes, allbeams, matrices, energy, conf, sdf, comp, sol_GP, to, T)  
                
            end
            
            @timeit_debug to "Compute multi freedom constrains contributions" begin
                
                # impose multi freedom constrains
                compute_multifreedom_constraints!(matrices, nodes, pncons, tⁿ⁺¹)             
                
            end 
            
            @timeit_debug to "Compute tangent matrix and residual" begin
                
                # eq142 in [2]: compute tangent matrix
                compute_Ktan_sparse!(nodes_sol, matrices, α, β, γ, Δt)

                # eq140 in [2]: compute residual 
                compute_res_corrector!(nodes_sol, matrices, solⁿ⁺¹, solⁿ, α)

            end
            
            @timeit_debug to "Update nodal solutions" begin
                
                # update global nodal solution at each iteration
                update_nodal_solution_corrector_loop!(nodes_sol, disp_dofs)
                
            end 
            
            @timeit_debug to "Impose BCs" begin
                
                
                # if necessary, move to cylindrical coordinates     
                if conf.bc.flag_cylindrical == 1  
                    dirichlet_global_to_local!(nodes_sol, nodes)          
                end     
                
                # impose single freedom constrains
                impose_BC_displacements!(nodes_sol, u_imp, fixed_dofs)
                
            end 
            
            @timeit_debug to "Compute free tangent matrix and residual" begin

                fill_r_free!(nodes_sol, free_dofs)                
                fill_Ktan_free!(nodes_sol)
        
            end 
            
            @timeit_debug to "Linear solve" begin
                
                solve!(ps, nodes_sol.ΔD_free, nodes_sol.Ktan_free, nodes_sol.r_free)

            end 
            
            @timeit_debug to "Update global and local variables" begin
    
                # fill whole dofs solution vector with free dofs solution
                fill_ΔD_free_dofs!(nodes_sol, free_dofs)
                
                # if necessary, return to cartesian coordinates        
                if conf.bc.flag_cylindrical == 1   
                    dirichlet_local_to_global!(nodes_sol,  nodes)            
                end     
                
                # update global displacement variables
                update_nodal_solutions_corrector!(nodes_sol, disp_dofs, γ, β, Δt)
                
                # update the nodes displacement, velocity, acceleration and rotation
                update_local_corrector!(nodes, nodes_sol.ΔD, Δt, nodes_sol, comp)
                
                # update the current solution vectors (used in the next time step)
                update_current_solution_corrector!(solⁿ⁺¹, ndofs, matrices)
            end 
            
            @timeit_debug to "Compute norms" begin
                
                # update norms (used to check while cycle)
                aux_tol = compute_norms_corrector(k, aux_tol, solⁿ⁺¹, nodes_sol, matrices, SHOW_COMP_TIME)
                
            end 

        end
        
        
    end 
    
    # -------------------------------------------------------------------------------------------
    # OUTPUT
    # -------------------------------------------------------------------------------------------
    
    return k 
    
end 