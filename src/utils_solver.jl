#----------------------------------
# FUNCTIONS USED IN THE TIME LOOP 
#----------------------------------

# Cleans the output folders from files of the precedent computation
function clean_folders(thisDirOutputPath)

    if thisDirOutputPath != ""
        dir = pwd()
        cd(thisDirOutputPath)
        foreach(rm, filter(endswith(".vtk"), readdir()))
        cd(dir)
    end 
    
end 

# Calls the functions saving the VTKs related to the nodes and beams positions i-snapshot
function save_VTK(i, allnodes, allbeams, sol_GP, int_pos, int_conn, dirOutput, SAVE_INTERPOLATION_VTK = false, SAVE_GP_VTK = false, SAVE_NODES_VTK = false)
    
    if SAVE_INTERPOLATION_VTK    
        write_VTK_beams(i, allnodes, allbeams, int_pos, int_conn, dirOutput)
    end 
    
    if SAVE_GP_VTK 
        write_VTK_GP(i, sol_GP, dirOutput) 
    end 
    
    if SAVE_NODES_VTK 
        write_VTK_nodes(i, allnodes, allbeams, dirOutput)
    end 

end 

# Cleans folders, pre-allocate and initialise the variables used during the simulation and save the VTKs of the initial configuration
function solver_initialisation(conf, allnodes, allbeams, int_pos, int_conn, comp, cons, thisDirOutputPath, SAVE_INTERPOLATION_VTK = false, SAVE_GP_VTK = false, SAVE_NODES_VTK = false, T=Float64)

    clean_folders(thisDirOutputPath)

    sol_n = constructor_solution(conf, T)
    sol_n1 = constructor_solution(conf, T)
    sol_GP = constructor_solution_GP(length(allbeams), T)
    energy = constructor_energy(T) 
    fixed_matrices = constructor_preallocated_matrices_fixed(allbeams, comp, T)
    uimp = zeros(length(conf.bc.fixed_dofs))
    matrices, nodes_sol = constructor_sparse_matrices!(allbeams, allnodes, cons, conf, T)

    save_VTK(0, allnodes, allbeams, sol_GP, int_pos, int_conn, thisDirOutputPath, SAVE_INTERPOLATION_VTK,  SAVE_GP_VTK, SAVE_NODES_VTK)

    return sol_n, sol_n1, sol_GP, nodes_sol, matrices, energy, fixed_matrices, uimp

end 

# Save energy evolution in text file
function save_energy(plot_Phi_energy, plot_K_energy, plot_C_energy, outputDir)

    open(outputDir * "/PHI ENERGY.txt", "w") do io
        writedlm(io, plot_Phi_energy)
    end

    open(outputDir * "/K ENERGY.txt", "w") do io
        writedlm(io, plot_K_energy)
    end

    open(outputDir * "/C ENERGY.txt", "w") do io
        writedlm(io, plot_C_energy)
    end

end

# Update nodes values @n if converged
function  update_nodes_converged!(allnodes)

    @inbounds for i in 1:length(allnodes)      

        allnodes.u_n[i] = deepcopy(allnodes.u[i])
        allnodes.udt_n[i] = deepcopy(allnodes.udt[i])
        allnodes.udtdt_n[i] = deepcopy(allnodes.udtdt[i])
        allnodes.w_n[i] = deepcopy(allnodes.w[i])
        allnodes.wdt_n[i] = deepcopy(allnodes.wdt[i])
        allnodes.wdtdt_n[i] = deepcopy(allnodes.wdtdt[i])
        allnodes.R_n[i] = deepcopy(allnodes.R[i])
        allnodes.Delt_n[i] = deepcopy(allnodes.Delt[i])

    end

end

# Update nodes values @n+1 if NOT converged
function update_nodes_not_converged!(allnodes)

    @inbounds for i in 1:length(allnodes)

        allnodes.u[i] = deepcopy(allnodes.u_n[i])
        allnodes.udt[i] = deepcopy(allnodes.udt_n[i])
        allnodes.udtdt[i] = deepcopy(allnodes.udtdt_n[i])
        allnodes.w[i] = deepcopy(allnodes.w_n[i])
        allnodes.wdt[i] = deepcopy(allnodes.wdt_n[i])
        allnodes.wdtdt[i] = deepcopy(allnodes.wdtdt_n[i])
        allnodes.R[i] = deepcopy(allnodes.R_n[i])
        allnodes.Delt[i] = deepcopy(allnodes.Delt_n[i])

    end

end

#----------------------------------------
# FUNCTIONS TO UPDATE THE EXTERNAL LOADS
#----------------------------------------

# Replaces fext with the external force @t
function get_current_external_force!(fext, t, conf)

    for i in conf.ext_forces.dof_load
        fext[i] = conf.ext_forces.Fext(t)     
    end

end 

# Replaces fext with the external force @t for the radial crimping
function get_current_external_force_crimping!(fext, t, conf, allnodes)
    
    update_local_to_global_matrix!(allnodes)
    
    for n in allnodes

        dof_disp_n = n.idof_disp
        fext_n =  [conf.ext_forces.Fext(t), 0, 0] 
        fext_n = (n.R_global_to_local)' * fext_n
        fext[dof_disp_n] .= fext_n   

    end 
   
end 

# Update the solution external force @t
function update_current_solution_external_force!(sol_n1, t, conf)

    get_current_external_force!(sol_n1.fext, t, conf)

end 

# Update the solution external force @t for the radial crimping
function update_current_solution_external_force_crimping!(sol_n1, t, conf, allnodes)

    get_current_external_force_crimping!(sol_n1.fext, t, conf, allnodes)
    
end 

#----------------------------------------
# FUNCTIONS TO UPDATE AND IMPOSED THE BCs
#----------------------------------------

# Convert the tangent matrix and the residual from carthesian to cylindrical coordinates
function dirichlet_global_to_local!(nodes_sol, allnodes)
    
    @inbounds for a = 1:length(allnodes)
        
        Ra = allnodes.R_global_to_local[a]    
        RaT = Ra'
        
        adof = 6*(a-1) .+ Vec3(1,2,3)

        nodes_sol.r[adof] .= Ra*nodes_sol.r[adof]
        nodes_sol.asol[adof] .= Ra*nodes_sol.asol[adof]
        
        @inbounds for b = 1:length(allnodes)
            
            bdof = 6*(b-1) .+ Vec3(1,2,3)

            Kab_loc = Mat33( 
            nodes_sol.Ktan[adof[1], bdof[1]], nodes_sol.Ktan[adof[2], bdof[1]], nodes_sol.Ktan[adof[3], bdof[1]],
            nodes_sol.Ktan[adof[1], bdof[2]], nodes_sol.Ktan[adof[2], bdof[2]], nodes_sol.Ktan[adof[3], bdof[2]],
            nodes_sol.Ktan[adof[1], bdof[3]], nodes_sol.Ktan[adof[2], bdof[3]], nodes_sol.Ktan[adof[3], bdof[3]])
                
            Kba_loc = Kab_loc'

            Kab = Ra*Kab_loc
            Kba = Kba_loc*RaT

            @inbounds for (i, x) in enumerate(adof)
                @inbounds for (j, y) in enumerate(bdof)
        
                    nodes_sol.Ktan[x, y] = Kab[i, j]
                    nodes_sol.Ktan[y, x] = Kba[i, j]
        
                end 
            end

        end   
        
    end
    
end

# Convert the tangent matrix and the residual from cylindrical to carthesian coordinates
function dirichlet_local_to_global!(nodes_sol, allnodes)
    
    @inbounds for a = 1:length(allnodes)
        
        Ra = allnodes.R_global_to_local[a]       
        RaT = Ra'
        
        adof = 6*(a-1) .+ Vec3(1,2,3) 
        nodes_sol.ΔD[adof] .= RaT*nodes_sol.ΔD[adof]
        
    end
    
end

# Replaces uimp with the boundary conditions @t
function update_current_boundary_conditions!(uimp, t, conf)
    
    @inbounds for i in conf.bc.dof_disps
        uimp[i] = conf.bc.Fdisp(t)     
    end  

end 

function update_current_boundary_conditions_vector!(uimp, t, conf)
    
    Tstar =  convert(Int, round(t, RoundDown))
    Tstar1 = Tstar+1
    unew3 = read_ICs("outputGeometricalMorphing/u$Tstar1.txt")

    nnodesStent = length(unew3)
    unew = zeros(nnodesStent*3)
    for (j,i) in enumerate(1:3:nnodesStent*3) 
        unew[i+0] = unew3[j][1]
        unew[i+1] = unew3[j][2] 
        unew[i+2] = unew3[j][3]
    end 

    @inbounds for i in 1:length(uimp)
        uimp[i] = unew[i] 
    end  
    
end 

# Imposes single freedom constrains at the current step
function impose_BC_displacements!(nodes_sol, uimp, fixed_dofs)
    
    @inbounds for i in 1:length(fixed_dofs)
        
        idof = fixed_dofs[i]   
        D_prev = nodes_sol.asol[idof]   
        D_imposed = uimp[i]        
        ΔD_imposed = D_imposed - D_prev            
        nodes_sol.ΔD[idof] = ΔD_imposed

        @inbounds for i in 1:length(nodes_sol.r)  
            nodes_sol.r[i] = nodes_sol.r[i] .- ΔD_imposed .* nodes_sol.Ktan[i, idof] 
        end 
        
    end
    
end

#------------------------------------------------------
# FUNCTIONS USE ALONG WITH ROTATION MATRIX IN CRIMPING
#-------------------------------------------------------

# Update R_global_to_local for each node
function update_local_to_global_matrix!(allnodes)
    
    @inbounds for i in 1:length(allnodes)
        
        x = allnodes.pos[i]+ allnodes.u[i]
        theta = atan(x[2], x[1])
        allnodes.R_global_to_local[i] = Mat33(cos(theta), -sin(theta), 0,  sin(theta), cos(theta), 0, 0, 0, 1)
        
    end 
    
end

#------------------------------------------------
# FUNCTIONS USED IN THE CORRECTOR AND PREDICTOR
#------------------------------------------------

# At the beginning of the corrector, updates the NodalSolution global vectors with the current Configuration local vectors
function update_global_corrector!(nodes_sol, allnodes, disp_dof)
    
    @inbounds for n in allnodes

        nodes_sol.D[n.idof_disp] .= n.u
        nodes_sol.Ddt[n.idof_disp] .= n.udt
        nodes_sol.Ddtdt[n.idof_disp] .= n.udtdt
        nodes_sol.D[n.idof_ang] .= n.w
        nodes_sol.Ddt[n.idof_ang] .= n.wdt
        nodes_sol.Ddtdt[n.idof_ang] .= n.wdtdt

    end
    
end 

# At the beginning of the predictor, updates the NodalSolution global vectors with the current Configuration local vectors(not updating angles)
function update_global_predictor!(nodes_sol, allnodes)
        
    @inbounds for n in allnodes

        nodes_sol.D[n.idof_disp] .= n.u
        nodes_sol.Ddt[n.idof_disp] .= n.udt
        nodes_sol.Ddtdt[n.idof_disp] .= n.udtdt
        nodes_sol.D[n.idof_ang] .= 0
        nodes_sol.Ddt[n.idof_ang] .= n.wdt
        nodes_sol.Ddtdt[n.idof_ang] .= n.wdtdt

    end
    
    
    @inbounds for i in 1:length(nodes_sol.asol)

        nodes_sol.asol[i] = nodes_sol.D[i] 
        nodes_sol.ΔD[i] = 0  

    end 
    
end 

# At the end of the corrector, updates the Configuration local vectors with the NodalSolution global vectors computed during the current iteration
function update_local_corrector!(allnodes, ΔD_k, dt, nodes_sol, comp)
    
    beta = comp.beta
    gamma = comp.gamma
    
    @inbounds for i in 1:length(allnodes)

        allnodes.u[i] = nodes_sol.D[allnodes.idof_disp[i]]
        allnodes.udt[i] = nodes_sol.Ddt[allnodes.idof_disp[i]]
        allnodes.udtdt[i] = nodes_sol.Ddtdt[allnodes.idof_disp[i]]

        Sw = get_skew_skymmetric_matrix_from_vector(ΔD_k[allnodes.idof_ang[i]])
        allnodes.Delt[i] = exp(Sw)*allnodes.Delt[i]

        wdt_n = allnodes.wdt_n[i]
        wdtdt_n = allnodes.wdtdt_n[i]
        w_n1 = get_angle_from_rotation_matrix(allnodes.Delt[i])
        allnodes.wdt[i] = allnodes.Delt[i] * (gamma/(beta*dt)*w_n1 + (beta-gamma)/beta*wdt_n + dt*(beta-gamma/2)/beta*wdtdt_n)     
        allnodes.wdtdt[i] = allnodes.Delt[i] * (1/(beta*dt^2)*w_n1 - 1/(beta*dt)*wdt_n - (0.5-beta)/beta*wdtdt_n)

        allnodes.R[i] = allnodes.Delt[i]*allnodes.R_n[i] 

    end
    
end 

# At the end of the predictor, updates the Configuration local vectors with the NodalSolution global vectors predicted for the current time step
function update_local_predictor!(allnodes, nodes_sol)
    
    @inbounds for i in 1:length(allnodes)

        allnodes.u[i] = nodes_sol.D[allnodes.idof_disp[i]]
        allnodes.udt[i] = nodes_sol.Ddt[allnodes.idof_disp[i]]
        allnodes.udtdt[i] = nodes_sol.Ddtdt[allnodes.idof_disp[i]]

        allnodes.w[i] = nodes_sol.D[allnodes.idof_ang[i]]
        allnodes.wdt[i] = nodes_sol.Ddt[allnodes.idof_ang[i]]
        allnodes.wdtdt[i] = nodes_sol.Ddtdt[allnodes.idof_ang[i]]

        Sw = get_skew_skymmetric_matrix_from_vector(nodes_sol.D[allnodes.idof_ang[i]])
        allnodes.Delt[i] = exp(Sw)*allnodes.Delt[i]
        allnodes.R[i] = allnodes.Delt[i]*allnodes.R_n[i] 
        
    end
    
end 

# Update the displacement vectors with the solution of the linear system in the corrector
function update_nodal_solutions_corrector!(nodes_sol, disp_dof, gamma, beta, dt)
    
    @inbounds for i in disp_dof
        nodes_sol.D[i]  =  nodes_sol.D[i]  + nodes_sol.ΔD[i]
        nodes_sol.Ddt[i] =  nodes_sol.Ddt[i] + (gamma/(beta*dt))*nodes_sol.ΔD[i]
        nodes_sol.Ddtdt[i] = nodes_sol.Ddtdt[i] + (1/(beta*dt^2))*nodes_sol.ΔD[i]
    end 
    
end 

# Update the displacement vectors with the solution of the linear system in the predictor
function update_nodal_solutions_predictor!(nodes_sol, beta, gamma, dt)

    @inbounds for i in 1:length(nodes_sol.D)

        nodes_sol.Ddt_n[i] = nodes_sol.Ddt[i] # need to save the last velocity vector 
        nodes_sol.D[i] = nodes_sol.D[i] + nodes_sol.ΔD[i]
        nodes_sol.Ddt[i] = nodes_sol.Ddt[i] + (gamma/(beta*dt))*nodes_sol.ΔD[i] - (gamma/beta).*nodes_sol.Ddt[i] + (dt*(2*beta-gamma)/(2*beta))*nodes_sol.Ddtdt[i]
        nodes_sol.Ddtdt[i] = nodes_sol.Ddtdt[i] + (1/(beta*dt^2))*(nodes_sol.ΔD[i] - dt*nodes_sol.Ddt_n[i] - (dt^2)/2*nodes_sol.Ddtdt[i])

    end 
        
end 

# Update current solutions in the correct loop
function update_current_solution_corrector!(sol_n1, ndofs, matrices)
    
    @inbounds for i in 1:ndofs

        sol_n1.Tint[i] = matrices.Tint[i]
        sol_n1.Tk[i] =  matrices.Tk[i]
        sol_n1.Tct[i] =  matrices.Tct[i]   
        sol_n1.Tconstr[i] =  matrices.Tconstr[i]   

    end 
    
end 

# Compute residual and increment vector residual in the corrector loop
function compute_norms_corrector(k, aux_tol, sol_n1, nodes_sol, matrices, SHOW_COMP_TIME::Bool = false)

    aux_tol_old = aux_tol
    nodes_sol.f_aux .= sol_n1.fext .+ sol_n1.Tct .+ matrices.Tconstr .- matrices.Tdamp
    norm_f = norm(nodes_sol.f_aux)
    norm_res = norm(nodes_sol.r_free)
    
    if norm_f > 1e-1
        aux_tol = norm_res/norm_f
    else
        aux_tol = norm_res 
    end
    
    if (k == 1)
        norm_ddk = norm(nodes_sol.ΔD_free)
        
        if SHOW_COMP_TIME
            println("iteration $k, ||res|| = $aux_tol   ||ΔD_k[free_dof]|| = $norm_ddk")
        end
        
    else
        norm_ddk = norm(nodes_sol.ΔD_free)
        
        if SHOW_COMP_TIME
            frac_norm = log10(aux_tol_old/aux_tol)
            println("iteration $k, ||res|| = $aux_tol ($frac_norm)  ||ΔD_k[free_dof]|| = $norm_ddk")
        end
        
    end
    
    return aux_tol
    
end

# Update the tangent matrix in the predictor and corrector
function compute_Ktan_sparse!(nodes_sol, matrices, alpha, beta, gamma, dt)
    
    @inbounds for i in 1:length(nodes_sol.Ktan.nzval)
        nodes_sol.Ktan.nzval[i] = (1+alpha) * (matrices.Kint.nzval[i] - matrices.Kct.nzval[i] - matrices.Kconstr.nzval[i]) + (1/(beta*dt^2)) *  matrices.M.nzval[i] + (gamma/(beta*dt)) * (matrices.Ck.nzval[i] - matrices.Cconstr.nzval[i])
    end 

end 

# Update residual in the corrector loop
function compute_res_corrector!(nodes_sol, matrices, sol_n1, sol_n, alpha)
    
    @inbounds for i in 1:length(nodes_sol.r)
        nodes_sol.r[i] = (1+alpha) * (sol_n1.fext[i] + matrices.Tconstr[i] + matrices.Tct[i] - matrices.Tint[i]) - alpha * (sol_n.fext[i]  + sol_n.Tconstr[i]  + sol_n.Tct[i] - sol_n.Tint[i]) - matrices.Tk[i]
    end 
    
end 

# Update nodal solutions in the corrector loop
function update_nodal_solution_corrector_loop!(nodes_sol, disp_dof)
    
    nodes_sol.ΔD .= 0    
    nodes_sol.asol .= 0   
    
    @inbounds for i in disp_dof
        nodes_sol.asol[i] = nodes_sol.D[i] 
    end 
    
end

# Fill the free dofs residual vector (preallocated) to be used to solve the linear system
function fill_r_free!(nodes_sol, free_dof)
    
    @inbounds for (index, value) in enumerate(free_dof)
        nodes_sol.r_free[index] = nodes_sol.r[value] 
    end 
    
end 

# Fill the free dofs tangent matrix (preallocated) to be used to solve the linear system
function fill_Ktan_free!(nodes_sol)
     
    @inbounds for (index, value) in enumerate(nodes_sol.sparsity_map_free)
        nodes_sol.Ktan_free.nzval[index] = nodes_sol.Ktan.nzval[value]
    end

end

# Fill the free dofs of the whole displacements vector (preallocated) with the solution of the linear sysyem
function fill_ΔD_free_dofs!(nodes_sol, free_dof)
    
    @inbounds for (index,value) in enumerate(free_dof)
        nodes_sol.ΔD[value] = nodes_sol.ΔD_free[index]
    end 
    
end 
