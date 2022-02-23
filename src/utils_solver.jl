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
function save_VTK(i, nodes, beams, sol_GP, int_pos, int_conn, dirOutput, SAVE_INTERPOLATION_VTK = false, SAVE_GP_VTK = false, SAVE_NODES_VTK = false)
    
    if SAVE_INTERPOLATION_VTK    
        write_VTK_beams(i, nodes, beams, int_pos, int_conn, dirOutput)
    end 
    
    if SAVE_GP_VTK 
        write_VTK_GP(i, sol_GP, dirOutput) 
    end 
    
    if SAVE_NODES_VTK 
        write_VTK_nodes(i, nodes, beams, dirOutput)
    end 

end 

# Cleans folders, pre-allocate and initialise the variables used during the simulation and save the VTKs of the initial configuration
function solver_initialisation(conf, nodes, beams, constraints, thisDirOutputPath, SAVE_INTERPOLATION_VTK = false, SAVE_GP_VTK = false, SAVE_NODES_VTK = false, T=Float64)

    clean_folders(thisDirOutputPath)

    solⁿ = constructor_solution(conf, T)
    solⁿ⁺¹ = deepcopy(solⁿ)
    # sol_GP = constructor_solution_GP(length(beams), T)
    sol_GP = 0
    energy = constructor_energy(T) 
    matrices, nodes_sol = constructor_sparse_matrices!(beams, nodes, constraints, conf)

    #save_VTK(0, nodes, beams, thisDirOutputPath, SAVE_INTERPOLATION_VTK,  SAVE_GP_VTK, SAVE_NODES_VTK)

    return solⁿ, solⁿ⁺¹, nodes_sol, matrices, energy

end 

# Save energy evolution in text file
function save_energy(plot_Phi_energy, plot_K_energy, plot_contact_energy, outputDir)

    open(outputDir * "/PHI ENERGY.txt", "w") do Iₒ
        writedlm(Iₒ, plot_Phi_energy)
    end

    open(outputDir * "/K ENERGY.txt", "w") do Iₒ
        writedlm(Iₒ, plot_K_energy)
    end

    open(outputDir * "/C ENERGY.txt", "w") do Iₒ
        writedlm(Iₒ, plot_contact_energy)
    end

end

# Update nodes values @n if converged
function  update_nodes_converged!(nodes, solⁿ, matrices)

    @inbounds for i in eachindex(nodes)      

        nodes.uⁿ[i] = nodes.u[i]
        nodes.u̇ⁿ[i] = nodes.u̇[i]
        nodes.üⁿ[i] = nodes.ü[i]
        nodes.wⁿ[i] = nodes.w[i]
        nodes.ẇⁿ[i] = nodes.ẇ[i]
        nodes.ẅⁿ[i] = nodes.ẅ[i]
        nodes.Rⁿ[i] = nodes.R[i]
        nodes.ΔRⁿ[i] = nodes.ΔR[i]

    end

    # update the current solution vectors (used in the next time step)
    solⁿ.Tⁱⁿᵗ .= matrices.Tⁱⁿᵗ
    solⁿ.Tᵏ .=  matrices.Tᵏ
    solⁿ.Tᶜ .=  matrices.Tᶜ
    solⁿ.Tᶜᵒⁿ .=  matrices.Tᶜᵒⁿ

end

# Update nodes values @n+1 if NOT converged
function update_nodes_not_converged!(nodes)

    @inbounds for i in eachindex(nodes)

        nodes.u[i] = nodes.uⁿ[i]
        nodes.u̇[i] = nodes.u̇ⁿ[i]
        nodes.ü[i] = nodes.üⁿ[i]
        nodes.w[i] = nodes.wⁿ[i]
        nodes.ẇ[i] = nodes.ẇⁿ[i]
        nodes.ẅ[i] = nodes.ẅⁿ[i]
        nodes.R[i] = nodes.Rⁿ[i]
        nodes.ΔR[i] = nodes.ΔRⁿ[i]

    end

end

#----------------------------------------
# FUNCTIONS TO UPDATE THE EXTERNAL LOADS
#----------------------------------------

# Replaces fᵉˣᵗ with the external force @t
function get_current_external_force!(fᵉˣᵗ, t, conf)

    for i in conf.ext_forces.loaded_dofs
        fᵉˣᵗ[i] = conf.ext_forces.F(t)     
    end

end 

# Replaces fᵉˣᵗ with the external force @t for the radial crimping
function get_current_external_force_crimping!(fᵉˣᵗ, t, conf, nodes)
    
    update_local_to_global_matrix!(nodes)
    
    for n in nodes

        dof_dispⁿ = n.idof_disp
        fᵉˣᵗⁿ =  [conf.ext_forces.F(t), 0, 0] 
        fᵉˣᵗⁿ = n.R_global_to_local' * fᵉˣᵗⁿ
        fᵉˣᵗ[dof_dispⁿ] .= fᵉˣᵗⁿ   

    end 
   
end 

# Update the solution external force @t
function update_current_solution_external_force!(solⁿ⁺¹, t, conf)

    get_current_external_force!(solⁿ⁺¹.fᵉˣᵗ, t, conf)

end 

# Update the solution external force @t for the radial crimping
function update_current_solution_external_force_crimping!(solⁿ⁺¹, t, conf, nodes)

    get_current_external_force_crimping!(solⁿ⁺¹.fᵉˣᵗ, t, conf, nodes)
    
end 

#----------------------------------------
# FUNCTIONS TO UPDATE AND IMPOSED THE BCs
#----------------------------------------

# Convert the tangent matrix and the residual from carthesian to cylindrical coordinates
function dirichlet_global_to_local!(nodes_sol, nodes)
    
    @inbounds for a = 1:length(nodes)
        
        Ra = nodes.R_global_to_local[a]    
        RaT = Ra'
        
        adof = 6*(a-1) .+ Vec3(1,2,3)

        nodes_sol.r[adof] .= Ra*nodes_sol.r[adof]
        nodes_sol.asol[adof] .= Ra*nodes_sol.asol[adof]
        
        @inbounds for b = 1:length(nodes)
            
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
function dirichlet_local_to_global!(nodes_sol, nodes)
    
    @inbounds for i in eachindex(nodes)
        
        R = nodes.R_global_to_local[i]       
        
        idof = 6*(i-1) .+ Vec3(1,2,3) 
        nodes_sol.ΔD[idof] .= R'*nodes_sol.ΔD[idof]
        
    end
    
end

# Replaces uimp with the boundary conditions @t
function update_current_boundary_conditions!(uimp, t, conf)
    
    @inbounds for i in conf.bc.dof_disps
        uimp[i] = conf.bc.Fdisp(t)     
    end  

end 

function update_current_boundary_conditions_vector!(uimp::AbstractVector{T}, t, conf) where T
    
    Tstar =  convert(Int, round(t, RoundDown))
    Tstar1 = Tstar+1
    unew3 = read_ICs("outputGeometricalMorphing/u$Tstar1.txt", T)

    nnodesStent = length(unew3)
    for (j,i) in enumerate(1:3:nnodesStent*3) 
        uimp[i+0] = unew3[j][1]
        uimp[i+1] = unew3[j][2] 
        uimp[i+2] = unew3[j][3]
    end 
    
end 

# Imposes single freedom constraints at the current step
function apply_BCs!(nodes_sol, bcs)
    
    @inbounds for i in eachindex(bcs.fixed_dofs)
        
        idof = fixed_dofs[i]   
        D_prev = nodes_sol.D[idof]   
        D_imposed = bcs.disp_vals[idof]
        ΔD_imposed = D_imposed - D_prev            
        nodes_sol.ΔD[idof] = ΔD_imposed

        for j in nzrange(nodes_sol.Ktan_mat, idof)
            nodes_sol.r[nodes_sol.Ktan_mat.rowval[j]] -= ΔD_imposed * nodes_sol.Ktan[j]
        end

        
    end

    
end

#------------------------------------------------------
# FUNCTIONS USE ALONG WITH ROTATION MATRIX IN CRIMPING
#-------------------------------------------------------

# Update R_global_to_local for each node
function update_local_to_global_matrix!(nodes)
    
    @inbounds for i in eachindex(nodes)
        
        x = nodes.X₀[i]+ nodes.u[i]
        theta = atan(x[2], x[1])
        nodes.R_global_to_local[i] = Mat33(cos(theta), -sin(theta), 0,  sin(theta), cos(theta), 0, 0, 0, 1)
        
    end 
    
end

#------------------------------------------------
# FUNCTIONS USED IN THE CORRECTOR AND PREDICTOR
#------------------------------------------------

# At the beginning of the corrector, updates the NodalSolution global vectors with the current Configuration local vectors
function update_global_corrector!(nodes_sol, nodes, disp_dof)
    
    @inbounds for n in nodes

        nodes_sol.D[n.idof_disp] .= n.u
        nodes_sol.Ḋ[n.idof_disp] .= n.u̇
        nodes_sol.D̈[n.idof_disp] .= n.ü
        nodes_sol.D[n.idof_rot] .= n.w
        nodes_sol.Ḋ[n.idof_rot] .= n.ẇ
        nodes_sol.D̈[n.idof_rot] .= n.ẅ

    end
    
end 



# At the end of the corrector, updates the Configuration local vectors with the NodalSolution global vectors computed during the current iteration
function update_local_corrector!(nodes, nodes_sol, Δt, β, γ)
    
    
    @inbounds for i in eachindex(nodes)

        nodes.u[i] = nodes_sol.D[nodes.idof_disp[i]]
        nodes.u̇[i] = nodes_sol.Ḋ[nodes.idof_disp[i]]
        nodes.ü[i] = nodes_sol.D̈[nodes.idof_disp[i]]

        nodes.ΔR[i] = rotation_matrix(nodes_sol.ΔD[nodes.idof_rot[i]]) * nodes.ΔR[i]

        ẇⁿ = nodes.ẇⁿ[i]
        ẅⁿ = nodes.ẅⁿ[i]
        wⁿ⁺¹ = toangle(nodes.ΔR[i])
        nodes.ẇ[i] = nodes.ΔR[i] * (γ/(β*Δt)*wⁿ⁺¹ + (β-γ)/β*ẇⁿ + Δt*(β-γ/2)/β*ẅⁿ)     
        nodes.ẅ[i] = nodes.ΔR[i] * (1/(β*Δt^2)*wⁿ⁺¹ - 1/(β*Δt)*ẇⁿ - (1/(2*β)-1)*ẅⁿ)

        nodes.R[i] = nodes.ΔR[i]*nodes.Rⁿ[i] 

    end
    
end 

# At the end of the predictor, updates the Configuration local vectors with the NodalSolution global vectors predicted for the current time step
function update_local_predictor!(nodes, nodes_sol, Δt, β, γ)

    
    @inbounds for i in eachindex(nodes)

        disp_dofs = nodes.idof_disp[i]
        dofs_rot = nodes.idof_rot[i]

        nodes.u[i] = nodes_sol.D[disp_dofs]
        nodes.u̇[i] = nodes_sol.Ḋ[disp_dofs]
        nodes.ü[i] = nodes_sol.D̈[disp_dofs]

        θ̃ = nodes_sol.D[dofs_rot]
        nodes.w[i] = θ̃
        nodes.ΔR[i] = rotation_matrix(θ̃)
        nodes.R[i] = nodes.ΔR[i]*nodes.Rⁿ[i]

        ẇⁿ = nodes.ẇⁿ[i]
        ẅⁿ = nodes.ẅⁿ[i]
        nodes.ẇ[i] = nodes.ΔR[i] * (γ/(β*Δt)*θ̃ + (β-γ)/β*ẇⁿ + Δt*(β-γ/2)/β*ẅⁿ)     
        nodes.ẅ[i] = nodes.ΔR[i] * (1/(β*Δt^2)*θ̃ - 1/(β*Δt)*ẇⁿ - (1/(2*β)-1)*ẅⁿ)
        
        
    end

    
end 

# Update the displacement vectors with the solution of the linear system in the corrector
function update_nodal_solutions_corrector!(nodes_sol, disp_dof, γ, β, Δt)
    
    @inbounds for i in disp_dof
        nodes_sol.D[i]  =  nodes_sol.D[i]  + nodes_sol.ΔD[i]
        nodes_sol.Ḋ[i] =  nodes_sol.Ḋ[i] + (γ/(β*Δt))*nodes_sol.ΔD[i]
        nodes_sol.D̈[i] = nodes_sol.D̈[i] + (1/(β*Δt^2))*nodes_sol.ΔD[i]
    end 
    
end 

# Update the displacement vectors with the solution of the linear system in the predictor
function update_nodal_solutions_predictor!(nodes_sol, β, γ, Δt)

    @inbounds for i in eachindex(nodes_sol.D)

        nodes_sol.Ḋⁿ[i] = nodes_sol.Ḋ[i] # need to save the last velocity vector 
        nodes_sol.D[i] = nodes_sol.D[i] + nodes_sol.ΔD[i]
        nodes_sol.Ḋ[i] = nodes_sol.Ḋ[i] + (γ/(β*Δt))*nodes_sol.ΔD[i] - (γ/β)*nodes_sol.Ḋ[i] + (Δt*(2*β-γ)/(2*β))*nodes_sol.D̈[i]
        nodes_sol.D̈[i] = nodes_sol.D̈[i] + (1/(β*Δt^2))*(nodes_sol.ΔD[i] - Δt*nodes_sol.Ḋⁿ[i] - (Δt^2)/2*nodes_sol.D̈[i])

    end 
        
end 

# Update current solutions in the correct loop
function update_current_solution_corrector!(solⁿ⁺¹, ndofs, matrices)
    
    @inbounds for i in 1:ndofs

        solⁿ⁺¹.Tⁱⁿᵗ[i] = matrices.Tⁱⁿᵗ[i]
        solⁿ⁺¹.Tᵏ[i] =  matrices.Tᵏ[i]
        solⁿ⁺¹.Tᶜ[i] =  matrices.Tᶜ[i]   
        solⁿ⁺¹.Tᶜᵒⁿ[i] =  matrices.Tᶜᵒⁿ[i]   

    end 
    
end 

# Compute residual and increment vector residual in the corrector loop
function compute_norms_corrector(k, solⁿ⁺¹, nodes_sol, matrices, SHOW_COMP_TIME = false)
    
    f_norm = mapreduce((f...) -> sum(f^2 for f in f), +, solⁿ⁺¹.fᵉˣᵗ, solⁿ⁺¹.Tᶜ, matrices.Tᶜᵒⁿ)
    f_norm = sqrt(f_norm)
    res_norm = norm(nodes_sol.r_free)
    
    res_norm = f_norm > 1e-1 ? res_norm/f_norm : res_norm
    
    ΔD_norm = norm(nodes_sol.ΔD_free)
    if SHOW_COMP_TIME
        k == 1 && @printf "%4s\t%8s\t%8s\n" "iter" "‖res‖" "‖ΔD‖"
        @printf "%4d\t%1.2e\t%1.2e\n"  k  res_norm  ΔD_norm 
    end
    
    return res_norm, ΔD_norm
    
end



function tangent_and_residuals_predictor!(nodes_sol, matrices, solⁿ, solⁿ⁺¹, Δt, α, β, γ)
        
        @. nodes_sol.Ktan = (1+α) * matrices.K + (1/(β*Δt^2)) * matrices.M + (γ/(β*Δt)) * matrices.C

        @. nodes_sol.r =  (1+α) * solⁿ⁺¹.fᵉˣᵗ - matrices.Tⁱⁿᵗ - matrices.Tᶜ - matrices.Tᶜᵒⁿ - matrices.Tᵏ - α * solⁿ.fᵉˣᵗ
        @. nodes_sol.temp = γ/β * nodes_sol.Ḋ - (Δt/2*(2*β-γ)/β) * nodes_sol.D̈
        mul!(nodes_sol.r, matrices.C_mat, nodes_sol.temp, 1, 1)
        @. nodes_sol.temp =  Δt * nodes_sol.Ḋ + Δt^2/2 * nodes_sol.D̈
        mul!(nodes_sol.r, matrices.M_mat, nodes_sol.temp, 1/(β*Δt^2), 1)

end



function tangent_and_residuals_corrector!(nodes_sol, matrices, solⁿ, solⁿ⁺¹, Δt, α, β, γ)
        
    @. nodes_sol.Ktan = (1+α) * matrices.K + (1/(β*Δt^2)) * matrices.M + (γ/(β*Δt)) * matrices.C
    @. nodes_sol.r = (1+α) * (solⁿ⁺¹.fᵉˣᵗ + matrices.Tᶜᵒⁿ + matrices.Tᶜ - matrices.Tⁱⁿᵗ) - α * (solⁿ.fᵉˣᵗ + solⁿ.Tᶜᵒⁿ + solⁿ.Tᶜ - solⁿ.Tⁱⁿᵗ) - matrices.Tᵏ

end

