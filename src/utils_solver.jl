#----------------------------------
# FUNCTIONS USED IN THE TIME LOOP 
#----------------------------------

# Cleans the output folders from files of the precedent computation
function clean_folders(output_dir)

    if output_dir != ""
        dir = pwd()
        cd(output_dir)
        foreach(rm, filter(endswith(".vtu"), readdir()))
        foreach(rm, filter(endswith(".pvd"), readdir()))
        cd(dir)
    end 
    
end 

# Cleans folders, pre-allocate and initialise the variables used during the simulation and save the VTKs of the initial configuration
function solver_initialisation(conf, nodes, beams, constraints, sdf, output_dir, T=Float64)

    solⁿ = Solution(conf, T)
    solⁿ⁺¹ = deepcopy(solⁿ)
    energy = Energy(T) 
    matrices, nodes_sol = constructor_sparse_matrices!(beams, nodes, constraints, conf)

    vtkdata = VTKData(length(beams), output_dir, sdf, T)

    return solⁿ, solⁿ⁺¹, nodes_sol, matrices, energy, vtkdata

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


# Imposes single freedom constraints at the current step
function apply_BCs!(nodes_sol, bcs)
    
    @inbounds for idof in bcs.fixed_dofs
         
        D_prev = nodes_sol.D[idof]   
        D_imposed = bcs.disp_vals[idof]
        ΔD_imposed = D_imposed - D_prev            
        nodes_sol.ΔD[idof] = ΔD_imposed

        for j in nzrange(nodes_sol.Ktan_mat, idof)
            nodes_sol.r[nodes_sol.Ktan_mat.rowval[j]] -= ΔD_imposed * nodes_sol.Ktan[j]
        end

        
    end

    
end


#------------------------------------------------
# FUNCTIONS USED IN THE CORRECTOR AND PREDICTOR
#------------------------------------------------



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



function residuals_corrector!(nodes_sol, matrices, solⁿ, solⁿ⁺¹, Δt, α, β, γ)
        
    @. nodes_sol.r = (1+α) * (solⁿ⁺¹.fᵉˣᵗ + matrices.Tᶜᵒⁿ + matrices.Tᶜ - matrices.Tⁱⁿᵗ) - α * (solⁿ.fᵉˣᵗ + solⁿ.Tᶜᵒⁿ + solⁿ.Tᶜ - solⁿ.Tⁱⁿᵗ) - matrices.Tᵏ

end
