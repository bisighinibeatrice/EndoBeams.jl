# ---------------------------------------------------------------
# This function ensures displacement, velocity, and acceleration states are 
# correctly set up for the next time step.
# ---------------------------------------------------------------

function initialize_global_corrector!(conf::BeamsConfiguration, state::SimulationState)
    for n in LazyRows(conf.nodes)
        # Reset incremental displacement for translational DOFs
        state.solⁿ⁺¹.ΔD[n.global_dofs_disp] .= 0  
        
        # Set displacement, velocity, and acceleration from nodal values
        state.solⁿ⁺¹.D[n.global_dofs_disp] .= n.u
        state.solⁿ⁺¹.Ḋ[n.global_dofs_disp] .= n.u̇
        state.solⁿ⁺¹.D̈[n.global_dofs_disp] .= n.ü
        
        # Reset incremental displacement for rotational DOFs
        state.solⁿ⁺¹.ΔD[n.global_dofs_rot] .= 0  
        
        # Set rotational displacement, velocity, and acceleration
        state.solⁿ⁺¹.D[n.global_dofs_rot] .= n.w
        state.solⁿ⁺¹.Ḋ[n.global_dofs_rot] .= n.ẇ
        state.solⁿ⁺¹.D̈[n.global_dofs_rot] .= n.ẅ
    end
end

# ---------------------------------------------------------------
# This function sets initial guesses for displacements, velocities, and accelerations.
# ---------------------------------------------------------------

function initialize_global_predictor!(conf::BeamsConfiguration, state::SimulationState)
    for n in LazyRows(conf.nodes)
        # Reset incremental displacement for translational DOFs
        state.solⁿ⁺¹.ΔD[n.global_dofs_disp] .= 0  
        
        # Set initial displacement, velocity, and acceleration for translational DOFs
        state.solⁿ⁺¹.D[n.global_dofs_disp] .= n.u
        state.solⁿ⁺¹.Ḋ[n.global_dofs_disp] .= n.u̇
        state.solⁿ⁺¹.D̈[n.global_dofs_disp] .= n.ü
        
        # Reset incremental displacement for rotational DOFs
        state.solⁿ⁺¹.ΔD[n.global_dofs_rot] .= 0  
        
        # Predict rotational displacement as zero but keep velocity and acceleration
        state.solⁿ⁺¹.D[n.global_dofs_rot] .= 0
        state.solⁿ⁺¹.Ḋ[n.global_dofs_rot] .= n.ẇ
        state.solⁿ⁺¹.D̈[n.global_dofs_rot] .= n.ẅ
    end
end 

# ---------------------------------------------------------------
# Computes the tangent matrix and residual forces in the predictor step.
# Uses Newmark parameters (α, β, γ) to construct the system.
# ---------------------------------------------------------------
function compute_tangent_and_residuals_predictor!(state::SimulationState, Δt, α, β, γ)

    # Compute tangent stiffness matrix (Ktan) at predictor step
    # Ktan = (1+α) * K + (γ/(βΔt)) * C + (1/(βΔt²)) * M
    @. state.solⁿ⁺¹.Ktan.nzval = (1 + α) * state.matricesⁿ⁺¹.K.nzval + (γ / (β * Δt)) * state.matricesⁿ⁺¹.C.nzval + (1 / (β * Δt^2)) * state.matricesⁿ⁺¹.M.nzval

    # Compute residual vector (r) at predictor step
    # r = (1+α) * external forces - (internal forces + contact forces) - α * previous external forces
    @. state.solⁿ⁺¹.r = (1 + α) * state.forcesⁿ⁺¹.fᵉˣᵗ - state.forcesⁿ.Tⁱⁿᵗ - state.forcesⁿ.Tᵏ - state.forcesⁿ.Tᶜ - state.forcesⁿ.Tᶜᵒⁿ - α * state.forcesⁿ.fᵉˣᵗ

    # Add contribution from damping matrix (C)
    # temp = (γ/β) * Ḋ - (Δt/2 * (2β - γ) / β) * D̈
    @. state.solⁿ⁺¹.temp = (γ / β) * state.solⁿ⁺¹.Ḋ - (Δt / 2 * (2β - γ) / β) * state.solⁿ⁺¹.D̈
    mul!(state.solⁿ⁺¹.r, state.matricesⁿ⁺¹.C, state.solⁿ⁺¹.temp, 1, 1)  # r += C * temp

    # Add contribution from mass matrix (M)
    # temp = Δt * Ḋ + (Δt² / 2) * D̈
    @. state.solⁿ⁺¹.temp = Δt * state.solⁿ⁺¹.Ḋ + (Δt^2 / 2) * state.solⁿ⁺¹.D̈
    mul!(state.solⁿ⁺¹.r, state.matricesⁿ⁺¹.M, state.solⁿ⁺¹.temp, 1 / (β * Δt^2), 1)  # r += (1 / (βΔt²)) * M * temp

end

# ---------------------------------------------------------------
# Computes the tangent matrix and residual forces in the corrector step.
# Used to refine the solution after the predictor step.
# ---------------------------------------------------------------
function compute_tangent_and_residuals_corrector!(state::SimulationState, Δt, α, β, γ)

    # Compute tangent stiffness matrix (Ktan) at corrector step (same as in predictor)
    @. state.solⁿ⁺¹.Ktan.nzval = (1 + α) * state.matricesⁿ⁺¹.K.nzval + (γ / (β * Δt)) * state.matricesⁿ⁺¹.C.nzval + (1 / (β * Δt^2)) * state.matricesⁿ⁺¹.M.nzval

    # Compute residual vector (r) at corrector step
    # r = (1+α) * (external + contact - internal forces) - α * (previous external + contact - internal forces) - previous contact force
    @. state.solⁿ⁺¹.r = (1 + α) * (state.forcesⁿ⁺¹.fᵉˣᵗ + state.forcesⁿ⁺¹.Tᶜ + state.forcesⁿ⁺¹.Tᶜᵒⁿ - state.forcesⁿ⁺¹.Tⁱⁿᵗ) - α * (state.forcesⁿ.fᵉˣᵗ + state.forcesⁿ.Tᶜ + state.forcesⁿ.Tᶜᵒⁿ - state.forcesⁿ.Tⁱⁿᵗ) - state.forcesⁿ⁺¹.Tᵏ

end

# ---------------------------------------------------------------
# Extracts the residual and tangent stiffness matrix for the free DOFs.
# ---------------------------------------------------------------
function extract_free_dofs!(conf::SimulationConfiguration, state::SimulationState)
    
    # Extract the residual vector (r) for free DOFs
    @views state.solⁿ⁺¹.r_free .= state.solⁿ⁺¹.r[conf.bcs.free_dofs]

    # Extract the corresponding entries of the tangent stiffness matrix (Ktan)
    @views state.solⁿ⁺¹.Ktan_free.nzval .= state.solⁿ⁺¹.Ktan.nzval[state.matricesⁿ⁺¹.sparsity_free]
    
end

# ---------------------------------------------------------------
# Solves for the displacement increments at the free DOFs.
# ---------------------------------------------------------------
function solve_free_dofs!(conf::SimulationConfiguration, state::SimulationState, solver)
    
    # 🔹 Perform symbolic analysis on the reduced stiffness matrix
    # This step prepares the solver for the system K ΔD = r.
    analyze!(solver, state.solⁿ⁺¹.Ktan_free, state.solⁿ⁺¹.r_free)
    
    # Solve for displacement increments at free DOFs
    # The solver finds ΔD_free by solving Ktan_free * ΔD_free = r_free
    solve!(solver, state.solⁿ⁺¹.ΔD_free, state.solⁿ⁺¹.Ktan_free, state.solⁿ⁺¹.r_free)

    # Store the computed displacement increments in the full displacement vector
    # Only the free DOFs are updated, keeping fixed DOFs unchanged.
    state.solⁿ⁺¹.ΔD[conf.bcs.free_dofs] .= state.solⁿ⁺¹.ΔD_free
end

# ---------------------------------------------------------------
# Computes the norms of the residual and displacement increment
# during the corrector loop for convergence checking.
# ---------------------------------------------------------------
function compute_norms_corrector(conf::SimulationConfiguration, state::SimulationState)
    
    # Extract the sets of free and fixed DOFs
    free_dofs = conf.bcs.free_dofs
    fixed_dofs = conf.bcs.fixed_dofs
    
    # Compute the norm of the displacement increment (ΔD)
    ΔD_norm = norm(state.solⁿ⁺¹.ΔD_free)
    
    # Compute the residual norm for free DOFs
    res_norm = norm(state.solⁿ⁺¹.r_free)
    
    # Compute the total external force norm for free DOFs
    f_norm = norm(state.forcesⁿ⁺¹.fᵉˣᵗ[free_dofs] .+ state.forcesⁿ⁺¹.Tᶜ[free_dofs] .+ state.forcesⁿ⁺¹.Tᵏ[free_dofs] .+ state.forcesⁿ⁺¹.Tᶜᵒⁿ[free_dofs])
    
    # Compute the total reaction force norm for fixed DOFs
    e_norm = norm(state.forcesⁿ⁺¹.fᵉˣᵗ[fixed_dofs] .+ state.forcesⁿ⁺¹.Tᶜ[fixed_dofs] .+ state.forcesⁿ⁺¹.Tᵏ[fixed_dofs] .+ state.forcesⁿ⁺¹.Tᶜᵒⁿ[fixed_dofs])
    
    # Normalize the residual norm for better convergence checks
    if f_norm + e_norm > 1e-12    
        res_norm = res_norm / (f_norm + e_norm)
    end
    
    # Return the computed norms for convergence assessment
    return res_norm, ΔD_norm   
end

# -------------------------------------------------------------------
# Updates the global displacement, velocity, and acceleration 
# at the end of the predictor step based on the computed increments.
# -------------------------------------------------------------------
function compute_global_predictor!(state::SimulationState, Δt, β, γ)
    
    # Precompute frequently used coefficients
    γ_over_β = γ / β
    γ_over_βΔt = γ / (β * Δt)
    inv_βΔt² = 1 / (β * Δt^2)
    Δt_over_2β = Δt * (2 * β - γ) / (2 * β)
    
    # Update global predictor values using Newmark-beta formulas
    @inbounds for i in eachindex(state.solⁿ⁺¹.ΔD)
        Ḋⁿ = deepcopy(state.solⁿ⁺¹.Ḋ[i])  # Store previous velocity to avoid overwriting during updates
        state.solⁿ⁺¹.D[i] += state.solⁿ⁺¹.ΔD[i] 
        state.solⁿ⁺¹.Ḋ[i] += γ_over_βΔt * state.solⁿ⁺¹.ΔD[i] - γ_over_β * Ḋⁿ + Δt_over_2β * state.solⁿ⁺¹.D̈[i]  
        state.solⁿ⁺¹.D̈[i] += inv_βΔt² * (state.solⁿ⁺¹.ΔD[i] - Δt * Ḋⁿ - (Δt^2 / 2) * state.solⁿ⁺¹.D̈[i])  
    end
end

# -------------------------------------------------------------------
# Updates the global displacement, velocity, and acceleration 
# at the end of the corrector step based on the computed increments.
# -------------------------------------------------------------------
function compute_global_corrector!(state::SimulationState, Δt, β, γ)
    
    # Precompute frequently used coefficients 
    inv_βΔt² = 1 / (β * Δt^2)
    γ_over_βΔt = γ / (β * Δt)

    # 🔹 Update global corrector values
    @inbounds for i in eachindex(state.solⁿ⁺¹.ΔD)
        state.solⁿ⁺¹.D[i] += state.solⁿ⁺¹.ΔD[i]
        state.solⁿ⁺¹.Ḋ[i] += γ_over_βΔt * state.solⁿ⁺¹.ΔD[i]  
        state.solⁿ⁺¹.D̈[i] += inv_βΔt² * state.solⁿ⁺¹.ΔD[i]  
    end
end 

# -------------------------------------------------------------------
# Updates local vectors for different configurations 
# based on the global predictor solution.
# -------------------------------------------------------------------

function compute_local_predictor!(conf::BeamsConfiguration, state::SimulationState, γ=nothing, β=nothing, Δt=nothing)
    
    @inbounds for i in eachindex(conf.nodes)
        # Update displacement, velocity, and acceleration
        conf.nodes.u[i] = state.solⁿ⁺¹.D[conf.nodes.global_dofs_disp[i]]
        conf.nodes.u̇[i] = state.solⁿ⁺¹.Ḋ[conf.nodes.global_dofs_disp[i]]
        conf.nodes.ü[i] = state.solⁿ⁺¹.D̈[conf.nodes.global_dofs_disp[i]]
        
        # Update rotational degrees of freedom
        θ̃ = state.solⁿ⁺¹.D[conf.nodes.global_dofs_rot[i]]
        conf.nodes.w[i] = θ̃
        conf.nodes.ΔR[i] = rotation_matrix(θ̃)  # Compute incremental rotation matrix
        conf.nodes.R[i] = conf.nodes.ΔR[i] * conf.nodes.Rⁿ[i]  # Update total rotation

        # Compute angular velocity and acceleration
        ẇⁿ = conf.nodes.ẇⁿ[i]
        ẅⁿ = conf.nodes.ẅⁿ[i]
        conf.nodes.ẇ[i] = conf.nodes.ΔR[i] * (γ/(β*Δt)*θ̃ + (β-γ)/β*ẇⁿ + Δt*(β-γ/2)/β*ẅⁿ)     
        conf.nodes.ẅ[i] = conf.nodes.ΔR[i] * (1/(β*Δt^2)*θ̃ - 1/(β*Δt)*ẇⁿ - (1/(2*β)-1)*ẅⁿ) 
    end  
end

# -------------------------------------------------------------------
# Updates local vectors for different configurations 
# based on the computed corrector solution.
# -------------------------------------------------------------------

function compute_local_corrector!(conf::BeamsConfiguration, state::SimulationState, γ=nothing, β=nothing, Δt=nothing)
    
    @inbounds for i in eachindex(conf.nodes)
        # Update displacement, velocity, and acceleration
        conf.nodes.u[i] = state.solⁿ⁺¹.D[conf.nodes.global_dofs_disp[i]]
        conf.nodes.u̇[i] = state.solⁿ⁺¹.Ḋ[conf.nodes.global_dofs_disp[i]]
        conf.nodes.ü[i] = state.solⁿ⁺¹.D̈[conf.nodes.global_dofs_disp[i]]
        
        # Compute the incremental rotation and apply it
        conf.nodes.ΔR[i] = rotation_matrix(state.solⁿ⁺¹.ΔD[conf.nodes.global_dofs_rot[i]]) * conf.nodes.ΔR[i]

        # Compute new angular velocity and acceleration
        ẇⁿ = conf.nodes.ẇⁿ[i]
        ẅⁿ = conf.nodes.ẅⁿ[i]
        wⁿ⁺¹ = toangle(conf.nodes.ΔR[i])
        conf.nodes.ẇ[i] = conf.nodes.ΔR[i] * (γ/(β*Δt)*wⁿ⁺¹ + (β-γ)/β*ẇⁿ + Δt*(β-γ/2)/β*ẅⁿ)     
        conf.nodes.ẅ[i] = conf.nodes.ΔR[i] * (1/(β*Δt^2)*wⁿ⁺¹ - 1/(β*Δt)*ẇⁿ - (1/(2*β)-1)*ẅⁿ)
        
        # Update total rotation
        conf.nodes.R[i] = conf.nodes.ΔR[i] * conf.nodes.Rⁿ[i] 
    end  
end 

# -------------------------------------------------------------------
# Updates nodal values and forces for the next time step 
# upon convergence for different configurations.
# -------------------------------------------------------------------

function update_converged!(conf::BeamsConfiguration, state::SimulationState)
    
    @inbounds for i in eachindex(conf.nodes)
        # Update nodal displacements, velocities, and accelerations to the converged values
        conf.nodes.uⁿ[i] = conf.nodes.u[i]
        conf.nodes.u̇ⁿ[i] = conf.nodes.u̇[i]
        conf.nodes.üⁿ[i] = conf.nodes.ü[i]
        conf.nodes.wⁿ[i] = conf.nodes.w[i]
        conf.nodes.ẇⁿ[i] = conf.nodes.ẇ[i]
        conf.nodes.ẅⁿ[i] = conf.nodes.ẅ[i]
        conf.nodes.Rⁿ[i] = conf.nodes.R[i]
        conf.nodes.ΔRⁿ[i] = conf.nodes.ΔR[i]
    end
    
    # Update the force vectors for the current time step
    state.forcesⁿ.Tⁱⁿᵗ .= state.forcesⁿ⁺¹.Tⁱⁿᵗ
    state.forcesⁿ.Tᵏ .= state.forcesⁿ⁺¹.Tᵏ
    state.forcesⁿ.Tᶜ .= state.forcesⁿ⁺¹.Tᶜ
    state.forcesⁿ.fᵉˣᵗ .= state.forcesⁿ⁺¹.fᵉˣᵗ
    state.forcesⁿ.Tᶜᵒⁿ .= state.forcesⁿ⁺¹.Tᶜᵒⁿ
    
end

# -------------------------------------------------------------------
# Reverts nodal values and forces for the next time step 
# if convergence is not achieved.
# -------------------------------------------------------------------

function update_not_converged!(conf::BeamsConfiguration, state::SimulationState)
    
    @inbounds for i in eachindex(conf.nodes)
        # Revert nodal displacements, velocities, and accelerations
        conf.nodes.u[i] = conf.nodes.uⁿ[i]
        conf.nodes.u̇[i] = conf.nodes.u̇ⁿ[i]
        conf.nodes.ü[i] = conf.nodes.üⁿ[i]
        conf.nodes.w[i] = conf.nodes.wⁿ[i]
        conf.nodes.ẇ[i] = conf.nodes.ẇⁿ[i]
        conf.nodes.ẅ[i] = conf.nodes.ẅⁿ[i]
        conf.nodes.R[i] = conf.nodes.Rⁿ[i]
        conf.nodes.ΔR[i] = conf.nodes.ΔRⁿ[i]
    end
    
    # Revert force vectors
    state.forcesⁿ⁺¹.Tⁱⁿᵗ .= state.forcesⁿ.Tⁱⁿᵗ
    state.forcesⁿ⁺¹.Tᵏ .= state.forcesⁿ.Tᵏ
    state.forcesⁿ⁺¹.Tᶜ .= state.forcesⁿ.Tᶜ
    state.forcesⁿ⁺¹.fᵉˣᵗ .= state.forcesⁿ.fᵉˣᵗ
    state.forcesⁿ⁺¹.Tᶜᵒⁿ .= state.forcesⁿ.Tᶜᵒⁿ
    
end