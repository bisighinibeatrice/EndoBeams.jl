#----------------------------------------------------
# FUNCTIONS FOR APPLYING BOUNDARY CONDITIONS 
#----------------------------------------------------

function apply_boundary_conditions!(bcs::BoundaryConditions, state::SimulationState)

    @inbounds for (index, dof) in enumerate(bcs.fixed_dofs) 
        prev_displacement = state.solⁿ⁺¹.D[dof]  # Previous displacement at DOF
        imposed_displacement = bcs.imposed_displacements[index] # Encastre (fixed) condition: displacement = 0

        imposed_displacement_change = imposed_displacement - prev_displacement  
        
        # Store the change in displacement
        state.solⁿ⁺¹.ΔD[dof] = imposed_displacement_change
        
        # Get the nonzero elements of the stiffness matrix for the current DOF
        nonzero_indices = nzrange(state.solⁿ⁺¹.Ktan, dof)
        
        # Loop through each nonzero entry in the stiffness matrix row corresponding to `dof`
        for j in nonzero_indices
            row_idx = state.solⁿ⁺¹.Ktan.rowval[j]  # Get the row index
            stiffness_value = state.solⁿ⁺¹.Ktan.nzval[j]     # Get the stiffness value
            
            # Modify the residual force vector
            state.solⁿ⁺¹.r[row_idx] -= imposed_displacement_change * stiffness_value 
        end
    end

end

#Applies fixed (encastre) boundary conditions at the current step
function apply_boundary_conditions!(conf::BeamsConfiguration, state::SimulationState)
    
    apply_boundary_conditions!(conf.bcs, state)

end

#----------------------------------------------------
# FUNCTIONS FOR UPDATING BOUNDARY CONDITIONS 
#----------------------------------------------------

function update_boundary_conditions!(conf::BeamsConfiguration, tⁿ⁺¹)
      
    # Update external forces at the current and intermediate time steps
    for i in conf.bcs.displaced_indices
        conf.bcs.imposed_displacements[i + conf.offset_dofs[]] = conf.bcs.imposed_displacements_function(tⁿ⁺¹)  # Apply force to the corresponding DOF
    end

end

#----------------------------------------------------
# FUNCTIONS FOR UPDATING LOADS
#----------------------------------------------------

# Updates external force values for a given mesh configuration
function update_loads!(conf::BeamsConfiguration, state::SimulationState, tⁿ⁺¹)
        
    # Update external forces at the current and intermediate time steps
    if conf.loads !== nothing && conf.loads.concentrated_force !== nothing
        for i in conf.loads.concentrated_force.loaded_dofs
            state.forcesⁿ⁺¹.fᵉˣᵗ[i + conf.offset_dofs[]] = conf.loads.concentrated_force.force_function(tⁿ⁺¹, i)  # Apply force to the corresponding DOF
        end
    end
    
end