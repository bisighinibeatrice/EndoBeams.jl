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

#----------------------------------------------------
# FUNCTIONS FOR CHANGING COORDINATE SYSTEM 
#----------------------------------------------------

# Transform Cartesian coordinates into cylindrical coordinates 
function apply_cylindrical_coordinate_system!(conf, state)
    
    @inbounds for node_a in conf.nodes
        Ra = node_a.R_carthesian_to_cylindrical      # Rotation matrix from global to local
        RaT = Ra'                          # Transpose for reverse transformation
        dofs_a = node_a.global_dofs_disp[1:3]          # Displacement DOFs for node_a

        # Transform displacement vectors to local cylindrical coordinates
        state.solⁿ⁺¹.r[dofs_a] .= Ra * state.solⁿ⁺¹.r[dofs_a]
        state.solⁿ⁺¹.D[dofs_a] .= Ra * state.solⁿ⁺¹.D[dofs_a]

        # Update stiffness matrix in cylindrical coordinates
        @inbounds for node_b in conf.nodes
            dofs_b = node_b.global_dofs_disp

            # Extract local submatrix from global stiffness matrix
            Kab_local = Mat33( 
                state.solⁿ⁺¹. Ktan[dofs_a[1], dofs_b[1]], state.solⁿ⁺¹. Ktan[dofs_a[2], dofs_b[1]], state.solⁿ⁺¹. Ktan[dofs_a[3], dofs_b[1]],
                state.solⁿ⁺¹. Ktan[dofs_a[1], dofs_b[2]], state.solⁿ⁺¹. Ktan[dofs_a[2], dofs_b[2]], state.solⁿ⁺¹. Ktan[dofs_a[3], dofs_b[2]],
                state.solⁿ⁺¹. Ktan[dofs_a[1], dofs_b[3]], state.solⁿ⁺¹. Ktan[dofs_a[2], dofs_b[3]], state.solⁿ⁺¹. Ktan[dofs_a[3], dofs_b[3]]
            )

            Kba_local = Kab_local'               # Transpose for symmetric update

            Kab = Ra * Kab_local                 # Rotate to cylindrical coords
            Kba = Kba_local * RaT

            # Store transformed submatrices back into global stiffness matrix
            @inbounds for (i, dof_i) in enumerate(dofs_a)
                @inbounds for (j, dof_j) in enumerate(dofs_b)
                    state.solⁿ⁺¹. Ktan[dof_i, dof_j] = Kab[i, j]
                    state.solⁿ⁺¹. Ktan[dof_j, dof_i] = Kba[i, j]  # Symmetric part
                end
            end
        end 
    end   

end

# Transform cylindrical coordinates back into Cartesian coordinates
function revert_to_carthesian_coordinate_system!(conf, state)

    @inbounds for node in conf.nodes
        RaT = node.R_carthesian_to_cylindrical'  # Transpose of rotation matrix

        dofs = node.global_dofs_disp

        # Rotate displacement increments and vectors back to global Cartesian coords
        state.solⁿ⁺¹.ΔD[dofs] = RaT * state.solⁿ⁺¹.ΔD[dofs]
        state.solⁿ⁺¹.D[dofs]  = RaT * state.solⁿ⁺¹.D[dofs]
    end

end