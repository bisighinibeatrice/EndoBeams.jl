#----------------------------------------------------
# Applies multi-degree-of-freedom constraint forces
# at the current simulation step for one constraint
#----------------------------------------------------
function compute_constraint!(state, nodes, constraint)

    k = constraint.k      # Stiffness coefficient
    η = constraint.η      # Damping coefficient

    i1 = constraint.node1 # Index of first node
    i2 = constraint.node2 # Index of second node

    dofs_a = nodes.global_dofs[i1]  # DOFs of node 1
    dofs_b = nodes.global_dofs[i2]  # DOFs of node 2

    xa = nodes.u[i1]   # Displacement of node 1
    xb = nodes.u[i2]   # Displacement of node 2
    va = nodes.u̇[i1]   # Velocity of node 1
    vb = nodes.u̇[i2]   # Velocity of node 2

    # Compute constraint force vector (ta)
    # Combines elastic (k*(xb - xa)) and damping (ηk*(vb - va)) terms
    ta = k * (xb - xa) + η * k * (vb - va)

    # Apply force to node 1 (positive contribution)
    @inbounds for (i, dof) in enumerate(dofs_a)
        state.forcesⁿ⁺¹.Tᶜᵒⁿ[dof] += ta[i]
    end

    # Apply equal and opposite force to node 2
    @inbounds for (i, dof) in enumerate(dofs_b)
        state.forcesⁿ⁺¹.Tᶜᵒⁿ[dof] -= ta[i]
    end

    # Define indices in the sparsity map for matrix updates
    # mkindices: main diagonal or coupling block (self terms)
    # pkindices: cross terms (interaction between nodes)
    mkindices = @SVector [i + 12*(i-1) for i in [1,2,3,7,8,9]]
    pkindices = @SVector [i + 12*(mod1(i+6, 12)-1) for i in [7,8,9,1,2,3]]

    # Add contributions to stiffness (K) and damping (C) matrices
    dofs = constraint.global_sparsity_map[mkindices]
    state.matricesⁿ⁺¹.K.nzval[dofs] .+= k
    state.matricesⁿ⁺¹.C.nzval[dofs] .+= η * k

    # Subtract coupling terms (off-diagonal interactions)
    dofs = constraint.global_sparsity_map[pkindices]
    state.matricesⁿ⁺¹.K.nzval[dofs] .-= k
    state.matricesⁿ⁺¹.C.nzval[dofs] .-= η * k

end

#----------------------------------------------------
# Assembles all constraints into global system
#----------------------------------------------------
function assembly_constraints!(conf::BeamsConfiguration, state::SimulationState)

    @unpack nodes, constraints = conf

    # Reset constraint force vector before assembly
    fill!(state.forcesⁿ⁺¹.Tᶜᵒⁿ, 0)

    # Compute and assemble contributions from each constraint
    @batch for constraint in LazyRows(constraints)
        compute_constraint!(state, nodes, constraint)
    end 
    
end
