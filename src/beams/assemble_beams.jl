function assemble!(conf::BeamsConfiguration, state::SimulationState, params::SimulationParams, only_beams=true) 

    # Initialize the force vectors and stiffness matrix for the entire structure
    if only_beams
        fill!(state.matricesⁿ⁺¹.K.nzval, 0)              # Global stiffness matrix
        fill!(state.matricesⁿ⁺¹.C.nzval, 0)              # Global damping matrix
        fill!(state.matricesⁿ⁺¹.M.nzval, 0)              # Global mass matrix
        fill!(state.forcesⁿ⁺¹.Tⁱⁿᵗ, 0)         # Internal force vector
        fill!(state.forcesⁿ⁺¹.Tᵏ, 0)           # Kinetic force vector
        fill!(state.forcesⁿ⁺¹.Tᶜ, 0)           # Kinetic force vector
    end

    # Unpack configuration to retrieve node and beam information
    @unpack nodes, beams = conf

    # Unpack information for Gauss points
    @unpack nᴳ_beams, ωᴳ_beams, zᴳ_beams = params
    gauss_params  = (nᴳ_beams, ωᴳ_beams, zᴳ_beams)   

    # Loop over all beam beams to compute their contributions
    @batch for beam in LazyRows(beams)

        # Retrieve node indices for the current beam beam
        n1 = beam.node1
        n2 = beam.node2

        # Extract nodal data for the two nodes connected by the beam
        X₁, X₂ = nodes.X₀[n1], nodes.X₀[n2]             # Initial positions
        u₁, u₂ = nodes.u[n1], nodes.u[n2]               # Displacements
        u̇₁, u̇₂ = nodes.u̇[n1], nodes.u̇[n2]               # Velocities
        ü₁, ü₂ = nodes.ü[n1], nodes.ü[n2]           # Accelerations
        ẇ₁, ẇ₂ = nodes.ẇ[n1], nodes.ẇ[n2]           # Rotational velocities
        ẅ₁, ẅ₂ = nodes.ẅ[n1], nodes.ẅ[n2]           # Rotational accelerations
        R₁, R₂ = nodes.R[n1], nodes.R[n2]               # Rotational transformations
        ΔR₁, ΔR₂ = nodes.ΔR[n1], nodes.ΔR[n2]           # Incremental rotations
    
        # Pack initialization constants for the beam computation
        init = (X₁, X₂, beam.l₀, beam.Rₑ⁰)  # Beam index, geometry, and initial rotation
        constants = (init, gauss_params, beam.properties)

        # Compute beam-level contributions
        strain_energy, kinetic_energy, 
        Tⁱⁿᵗ, Tᵏ, Kⁱⁿᵗ, M, Cᵏ = compute_beams(
            u₁, u₂, R₁, R₂, ΔR₁, ΔR₂, 
            u̇₁, u̇₂, ẇ₁, ẇ₂, 
            ü₁, ü₂, ẅ₁, ẅ₂, 
            constants
        )
        
        # Assemble contributions into global vectors and matrices
        dofs1 = nodes.global_dofs[n1]               # dofs for node 1
        dofs2 = nodes.global_dofs[n2]               # dofs for node 2
        dofs = vcat(dofs1, dofs2)              # Combine dofs for the beam

        # Accumulate beam-level energies
        state.energyⁿ⁺¹.strain_energy += strain_energy
        state.energyⁿ⁺¹.kinetic_energy += kinetic_energy

        # Add beam-level forces to global force vectors
        state.forcesⁿ⁺¹.Tᵏ[dofs] += Tᵏ               # Kinetic forces
        state.forcesⁿ⁺¹.Tⁱⁿᵗ[dofs] += Tⁱⁿᵗ           # Internal forces

        # Add beam-level matrices to global matrices
        state.matricesⁿ⁺¹.K.nzval[beam.global_sparsity_map] += vec(Kⁱⁿᵗ)  # Stiffness matrix
        state.matricesⁿ⁺¹.C.nzval[beam.global_sparsity_map] += vec(Cᵏ)  # Damping matrix
        state.matricesⁿ⁺¹.M.nzval[beam.global_sparsity_map] += vec(M)  # Mass matrix

    end

end
