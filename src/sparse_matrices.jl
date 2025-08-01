#--------------------------------------------------------------------------------------------------
# FUNCTIONS TO COMPUTE SPARSITY PATTERN FOR BEAM ELEMENTS (INCLUDING CONSTRAINTS)
#--------------------------------------------------------------------------------------------------

"""
    sparsity_beams(nodes, beams, constraints, fixed_dofs)

Computes the sparsity pattern (non-zero structure) of the global stiffness matrix for a beam system
that may include constraints between nodes.

Returns:
- `I`, `J`: Row and column indices for the sparse matrix.
- `sparsity_free`: Linear indices of free (non-fixed) degrees of freedom in the global matrix.
"""
function sparsity_beams(nodes, beams, constraints, fixed_dofs)

    beam_ndofs = 12                        # DOFs per beam element (6 per node × 2 nodes)
    N = length(beams) * beam_ndofs^2       # Estimate initial non-zero count from beams

    # Account for additional non-zero entries from constraints
    if !isnothing(constraints)
        N += length(constraints) * beam_ndofs^2
    end

    # Preallocate sparsity pattern arrays
    I = Vector{Int}(undef, N)                              # Row indices
    J = Vector{Int}(undef, N)                              # Column indices
    infos = Vector{Vector{Vec4{Int}}}(undef, N)            # Metadata for each non-zero entry

    k = 1  # Sparsity index counter

    # ----------------------------
    # Process beam elements
    # ----------------------------
    for bi in eachindex(beams)
        n1, n2 = beams.node1[bi], beams.node2[bi]                   # Nodes of beam
        dofs = vcat(nodes.local_dofs[n1], nodes.local_dofs[n2])    # DOFs for both nodes
        local_dofs = 1                                              # Local DOF index

        for j in dofs, i in dofs
            I[k] = i
            J[k] = j
            # Store beam info: (type = 1, beam index, local DOF, is_fixed flag)
            infos[k] = [Vec4(1, bi, local_dofs, i in fixed_dofs || j in fixed_dofs)]
            k += 1
            local_dofs += 1
        end
    end

    # ----------------------------
    # Process constraints (if present)
    # ----------------------------
    if !isnothing(constraints)
        for ci in eachindex(constraints)
            n1 = constraints.node1[ci]
            n2 = constraints.node2[ci]
            dofs = vcat(nodes.local_dofs[n1], nodes.local_dofs[n2])  # 6 DOFs per node
            local_dof = 1

            for j in dofs, i in dofs
                I[k] = i
                J[k] = j
                # Store constraint info: (type = 2, constraint index, local DOF, is_fixed flag)
                infos[k] = [Vec4(2, ci, local_dof, i in fixed_dofs || j in fixed_dofs)]
                k += 1
                local_dof += 1
            end
        end
    end

    # Create a sparse matrix holding metadata (Vec4 info for each entry)
    K = sparse(I, J, infos, maximum(I), maximum(J), vcat)

    # ----------------------------
    # Identify and collect free DOF indices
    # ----------------------------
    sparsity_free = Vector{Int}(undef, length(K.nzval))
    ksf = 0  # Counter for free DOFs

    # Create sparsity maps for beams and constraints
    beams_spmap = [MVector{144, Int}(undef) for _ in 1:length(beams)]
    if !isnothing(constraints)
        constraints_spmap = [MVector{144, Int}(undef) for _ in 1:length(constraints)]
    end

    for (i, infos) in enumerate(K.nzval)
        for (type, idx, local_dof, _) in infos
            if type == 1
                beams_spmap[idx][local_dof] = i
            elseif !isnothing(constraints) && type == 2
                constraints_spmap[idx][local_dof] = i
            end
        end

        # Mark as free if neither DOF is fixed
        if infos[1][4] == 0
            ksf += 1
            sparsity_free[ksf] = i
        end
    end

    # Trim unused entries from sparsity_free
    resize!(sparsity_free, ksf)

    # Assign sparsity maps to beam and constraint objects
    for i in eachindex(beams)
        beams.local_sparsity_map[i] = beams_spmap[i]
        beams.global_sparsity_map[i] = beams_spmap[i]
    end

    if !isnothing(constraints)
        for i in eachindex(constraints)
            constraints.local_sparsity_map[i] = constraints_spmap[i]
            constraints.global_sparsity_map[i] = constraints_spmap[i]
        end
    end

    return I, J, sparsity_free
end

"""
    sparse_matrices_beams!(conf::BeamsConfiguration)

Constructs the sparse matrix structure used for stiffness and damping in the beam system.

Returns:
- `matrices`: A custom `Matrices` structure holding the global sparsity pattern.
- `sol`: A `Solution` object initialized with the global and reduced stiffness matrices.
"""
function sparse_matrices_beams!(conf::BeamsConfiguration)

    @unpack beams, nodes, constraints, bcs, ndofs = conf

    free_dofs = bcs.free_dofs        # Indices of non-fixed DOFs
    fixed_dofs = bcs.fixed_dofs      # Indices of fixed DOFs
    nfreedofs = length(free_dofs)    # Total number of free DOFs

    # Compute sparsity pattern (I, J) and free entry indices
    I, J, sparsity_free = sparsity_beams(nodes, beams, constraints, fixed_dofs)

    # Initialize global stiffness matrix with structural sparsity
    Ktan = sparse(I, J, 0.0)

    # Extract reduced system matrix (free DOFs only)
    Ktan_free = Ktan[free_dofs, free_dofs]

    # Wrap up matrices and solution structure for simulation
    matrices = Matrices(I, J, sparsity_free)
    sol = Solution(Ktan, Ktan_free, ndofs, nfreedofs)

    return matrices, sol
end
