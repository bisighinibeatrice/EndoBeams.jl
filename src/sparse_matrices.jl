#--------------------------------------------------------------------------------------------------
# FUNCTIONS TO COMPUTE SPARSITY PATTERN FOR BEAM ELEMENTS
#------------------------------------------------------------------------------------------------------

# Function to compute sparsity pattern for beam elements
function sparsity_beams(nodes, beams, fixed_dofs)

    beam_ndofs = 12  # Number of dofs per beam element
    N = length(beams) * beam_ndofs^2  # Total non-zero entries in sparsity pattern for beams

    # # Add constraints to the sparsity pattern if any
    # if !isnothing(constraints)
    #     N += length(constraints) * beam_ndofs^2
    # end

    # Initialize sparsity pattern
    I = Vector{Int}(undef, N)    # Row indices for the sparse matrix
    J = Vector{Int}(undef, N)    # Column indices for the sparse matrix
    infos = Vector{Vector{Vec4{Int}}}(undef, N)  # Store information about non-zero entries

    k = 1  # Index for sparsity pattern
    # Loop through all beams to calculate the sparsity pattern
    for bi in eachindex(beams)
        n1, n2 = beams.node1[bi], beams.node2[bi]  # Get node indices for the beam
        dofs = vcat(nodes.local_dofs[n1], nodes.local_dofs[n2])  # Get dofs for both nodes
        local_dofs = 1  # Local DOF index for the beam element

        # Loop through all combinations of dofs for the beam
        for j in dofs, i in dofs
            I[k] = i
            J[k] = j
            infos[k] = [Vec4(1, bi, local_dofs, i in fixed_dofs || j in fixed_dofs)]  # Store info (type=1 for beam)
            k += 1
            local_dofs += 1
        end
    end

    # Create the sparse matrix using the computed sparsity pattern
    K = sparse(I, J, infos, maximum(I), maximum(J), vcat)

    # Initialize sparsity free vector to store free dofs
    sparsity_free = Vector{Int}(undef, length(K.nzval))
    ksf = 0  # Counter for free dofs

    # Create maps for beams and constraints sparsity
    beams_spmap = [MVector{144, Int}(undef) for _ in 1:length(beams)]
    # constraints_spmap = isnothing(constraints) ? [] : [MVector{144, Int}(undef) for _ in 1:length(constraints)]

    # Loop through non-zero entries in the sparse matrix
    for (i, infos) in enumerate(K.nzval)
        for (type, idx, local_dofs) in infos
            if type == 1
                beams_spmap[idx][local_dofs] = i  # For beams
            end
        end
        # Identify free dofs
        if infos[1][4] == 0
            ksf += 1
            sparsity_free[ksf] = i
        end
    end

    # Resize sparsity_free vector to remove unused entries
    resize!(sparsity_free, ksf)

    # Assign the sparsity map to each beam and constraint (if present)
    max_dofs = 0
    for i in eachindex(beams)
        beams.local_sparsity_map[i] = beams_spmap[i]
        beams.global_sparsity_map[i] = beams_spmap[i]
        if maximum(beams.global_sparsity_map[i])>max_dofs
            max_dofs = maximum(beams.global_sparsity_map[i])
        end
    end

    return I, J, sparsity_free, max_dofs
end

# Function to construct sparse matrices for beam configuration
function sparse_matrices_beams!(conf::BeamsConfiguration)

    @unpack beams, nodes, bcs, ndofs = conf

    free_dofs = bcs.free_dofs  # Free dofs from boundary conditions
    fixed_dofs = bcs.fixed_dofs  # Fixed dofs from boundary conditions
    nfreedofs = length(free_dofs)  # Number of free dofs

    # Compute sparsity pattern for beams and constraints (if present)
    I, J, sparsity_free = sparsity_beams(nodes, beams, fixed_dofs)

    # Create the tangent stiffness matrix (initially empty)
    Ktan = sparse(I, J, 0.0)

    # Extract the sub-matrix for free dofs
    Ktan_free = Ktan[free_dofs, free_dofs]

    # Create matrices structure for further use
    matrices = Matrices(I, J, sparsity_free)

    # Create nodal solution to hold the stiffness matrix and other info
    sol = Solution(Ktan, Ktan_free, ndofs, nfreedofs)

    return matrices, sol
end
