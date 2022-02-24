#------------------------------------------
# FUNCTIONS TO COMPUTE SPARSITY PATTERN
#------------------------------------------

# Compute a structure where we associated to each node the nodes that share a common beam
compute_connected_nodes(constraints, beams, nnodes) = [connections(n, beams, constraints) for n in 1:nnodes]


function connections(n, beams, constraints)

    connected_nodes_i = Int[]
    
    for b in LazyRows(beams) 
        
        i1 = b.node1
        i2 = b.node2
        
        if i1 == n 
            push!(connected_nodes_i, i2)
        elseif i2 == n
            push!(connected_nodes_i, i1)
        end
        
    end
    
    if !isnothing(constraints)
        for c in LazyRows(constraints)
            
            i1 = c.node1
            i2 = c.node2
            
            if i1 == n 
                push!(connected_nodes_i, i2)
            elseif i2 == n
                push!(connected_nodes_i, i1)
                
            end
        end 
    end
    
    sort!(connected_nodes_i)

    return connected_nodes_i
    
end 

# Compute the sparsity map for the beams
function compute_beams_sparsity_map!(beams, nodes, N, cols, rows)
    
    for b in LazyRows(beams)
        
        i1 = b.node1
        i2 = b.node2
        idof1 = nodes.idof_6[i1]
        idof2 = nodes.idof_6[i2]
        
        if i1 > i2
            error("Error in the structure connectivity: i1 > i2")
        end 
        
        idof = vcat(idof1, idof2)
        
        k = 0
        for n in 1:N
            
            col = cols[n]
            row = rows[n]
            
            if col in idof && row in idof
                
                k += 1
                b.sparsity_map[k] = n
                
            end 
        end     
        
    end 
    
end 

# Compute the sparsity map for the constraints
function compute_constraints_sparsity_map!(constraints, nodes, N, cols, rows)
    
    for c in constraints
        
        i1 = c.node1
        i2 = c.node2
        idof1 = nodes.idof_6[i1]
        idof2 = nodes.idof_6[i2]
        
        if i1 > i2
            error("Error in the structure connectivity: i1 > i2")
        end 
        
        idof = vcat(idof1, idof2)
        
        k = 0
        for n in 1:N
            
            col = cols[n]
            row = rows[n]
            
            if col in idof && row in idof
                
                k += 1
                c.sparsity_map[k] = n
                
            end 
        end     
        
    end 
    
end 

compute_constraints_sparsity_map!(constraints::Nothing, nodes, N, cols, rows) = nothing

#------------------------------------------
# FUNCTIONS TO COMPUTE SPARSITY PATTERN
#------------------------------------------

# Compute the sparsity pattern for the tangent matrix
function compute_sparsity_pattern_tan!(ndofs_per_node, nnodes, connected_nodes)
    
    mat = zeros(Bool, nnodes*ndofs_per_node, nnodes*ndofs_per_node)
    
    for n in 1:nnodes
        idofⁿ = 6*(n-1) .+ @SVector [1,2,3,4,5,6]
        for j in idofⁿ
            for cn in connected_nodes[n] 
                idof_cn = 6*(cn-1) .+ @SVector [1,2,3,4,5,6]
                for i in idof_cn
                    mat[i,j] = true
                end
            end 
        end 
    end 
    
    nrows = size(mat,1)
    ncols = size(mat,2)
    rows = Int[]
    cols = Int[]
    
    for i in 1:nrows 
        for j in 1:ncols
            
            if mat[i,j] == true
                push!(rows, i)
                push!(cols, j)
            end
        end    
    end 
    
    return (rows, cols, ones(size(rows)))
    
end 

# Compute the sparsity pattern for the free total matrix
function compute_sparsity_pattern_free!(ntot, connected_nodes, free_dofs_all, nnodes, T=Float64)
    
    mat = zeros(Int, ntot,ntot)

    for n in 1:nnodes
        idofⁿ = 6*(n-1) .+ @SVector [1,2,3,4,5,6]
        for j in idofⁿ
            for cn in connected_nodes[n] 
                idof_cn = 6*(cn-1) .+ @SVector [1,2,3,4,5,6]
                for i in idof_cn
                    mat[i,j] = 2
                end
            end 
        end 
    end 

    #-------------------------------------------------

    nrows = size(mat,1)
    ncols = size(mat,2)
    rows = Int[]
    cols = Int[]
    zvals = T[]

    for i in 1:nrows 
        for j in 1:ncols
            if mat[i,j] != 0
                push!(rows, i)
                push!(cols, j)
                push!(zvals, mat[i,j])
            end
        end    
    end 

    #-------------------------------------------------
    
    mat_free = mat[free_dofs_all, free_dofs_all]

    nrows_free = size(mat_free,1)
    ncols_free = size(mat_free,2)
    rows_free = Int[]
    cols_free = Int[]
    zvals_free  = T[]

    k = 0
    for i in 1:nrows_free 
        for j in 1:ncols_free
            k+=1
            if mat_free[i,j] != 0
                push!(rows_free, i)
                push!(cols_free, j)
                push!(zvals_free, mat_free[i,j])
            end
            
        end    
    end 

    if size(sparse(rows_free, rows_free, zvals_free)) != size(mat_free)
        @error "SparseMatrixError: error in the construction of the free sparse matrix"
    end 

    #-------------------------------------------------
    
    spmap_free = Int[]
    k = 0 
    for n in 1:length(rows)
        i = rows[n]
        j = cols[n]
        k+=1
        if mat[i,j] != 0 && i in free_dofs_all && j in free_dofs_all
            push!(spmap_free, k)
        end
    end 
       
    return (rows_free, cols_free, zvals_free), spmap_free
    
end 

# Compute the sparsity pattern which allows to preallocate the sparse matrices and the sparsity map of the beams and constraints (for the assemby)
function compute_sparsity!(beams, nodes, constraints, conf)
    


    @infiltrate

    # compute linked nodes from connectivity 
    connected_nodes = compute_connected_nodes(constraints, beams, nnodes)
    
    # tangent matrix
    (rows_tan, cols_tan, zval_tan) = compute_sparsity_pattern_tan!(ndofs_per_node, nnodes, connected_nodes)
    n_tan = length(rows_tan) 
    
    # free total matrix
    (rows_free, cols_free, zvals_free), spmap_free = compute_sparsity_pattern_free!(ntot, connected_nodes, free_dofs, nnodes, T)

    # compute beam sparsity map 
    compute_beams_sparsity_map!(beams, nodes, n_tan, cols_tan, rows_tan)
    compute_constraints_sparsity_map!(constraints, nodes, n_tan, cols_tan, rows_tan)
    
    return (rows_tan, cols_tan, zval_tan), (rows_free, cols_free, zvals_free), spmap_free 
    
end 


function sparsity(nodes, beams, constraints, fixed_dofs)

    beam_ndofs = 12
    N = length(beams)*beam_ndofs^2
    if !isnothing(constraints)
        N += length(constraints)*beam_ndofs^2
    end
    I = Vector{Int}(undef, N)
    J = Vector{Int}(undef, N)
    infos = Vector{Vector{Vec4{Int}}}(undef, N)
    k = 1
    for bi in eachindex(beams)
        n1 = beams.node1[bi]
        n2 = beams.node2[bi]
        dofs = vcat(nodes.idof_6[n1], nodes.idof_6[n2])
        local_dof = 1
        for j in dofs, i in dofs
            I[k] = i
            J[k] = j
            infos[k] = [Vec4(1, bi, local_dof, i in fixed_dofs || j in fixed_dofs)]
            k += 1
            local_dof += 1
        end
    end

    if !isnothing(constraints)
        for ci in eachindex(constraints)
            n1 = constraints.node1[ci]
            n2 = constraints.node2[ci]
            dofs = vcat(nodes.idof_6[n1], nodes.idof_6[n2])
            local_dof = 1
            for j in dofs, i in dofs
                I[k] = i
                J[k] = j
                infos[k] = [Vec4(2, ci, local_dof, i in fixed_dofs || j in fixed_dofs)]
                k += 1
                local_dof += 1
            end
        end
    end


    K = sparse(I, J, infos, maximum(I), maximum(J), vcat)

    sparsity_free = Vector{Int}(undef, length(K.nzval))

    ksf = 0

    for (i, infos) in enumerate(K.nzval)
        for (type, idx, local_dof) in infos
            if type==1
                beams.sparsity_map[idx][local_dof] = i
            elseif type==2
                constraints.sparsity_map[idx][local_dof] = i
            end
        end
        if infos[1][4] == 0
            ksf += 1
            sparsity_free[ksf] = i
        end
    end

    resize!(sparsity_free, ksf)

    return I, J, sparsity_free
end


