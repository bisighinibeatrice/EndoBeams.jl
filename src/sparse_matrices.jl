#------------------------------------------
# FUNCTIONS TO COMPUTE SPARSITY PATTERN
#------------------------------------------

# Compute a structure where we associated to each node the nodes that share a common beam
function compute_connected_nodes!(connected_nodes, allconstraints, allbeams, nnodes)
    
    for n in 1:nnodes 
        
        connected_nodes_i = Int[]
        
        push!(connected_nodes_i, n)
        
        for b in allbeams 
            
            i1 = b.node1
            i2 = b.node2
            
            if i1 == n 
                
                push!(connected_nodes_i, i2)
                
            elseif i2 == n
                
                push!(connected_nodes_i, i1)
                
            end
            
        end
        
        for c in allconstraints
            
            i1 = c.node1
            i2 = c.node2
            
            if i1 == n 
                
                push!(connected_nodes_i, i2)
                
            elseif i2 == n
                
                push!(connected_nodes_i, i1)
                
            end
        end 
        
        sort!(connected_nodes_i)
        push!(connected_nodes, connected_nodes_i)
        
    end 
    
end

# Compute the sparsity map for the beams
function compute_beams_sparsity_map!(allbeams, allnodes, N, cols, rows)
    
    for b in allbeams
        
        i1 = b.node1
        i2 = b.node2
        idof1 = allnodes.idof_6[i1]
        idof2 = allnodes.idof_6[i2]
        
        if i1 > i2
            error("Error in the structure connectivity: i1 > i2")
        end 
        
        idof = vcat(idof1, idof2)
        
        cnt = 0
        for n in 1:N
            
            col = cols[n]
            row = rows[n]
            
            if col in idof && row in idof
                
                cnt += 1
                b.sparsity_map[cnt] = n
                
            end 
        end     
        
    end 
    
end 

# Compute the sparsity map for the constraints
function compute_constraints_sparsity_map!(allconstraints, allnodes, N, cols, rows)
    
    for c in allconstraints
        
        i1 = c.node1
        i2 = c.node2
        idof1 = allnodes.idof_6[i1]
        idof2 = allnodes.idof_6[i2]
        
        if i1 > i2
            error("Error in the structure connectivity: i1 > i2")
        end 
        
        idof = vcat(idof1, idof2)
        
        cnt = 0
        for n in 1:N
            
            col = cols[n]
            row = rows[n]
            
            if col in idof && row in idof
                
                cnt += 1
                c.sparsity_map[cnt] = n
                
            end 
        end     
        
    end 
    
end 

#------------------------------------------
# FUNCTIONS TO COMPUTE SPARSITY PATTERN
#------------------------------------------

# Compute the sparsity pattern for the tangent matrix
function compute_sparsity_pattern_tan!(ndofs_per_node, nnodes, connected_nodes)
    
    mat = zeros(Bool, nnodes*ndofs_per_node, nnodes*ndofs_per_node)
    
    for n in 1:nnodes
        idof_n = 6*(n-1) .+ [1,2,3,4,5,6]
        for j in idof_n
            for cn in connected_nodes[n] 
                idof_cn = 6*(cn-1) .+ [1,2,3,4,5,6]
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
    
    mat = zeros(ntot,ntot)

    for n in 1:nnodes
        idof_n = 6*(n-1) .+ [1,2,3,4,5,6]
        for j in idof_n
            for cn in connected_nodes[n] 
                idof_cn = 6*(cn-1) .+ [1,2,3,4,5,6]
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

    cnt = 0
    for i in 1:nrows_free 
        for j in 1:ncols_free
            cnt+=1
            if mat_free[i,j] != 0
                push!(rows_free, i)
                push!(cols_free, j)
                push!(zvals_free, mat_free[i,j])
            end
            
        end    
    end 

    if size(sparse(rows_free, rows_free, zvals_free)) != size(mat_free)
        @debug "SparseMatrixError: error in the construction of the free sparse matrix"
    end 

    #-------------------------------------------------
    
    spmap_free = Int[]
    cnt = 0 
    for n in 1:length(rows)
        i = rows[n]
        j = cols[n]
        cnt+=1
        if mat[i,j] != 0 && i in free_dofs_all && j in free_dofs_all
            push!(spmap_free, cnt)
        end
    end 
       
    return (rows_free, cols_free, zvals_free), spmap_free
    
end 

# Compute the sparsity pattern which allows to preallocate the sparse matrices and the sparsity map of the beams and constraints (for the assemby)
function compute_sparsity!(allbeams, allnodes, pncons, conf, T=Float64)
    
    # initialization
    free_dofs = conf.bc.free_dofs

    nnodes = length(allnodes)
    ndofs_per_node = 6
    ndofs = conf.ndofs 
    ntot = ndofs 

    # compute linked nodes from connectivity 
    connected_nodes = Vector{Vector{Int}}() #TODO: chack if there's another way to do this with less allocations
    compute_connected_nodes!(connected_nodes, pncons, allbeams, nnodes)
    
    # tangent matrix
    (rows_tan, cols_tan, zval_tan) = compute_sparsity_pattern_tan!(ndofs_per_node, nnodes, connected_nodes)
    n_tan = length(rows_tan) 
    
    # free total matrix
    (rows_free, cols_free, zvals_free), spmap_free = compute_sparsity_pattern_free!(ntot, connected_nodes, free_dofs, nnodes, T)

    # compute beam sparsity map 
    compute_beams_sparsity_map!(allbeams, allnodes, n_tan, cols_tan, rows_tan)
    compute_constraints_sparsity_map!(pncons, allnodes, n_tan, cols_tan, rows_tan)
    
    return (rows_tan, cols_tan, zval_tan), (rows_free, cols_free, zvals_free), spmap_free 
    
end 

#------------------------------------------------
# FUNCTIONS TO UPDATE SPARSE MATRICES AND VECTORS
#------------------------------------------------

# Update nzval of a global sparse matrix given a local matrix (sum)
function update_spmat_sum(sp, sparsity_map, newz) 
    
    for (k, i) in enumerate(sparsity_map)
        sp.nzval[i] += newz[k]
    end
    
end 

# Update nzval of a global sparse matrix given a local matrix 
function update_spmat(sp, sparsity_map, newz) 
    
    for (k, i) in enumerate(sparsity_map)
        sp.nzval[i] = newz[k]
    end
    
end 

# Update a sparse vector given a local vector (sum)
function update_vec_sum(sp, idof, newK) 
    
    @inbounds for (k, i) in enumerate(idof)
        sp[i] += newK[k]     
    end 
    
end 

# Update a sparse vector given a local vector 
function update_vec(sp, idof, newK) 
    
    @inbounds for (k, i) in enumerate(idof)
        sp[i] = newK[k]     
    end 
    
end 