#------------------------------------------
# FUNCTIONS TO COMPUTE SPARSITY PATTERN
#------------------------------------------


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

    beams_spmap = [MVector{144, Int}(undef) for _ in 1:length(beams)]
    if !isnothing(constraints)
        constraints_spmap = [MVector{144, Int}(undef) for _ in 1:length(constraints)]
    end

    for (i, infos) in enumerate(K.nzval)
        for (type, idx, local_dof) in infos
            if type==1
                beams_spmap[idx][local_dof] = i
            elseif type==2 && !isnothing(constraints)
                constraints_spmap[idx][local_dof] = i
            end
        end
        if infos[1][4] == 0
            ksf += 1
            sparsity_free[ksf] = i
        end
    end

    resize!(sparsity_free, ksf)

    for i in eachindex(beams)
        beams.sparsity_map[i] = beams_spmap[i]
    end
    if !isnothing(constraints)
        for i in eachindex(constraints)
            constraints.sparsity_map[i] = constraints_spmap[i]
        end
    end

    return I, J, sparsity_free
end


