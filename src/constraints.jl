#----------------------------------
# STRUCTURE
#----------------------------------

struct Constraint
    node1::Int # index of node 1   
    node2::Int # index of node 2
    k::Float64
    η::Float64
    sparsity_map::SVector{144, Int} 
    Req::Float64
end



#----------------------------------
# CONSTRUCTOR
#----------------------------------

"Constructor of the constraints StructArray"
function build_constraints(nodespairs, k, η, Req=0)
    
    constraints = StructArray(Constraint(nodespairs[i][1], nodespairs[i][2], k isa AbstractVector ? k[i] : k, η isa AbstractVector ? η[i] : η, Req) for i in 1:size(nodespairs, 1))
    
    return constraints
    
end 

function build_constraints(nodespairs::Matrix, k, η, Req=0)
    
    constraints = StructArray(Constraint(nodespairs[i, 1], nodespairs[i, 2], k isa AbstractVector ? k[i] : k, η isa AbstractVector ? η[i] : η, Req) for i in 1:size(nodespairs, 1))
    
    return constraints
    
end 


# Constructor of one Constraint structure
function Constraint(node1, node2, k, η, Req)
    
    return Constraint(node1, node2, k, η, zeros(Int, 144), Req)
   
end 

#------------------------------------------
# FUNCTIONS 
#------------------------------------------

# Imposes multi freedom constraints at the current step
function constraint_zero_distance!(matrices, nodes, c)

    k = c.k
    η = c.η

    i1 = c.node1
    i2 = c.node2
    dofs_a = nodes.idof_disp[i1]
    dofs_b = nodes.idof_disp[i2]

    xa = nodes.u[i1]
    xb = nodes.u[i2]
    va = nodes.u̇[i1]
    vb = nodes.u̇[i2]

    ta = k * (xb - xa) + η * k * (vb - va)

    @inbounds for (i, dof) in enumerate(dofs_a)
        matrices.Tᶜᵒⁿ[dof] += ta[i]
    end

    @inbounds for (i, dof) in enumerate(dofs_b)
        matrices.Tᶜᵒⁿ[dof] -= ta[i]
    end

    mkindices = @SVector [i+12*(i-1) for i in [1,2,3,7,8,9]]
    pkindices = @SVector [i+12*(mod1(i+6, 12)-1) for i in [7,8,9,1,2,3]]

    dofs = c.sparsity_map[mkindices]
    matrices.K[dofs] .+= k
    matrices.C[dofs] .+= η*k

    dofs = c.sparsity_map[pkindices]
    matrices.K[dofs] .-= k
    matrices.C[dofs] .-= η*k

end

function constraint_fixed_distance!(matrices, nodes, c, t)

    k = c.k 
    η = c.η
    r0 = c.Req

    i1 = c.node1
    i2 = c.node2
    dofs_a = nodes.idof_disp[i1]
    dofs_b = nodes.idof_disp[i2]

    xa = nodes.u[i1]
    xb = nodes.u[i2]
    va = nodes.u̇[i1]
    vb = nodes.u̇[i2]

    ra = xb - xa
    rb = xa - xb
    rnorm = norm(ra)
    rnorm3 = rnorm^3
    kr = k*(1-r0/rnorm)
    Ka = k*r0/rnorm3 * (ra * ra')
    Kb = k*r0/rnorm3 * (rb * rb')

    ta = kr * (rb - ra) + η * kr * (vb - va)

    @inbounds for (i, dof) in enumerate(dofs_a)
        matrices.Tᶜᵒⁿ[dof] += ta[i]
    end

    @inbounds for (i, dof) in enumerate(dofs_b)
        matrices.Tᶜᵒⁿ[dof] -= ta[i]
    end

    Kconstr =  [
        -kr-Ka[1,1], -Ka[2,1], -Ka[3,1], 0, 0, 0, kr+Kb[1,1], Kb[2,1], Kb[3,1], 0, 0, 0,
        -Ka[1,2], -kr-Ka[2,2], -Ka[3,2], 0, 0, 0, Kb[1,2], kr+Kb[2,2], Kb[3,2], 0, 0, 0,
        -Ka[1,3], -Ka[2,3], -kr-Ka[3,3], 0, 0, 0, Kb[1,3], Kb[2,3], kr+Kb[3,3], 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        kr+Ka[1,1], Ka[2,1], Ka[3,1], 0, 0, 0, -kr-Kb[1,1], -Kb[2,1], -Kb[3,1], 0, 0, 0,
        Ka[1,2], kr+Ka[2,2], Ka[3,2], 0, 0, 0, -Kb[1,2], -kr-Kb[2,2], -Kb[3,2], 0, 0, 0,
        Ka[1,3], Ka[2,3], kr+Ka[3,3], 0, 0, 0, -Kb[1,3], -Kb[2,3], -kr-Kb[3,3], 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    matrices.K[c.sparsity_map] .-= Kconstr
    matrices.C[c.sparsity_map] .-= η*Kconstr

end

function constraints!(matrices, nodes, constraints, tⁿ⁺¹)

    fill!(matrices.Tᶜᵒⁿ, 0)

    for c in LazyRows(constraints)

        if c.Req == 0 
            constraint_zero_distance!(matrices, nodes, c)
        elseif c.Req !=0
            constraint_fixed_distance!(matrices, nodes, c, tⁿ⁺¹)
        end 

    end 

end

# constraints!(matrices, nodes, constraints::Nothing) = nothing
# constraint!(matrices, nodes, c::LazyRow{Constraint}) = constraint_zero_distance!(matrices, nodes, c)
# constraint!(matrices, nodes, c::LazyRow{ConstraintFixedDistance}) = constraint_fixed_distance!(matrices, nodes, c)

