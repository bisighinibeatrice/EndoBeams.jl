#----------------------------------
# STRUCTURE
#----------------------------------

struct Constraint
    node1::Int # index of node 1   
    node2::Int # index of node 2
    k::Float64
    η::Float64
    sparsity_map::SVector{144, Int} 
end

#----------------------------------
# CONSTRUCTOR
#----------------------------------

"Constructor of the constraints StructArray"
function build_constraints(connectivity, k, η)
    
    constraints = StructArray(Constraint(connectivity[i, 1], connectivity[i, 2], k, η) for i in 1:size(connectivity, 1))
    
    return constraints
    
end 

# Constructor of one Constraint structure
function Constraint(node1, node2, k, η)
    
    return Constraint(node1, node2, k, η, zeros(Int, 144))
   
end 


#------------------------------------------
# FUNCTIONS 
#------------------------------------------

# Imposes multi freedom constraints at the current step
function constraints!(matrices, nodes, constraints)

    fill!(matrices.Tᶜᵒⁿ, 0)
    
    for c in LazyRows(constraints)

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

        # Kᶜᵒⁿ :
        #     -k, 0,  0,  0, 0, 0, k,  0,  0,  0, 0, 0,
        #     0,  -k, 0,  0, 0, 0, 0,  k,  0,  0, 0, 0,
        #     0,  0,  -k, 0, 0, 0, 0,  0,  k,  0, 0, 0,
        #     0,  0,  0,  0, 0, 0, 0,  0,  0,  0, 0, 0,
        #     0,  0,  0,  0, 0, 0, 0,  0,  0,  0, 0, 0,
        #     0,  0,  0,  0, 0, 0, 0,  0,  0,  0, 0, 0,
        #     k,  0,  0,  0, 0, 0, -k, 0,  0,  0, 0, 0,
        #     0,  k,  0,  0, 0, 0, 0,  -k, 0,  0, 0, 0,
        #     0,  0,  k,  0, 0, 0, 0,  0,  -k, 0, 0, 0,
        #     0,  0,  0,  0, 0, 0, 0,  0,  0,  0, 0, 0,
        #     0,  0,  0,  0, 0, 0, 0,  0,  0,  0, 0, 0,
        #     0,  0,  0,  0, 0, 0, 0,  0,  0,  0, 0, 0

        mkindices = @SVector [i+12*(i-1) for i in [1,2,3,7,8,9]]
        pkindices = @SVector [i+12*(mod1(i+6, 12)-1) for i in [7,8,9,1,2,3]]


        dofs = c.sparsity_map[mkindices]
        matrices.K[dofs] .+= k
        matrices.C[dofs] .+= η*k

        dofs = c.sparsity_map[pkindices]
        matrices.K[dofs] .-= k
        matrices.C[dofs] .-= η*k

    end 

end



constraints!(matrices, nodes, constraints::Nothing) = nothing
