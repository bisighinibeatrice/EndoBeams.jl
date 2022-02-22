#----------------------------------
# STRUCTURE
#----------------------------------

struct MyConstraint{T}
    node1::Int # index of node 1   
    node2::Int # index of node 2
    stiffness::T
    damping::T
    sparsity_map::Vector{Int} 
end

#----------------------------------
# CONSTRUCTOR
#----------------------------------

"Constructor of the constraints StructArray"
function constructor_constraints(connectivity, stiffness, damping)
    
    constraints = StructArray(constructor_constraint(connectivity[i][1], connectivity[i][2], stiffness, damping) for i in eachindex(connectivity))
    
    return constraints
    
end 

# Constructor of one Constraint structure
function constructor_constraint(node1, node2, stiffness, damping, T=Float64)
    
    return MyConstraint{T}(node1, node2, stiffness, damping, zeros(Int, 144))
   
end 


#------------------------------------------
# FUNCTIONS 
#------------------------------------------

# Imposes multi freedom constraints at the current step
function constraints!(matrices, nodes, allconstraints)

    fill!(matrices.Tᶜᵒⁿ, 0)
    
    for c in LazyRows(allconstraints)

        k = c.stiffness
        α = c.damping

        i1 = c.node1
        i2 = c.node2
        dofs_a = nodes.idof_disp[i1]
        dofs_b = nodes.idof_disp[i2]

        xa = nodes.u[i1]
        xb = nodes.u[i2]
        va = nodes.u̇[i1]
        vb = nodes.u̇[i2]

        ta = k * (xb - xa) + α * k * (vb - va)

        @inbounds for (i, dof) in enumerate(dofs_a)
            matrices.Tᶜᵒⁿ[dof] += ta[i]
        end

        @inbounds for (i, dof) in enumerate(dofs_b)
            matrices.Tᶜᵒⁿ[dof] -= ta[i]
        end

        # Kᶜᵒⁿ = Mat1212(
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
        #     0,  0,  0,  0, 0, 0, 0,  0,  0,  0, 0, 0)

        mkindices = @SVector [i+12*(i-1) for i in [1,2,3,7,8,9]]
        pkindices = @SVector [i+12*(mod1(i+6, 12)-1) for i in [7,8,9,1,2,3]]


        @inbounds for dof in c.sparsity_map[mkindices]
            matrices.K[dof] += k
            matrices.C[dof] += α*k
        end

        @inbounds for dof in c.sparsity_map[pkindices]
            matrices.K[dof] -= k
            matrices.C[dof] -= α*k
        end

    end 

end



constraints!(matrices, nodes, allconstraints::Nothing) = nothing
