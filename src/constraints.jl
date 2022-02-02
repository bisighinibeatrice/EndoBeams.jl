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
    
    constraints = StructArray(constructor_constraint(connectivity[i][1], connectivity[i][2], stiffness, damping) for i in 1:length(connectivity))
    
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
function compute_multifreedom_constraints!(matrices, nodes, allconstraints, t)

    matrices.Tᶜᵒⁿ .= 0
    matrices.Kᶜᵒⁿ.nzval .= 0
    matrices.Cᶜᵒⁿ.nzval .= 0
    
    for c in allconstraints
        compute_constraint_contribution!(matrices, c, nodes)
    end 

end

function compute_constraint_contribution!(matrices, c, nodes)

    k = c.stiffness
    α = c.damping

    i1 = c.node1
    i2 = c.node2
    a = nodes[i1]
    b = nodes[i2]
    dofs_a = a.idof_6
    dofs_b = b.idof_6 

    xa = a.u
    xb = b.u
    va = a.u̇
    vb = b.u̇

    ta = k * (xb - xa) + α * k * (vb - va)

    update_vec(matrices.Tᶜᵒⁿ, dofs_a, ta)
    update_vec(matrices.Tᶜᵒⁿ, dofs_b, -ta)

    Kᶜᵒⁿ = Mat1212(
        -k, 0, 0, 0, 0, 0, k, 0, 0, 0, 0, 0,
        0, -k, 0, 0, 0, 0, 0, k, 0, 0, 0, 0,
        0, 0, -k, 0, 0, 0, 0, 0, k, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        k, 0, 0, 0, 0, 0, -k, 0, 0, 0, 0, 0,
        0, k, 0, 0, 0, 0, 0, -k, 0, 0, 0, 0,
        0, 0, k, 0, 0, 0, 0, 0, -k, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

    update_spmat(matrices.Kᶜᵒⁿ, c.sparsity_map, Kᶜᵒⁿ)
    update_spmat(matrices.Cᶜᵒⁿ, c.sparsity_map, α.*Kᶜᵒⁿ)

end 
