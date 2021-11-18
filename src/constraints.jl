#----------------------------------
# STRUCTURE
#----------------------------------

struct MyConstraint{T}
    node1::Int # index of node 1   
    node2::Int # index of node 2
    stiffness::T
    damping::T
    sparsity_map::Array{Int,1} 
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
function compute_multifreedom_constraints!(matrices, allnodes, allconstraints, t)

    matrices.Tconstr .= 0
    matrices.Kconstr.nzval .= 0
    matrices.Cconstr.nzval .= 0
    
    for c in allconstraints
        compute_constraint_contribution!(matrices, c, allnodes)
    end 

end

function compute_constraint_contribution!(matrices, c, allnodes)

    k = c.stiffness
    alpha = c.damping

    i1 = c.node1
    i2 = c.node2
    a = allnodes[i1]
    b = allnodes[i2]
    dofs_a = a.idof_6
    dofs_b = b.idof_6 

    xa = a.u
    xb = b.u
    va = a.udt
    vb = b.udt

    ta = k * (xb - xa) + alpha * k * (vb - va)

    update_vec(matrices.Tconstr, dofs_a, ta)
    update_vec(matrices.Tconstr, dofs_b, -ta)

    Kconstr = Mat1212(
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

    update_spmat(matrices.Kconstr, c.sparsity_map, Kconstr)
    update_spmat(matrices.Cconstr, c.sparsity_map, alpha.*Kconstr)

end 
