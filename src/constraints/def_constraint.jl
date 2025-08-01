#----------------------------------
# CONSTRAINT STRUCTURE DEFINITION
#----------------------------------

"""
    Constraint

A structure that defines a constraint between two nodes with specified stiffness (`k`), damping (`η`),  and a sparsity map for optimization purposes.
"""
struct Constraint
    node1::Int                  # Index of the first node
    node2::Int                  # Index of the second node
    k::Float64                  # Stiffness coefficient
    η::Float64                 # Damping coefficient
    local_sparsity_map::SVector{144, Int}  # Pre-allocated local sparsity map 
    global_sparsity_map::SVector{144, Int}  # Pre-allocated global sparsity map 
end

#----------------------------------
# FUNCTIONS TO BUILD CONSTRAINTS
#----------------------------------

"""
    build_constraints(nodespairs::AbstractVector, k, η)

Constructs a vector of `Constraint` objects from a vector of node index pairs.

- `nodespairs`: A vector of tuples or arrays, each containing a pair of node indices.
- `k`: Stiffness value(s), can be a scalar or vector.
- `η`: Damping value(s), can be a scalar or vector.
"""
function Constraints(nodespairs::AbstractVector, k, η)
    constraints = StructArray(
        Constraint(nodespairs[i][1], nodespairs[i][2], 
                   k isa AbstractVector ? k[i] : k, 
                   η isa AbstractVector ? η[i] : η) 
        for i in 1:size(nodespairs, 1)
    )
    return constraints
end 

"""
    build_constraints(nodespairs::AbstractMatrix, k, η)

Constructs a vector of `Constraint` objects from a 2-column matrix of node index pairs.

- `nodespairs`: A matrix where each row contains a pair of node indices.
- `k`: Stiffness value(s), can be a scalar or vector.
- `η`: Damping value(s), can be a scalar or vector.
"""
function Constraints(nodespairs::AbstractMatrix, k, η)
    constraints = StructArray(
        Constraint(nodespairs[i, 1], nodespairs[i, 2], 
                   k isa AbstractVector ? k[i] : k, 
                   η isa AbstractVector ? η[i] : η) 
        for i in 1:size(nodespairs, 1)
    )
    return constraints
end 

"""
    Constraint(node1, node2, k, η)

Constructor of Constraint type that initializes `sparsity_map` with zeros.
"""
function Constraint(node1, node2, k, η)
    return Constraint(node1, node2, k, η, zeros(Int, 144), zeros(Int, 144))
end
