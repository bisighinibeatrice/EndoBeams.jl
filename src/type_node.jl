#----------------------------------
# STRUCTURE
#----------------------------------
struct Node{T}
   
    i::Int # index
    X₀::Vec3{T} # initial position

    # dof: total, displacement, angular
    idof_6::Vec6{Int}
    idof_disp::Vec3{Int}
    idof_rot::Vec3{Int}

    # current configuration of the node (@n+1)
    u::Vec3{T} # displacement
    u̇::Vec3{T} # velocity
    ü::Vec3{T} # acceleration
    w::Vec3{T} # angle (spin vecotr)
    ẇ::Vec3{T} # angular velocity
    ẅ::Vec3{T} # angular acceleration
    R::Mat33{T} # local rotation matrix
    ΔR::Mat33{T} # local rotation matrix variation

    # last configuration of the node (@n)
    uⁿ::Vec3{T}
    u̇ⁿ::Vec3{T}
    üⁿ::Vec3{T}
    wⁿ::Vec3{T}
    ẇⁿ::Vec3{T}
    ẅⁿ::Vec3{T}
    Rⁿ::Mat33{T}
    ΔRⁿ::Mat33{T}

end

#----------------------------------
# CONSTRUCTOR
#----------------------------------

"""
nodes = constructor_nodes(X, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, Rₑ⁰=nothing, T=Float64) 

Constructor of the nodes StructArray:
- `X`: nodes StructArray (created with constructor_nodes);
- `u⁰`: initial displacements;
- `u̇⁰`: initial velocities;
- `ü⁰`: initial accelerations;
- `w⁰`: initial rotations;
- `ẇ⁰`: initial rotation velocities;
- `ẅ⁰`: initial rotation acceleration;
- `Rₑ⁰`: (not mandatory) initial rotation of the nodes.

Returns a StructArray{Node}, structure containing the information of the nodes. 
"""
function constructor_nodes(X::AbstractMatrix, u⁰::AbstractMatrix, u̇⁰::AbstractMatrix, ü⁰::AbstractMatrix, w⁰::AbstractMatrix, ẇ⁰::AbstractMatrix, ẅ⁰::AbstractMatrix, Rₑ⁰=nothing, T=Float64) 

    nodes = StructArray(Node{T}(
            i, 
            X[i, :], 
            Vec6(6*(i-1).+(1,2,3,4,5,6)),
            Vec3(6*(i-1).+(1,2,3)), 
            Vec3(6*(i-1).+(4,5,6)), 
            Vec3(u⁰[i,:]), 
            Vec3(u̇⁰[i,:]), 
            Vec3(ü⁰[i,:]), 
            Vec3(w⁰[i,:]), 
            Vec3(ẇ⁰[i,:]), 
            Vec3(ẅ⁰[i,:]), 
            Rₑ⁰ isa AbstractMatrix ? Rₑ⁰[i,:] : ID3, 
            ID3,
            Vec3(u⁰[i,:]), 
            Vec3(u̇⁰[i,:]), 
            Vec3(ü⁰[i,:]), 
            Vec3(w⁰[i,:]), 
            Vec3(ẇ⁰[i,:]), 
            Vec3(ẅ⁰[i,:]), 
            Rₑ⁰ isa AbstractMatrix ? Rₑ⁰[i,:] : ID3,
            ID3) for i in 1:size(X, 1))

    return nodes

end


function constructor_nodes(X::AbstractVector, u⁰::AbstractVector, u̇⁰::AbstractVector, ü⁰::AbstractVector, w⁰::AbstractVector, ẇ⁰::AbstractVector, ẅ⁰::AbstractVector, Rₑ⁰=nothing, T=Float64) 

    nnodes = length(X)
    nodes = StructArray{Node{T}}((
            1:nnodes, 
            convert(Vector{Vec3{T}}, X), 
            [Vec6(6*(i-1).+(1,2,3,4,5,6)) for i in 1:nnodes],
            [Vec3(6*(i-1).+(1,2,3)) for i in 1:nnodes],
            [Vec3(6*(i-1).+(4,5,6)) for i in 1:nnodes],
            convert(Vector{Vec3{T}}, u⁰),  
            convert(Vector{Vec3{T}}, u̇⁰),  
            convert(Vector{Vec3{T}}, ü⁰),  
            convert(Vector{Vec3{T}}, w⁰),  
            convert(Vector{Vec3{T}}, ẇ⁰),  
            convert(Vector{Vec3{T}}, ẅ⁰),  
            Rₑ⁰ isa AbstractVector ? convert(Vector{Mat33{T}}, Rₑ⁰) : [Mat33{T}(ID3) for i in 1:nnodes], 
            [Mat33{T}(ID3) for i in 1:nnodes],
            convert(Vector{Vec3{T}}, u⁰), 
            convert(Vector{Vec3{T}}, u̇⁰), 
            convert(Vector{Vec3{T}}, ü⁰), 
            convert(Vector{Vec3{T}}, w⁰), 
            convert(Vector{Vec3{T}}, ẇ⁰), 
            convert(Vector{Vec3{T}}, ẅ⁰), 
            Rₑ⁰ isa AbstractVector ? convert(Vector{Mat33{T}}, Rₑ⁰) : [Mat33{T}(ID3) for i in 1:nnodes],
            [Mat33{T}(ID3) for i in 1:nnodes]))

    return nodes

end


