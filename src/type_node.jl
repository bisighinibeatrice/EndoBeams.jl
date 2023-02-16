#----------------------------------
# STRUCTURE
#----------------------------------
struct Node
   
    i::Int # index
    X₀::Vec3{Float64} # initial position

    # dof: total, displacement, angular
    idof_6::Vec6{Int}
    idof_disp::Vec3{Int}
    idof_rot::Vec3{Int}

    # current configuration of the node (@n+1)
    u::Vec3{Float64} # displacement
    u̇::Vec3{Float64} # velocity
    ü::Vec3{Float64} # acceleration
    w::Vec3{Float64} # angle (spin vecotr)
    ẇ::Vec3{Float64} # angular velocity
    ẅ::Vec3{Float64} # angular acceleration
    R::Mat33{Float64} # local rotation matrix
    ΔR::Mat33{Float64} # local rotation matrix variation

    # last configuration of the node (@n)
    uⁿ::Vec3{Float64}
    u̇ⁿ::Vec3{Float64}
    üⁿ::Vec3{Float64}
    wⁿ::Vec3{Float64}
    ẇⁿ::Vec3{Float64}
    ẅⁿ::Vec3{Float64}
    Rⁿ::Mat33{Float64}
    ΔRⁿ::Mat33{Float64}

    R_global_to_local::Mat33{Float64}

end

#----------------------------------
# CONSTRUCTOR
#----------------------------------

"""
nodes = nodes(X, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, R=nothing) 

Constructor of the nodes StructArray:
- `X`: nodes StructArray (created with nodes);
- `u⁰`: initial displacements;
- `u̇⁰`: initial velocities;
- `ü⁰`: initial accelerations;
- `w⁰`: initial rotations;
- `ẇ⁰`: initial rotation velocities;
- `ẅ⁰`: initial rotation acceleration;
- `R`: (not mandatory) initial rotation of the nodes.

Returns a StructArray{Node}, structure containing the information of the nodes. 
"""
function build_nodes(X::AbstractMatrix, u⁰::AbstractMatrix, u̇⁰::AbstractMatrix, ü⁰::AbstractMatrix, w⁰::AbstractMatrix, ẇ⁰::AbstractMatrix, ẅ⁰::AbstractMatrix, plane=nothing, R=nothing) 

    nodes = StructArray(Node(
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
            R isa AbstractMatrix ? R[i,:] : ID3, 
            ID3,
            Vec3(u⁰[i,:]), 
            Vec3(u̇⁰[i,:]), 
            Vec3(ü⁰[i,:]), 
            Vec3(w⁰[i,:]), 
            Vec3(ẇ⁰[i,:]), 
            Vec3(ẅ⁰[i,:]), 
            R isa AbstractMatrix ? R[i,:] : ID3,
            ID3,
            plane isa String ? compute_local_to_global_matrix(X[i, :], plane) : ID3) for i in 1:size(X, 1))
    return nodes

end


function build_nodes(X::AbstractVector, u⁰::AbstractVector, u̇⁰::AbstractVector, ü⁰::AbstractVector, w⁰::AbstractVector, ẇ⁰::AbstractVector, ẅ⁰::AbstractVector, plane=nothing, R=nothing) 

    nnodes = length(X)
    nodes = StructArray{Node}((
            1:nnodes, 
            convert(Vector{Vec3{Float64}}, X), 
            [Vec6(6*(i-1).+(1,2,3,4,5,6)) for i in 1:nnodes],
            [Vec3(6*(i-1).+(1,2,3)) for i in 1:nnodes],
            [Vec3(6*(i-1).+(4,5,6)) for i in 1:nnodes],
            convert(Vector{Vec3{Float64}}, u⁰),  
            convert(Vector{Vec3{Float64}}, u̇⁰),  
            convert(Vector{Vec3{Float64}}, ü⁰),  
            convert(Vector{Vec3{Float64}}, w⁰),  
            convert(Vector{Vec3{Float64}}, ẇ⁰),  
            convert(Vector{Vec3{Float64}}, ẅ⁰),  
            R isa AbstractVector ? convert(Vector{Mat33{Float64}}, R) : [Mat33{Float64}(ID3) for i in 1:nnodes], 
            [Mat33{Float64}(ID3) for i in 1:nnodes],
            convert(Vector{Vec3{Float64}}, u⁰), 
            convert(Vector{Vec3{Float64}}, u̇⁰), 
            convert(Vector{Vec3{Float64}}, ü⁰), 
            convert(Vector{Vec3{Float64}}, w⁰), 
            convert(Vector{Vec3{Float64}}, ẇ⁰), 
            convert(Vector{Vec3{Float64}}, ẅ⁰), 
            R isa AbstractVector ? convert(Vector{Mat33{Float64}}, R) : [Mat33{Float64}(ID3) for i in 1:nnodes],
            [Mat33{Float64}(ID3) for i in 1:nnodes],
            plane isa String ? [Mat33(compute_local_to_global_matrix(X[i][:], plane)) for i in 1:nnodes] : [Mat33{Float64}(ID3) for i in 1:nnodes]))

    return nodes

end


function compute_local_to_global_matrix(Xi, plane)

    if plane == "xy"
        thetai = atan(Xi[2], Xi[1])
        return Mat33(cos(thetai), -sin(thetai), 0, sin(thetai), cos(thetai), 0, 0, 0, 1)
    elseif plane == "yz"
        thetai = atan(Xi[3], Xi[2])
        return Mat33(1, 0, 0, 0, cos(thetai), -sin(thetai), 0, sin(thetai), cos(thetai))
    elseif plane == "xz"
        thetai = atan(Xi[3], Xi[1])
        return Mat33(cos(thetai), 0, -sin(thetai), 0, 1, 0, sin(thetai), 0, cos(thetai))
    end 
    
end
