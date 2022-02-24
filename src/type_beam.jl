#----------------------------------
# STRUCTURE
#----------------------------------

struct Beam{T, Tᵏ}
    
    ind::Int # index of this beam 
    node1::Int # index of node 1   
    node2::Int # index of node 2
    l₀::T # initial beam length
    Rₑ⁰::Mat33{T} # initial beam rotation matrix
    K̄ⁱⁿᵗ::Tᵏ # beam internal matrix
    sparsity_map::Vector{Int} # sparsity map from local indices (beam matrix) to global indices (gloabl sparse matrix) -> computed in the constructor of the sparse matrices

end

#----------------------------------
# CONSTRUCTORS
#----------------------------------
"""
    beams = constructor_beams(nodes, connectivity, mat, geom, numberInterpolationPoints, Rₑ⁰=nothing, T=Float64)

Constructor of the beams StructArray:
- `nodes`: nodes StructArray (created with constructor_nodes);
- `connectivity`: connectivity of the mesh (Vec2{Int});
- `mat`: struct containing the material properties of the mesh (Material{T});
- `geom`: struct containing the geomtrical properties of the mesh (Geometry{T});
- `numberInterpolationPoints`: pnumber of points used for the interpolation of the beam centreline (Int);
- `Rₑ⁰`: (not mandatory) initial rotation of the beam elements.

Returns a StructArray{Beam}, structure containing the information of the beam elements. 
"""

function constructor_beams(nodes, connectivity::AbstractMatrix, mat, geom, Rₑ⁰=nothing)
    
    beams = StructArray(constructor_beam(i, nodes[connectivity[i, 1]], nodes[connectivity[i, 2]], mat, geom, Rₑ⁰ isa AbstractMatrix ? Rₑ⁰[i,:] : Rₑ⁰) for i in 1:size(connectivity, 1))

    return beams
    
end 

function constructor_beams(nodes, connectivity::AbstractVector, mat, geom, Rₑ⁰=nothing)
    
    beams = StructArray(constructor_beam(i, nodes[connectivity[i][1]], nodes[connectivity[i][2]], mat, geom, Rₑ⁰ isa AbstractVector ? Rₑ⁰[i] : Rₑ⁰) for i in 1:length(connectivity))

    return beams
    
end 

    
# Constructor of the one beam (Beam) given initial rotation
function constructor_beam(ind, node1::Node{T}, node2::Node{T}, mat, geom, Rₑ⁰) where T
    
    i1 = node1.i   
    i2 = node2.i  
    l₀ = norm(node1.X₀ - node2.X₀)   
    K̄ⁱⁿᵗ = K̄ⁱⁿᵗ_beam(mat, geom, l₀)
    
    return Beam{T, typeof(K̄ⁱⁿᵗ)}(ind, i1, i2, l₀, Rₑ⁰, K̄ⁱⁿᵗ, zeros(Int, 144))
    
end 

constructor_beam(ind, node1::Node{T}, node2::Node{T}, mat, geom, Rₑ⁰::Nothing) where T = constructor_beam(ind, node1, node2, mat, geom, local_R⁰(node1.X₀, node2.X₀))

