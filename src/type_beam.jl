#----------------------------------
# STRUCTURE
#----------------------------------

struct Beam{T, Tp}
    
    ind::Int # index of this beam 
    node1::Int # index of node 1   
    node2::Int # index of node 2
    l₀::T # initial beam length
    Rₑ⁰::Mat33{T} # initial beam rotation matrix
    properties::Tp # beam internal matrix
    sparsity_map::SVector{144, Int} # sparsity map from local indices (beam matrix) to global indices (gloabl sparse matrix) -> computed in the constructor of the sparse matrices

end

#----------------------------------
# CONSTRUCTORS
#----------------------------------
"""
    beams = beams(nodes, connectivity, beamprops, Rₑ⁰=nothing, T=Float64)

Constructor of the beams StructArray:
- `nodes`: nodes StructArray (created with nodes);
- `connectivity`: connectivity of the mesh (Vec2{Int});
- `material`: struct containing the material properties of the mesh (Material{T});
- `geometry`: struct containing the geomtrical properties of the mesh (Geometry{T});
- `Rₑ⁰`: (not mandatory) initial rotation of the beam elements.

Returns a StructArray{Beam}, structure containing the information of the beam elements. 
"""

function build_beams(nodes, connectivity::AbstractMatrix, E, ν, ρ, radius, damping, Rₑ⁰=nothing)
    
    beams = StructArray(Beam(i, 
                            nodes[connectivity[i, 1]], 
                            nodes[connectivity[i, 2]], 
                            E isa AbstractVector ? E[i] : E,
                            ν isa AbstractVector ? ν[i] : ν,
                            ρ isa AbstractVector ? ρ[i] : ρ,
                            radius isa AbstractVector ? radius[i] : radius,
                            damping isa AbstractVector ? damping[i] : damping,
                            Rₑ⁰ isa AbstractMatrix ? Rₑ⁰[i,:] : Rₑ⁰) for i in 1:size(connectivity, 1))

    return beams
    
end 

function build_beams(nodes, connectivity::AbstractVector, E, ν, ρ, radius, damping, Rₑ⁰=nothing)
    
    beams = StructArray(Beam(i, 
                            nodes[connectivity[i][1]], 
                            nodes[connectivity[i][2]], 
                            E isa AbstractVector ? E[i] : E,
                            ν isa AbstractVector ? ν[i] : ν,
                            ρ isa AbstractVector ? ρ[i] : ρ,
                            radius isa AbstractVector ? radius[i] : radius,
                            damping isa AbstractVector ? damping[i] : damping,
                            Rₑ⁰ isa AbstractVector ? Rₑ⁰[i] : Rₑ⁰) for i in 1:length(connectivity))

    return beams
    
end 


    
# Constructor of the one beam (Beam) given initial rotation
function Beam(ind, node1::Node{T}, node2::Node{T}, E, ν, ρ, radius, damping, Rₑ⁰) where T
    
    i1 = node1.i   
    i2 = node2.i  
    l₀ = norm(node1.X₀ - node2.X₀)   
    beamprops = BeamProperties(l₀, E, ν, ρ, radius, damping, T)
    
    return Beam{T, typeof(beamprops)}(ind, i1, i2, l₀, Rₑ⁰, beamprops, zeros(Int, 144))
    
end 

Beam(ind, node1::Node{T}, node2::Node{T}, E, ν, ρ, radius, damping, Rₑ⁰::Nothing) where T = Beam(ind, node1, node2, E, ν, ρ, radius, damping, local_R⁰(node1.X₀, node2.X₀))

