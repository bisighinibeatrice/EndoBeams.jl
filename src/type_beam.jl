#----------------------------------
# STRUCTURE
#----------------------------------

struct Beam{Tp}
    
    ind::Int # index of this beam 
    node1::Int # index of node 1   
    node2::Int # index of node 2
    l₀::Float64 # initial beam length
    Rₑ⁰::Mat33{Float64} # initial beam rotation matrix
    properties::Tp # beam internal matrix
    sparsity_map::SVector{144, Int} # sparsity map from local indices (beam matrix) to global indices (gloabl sparse matrix) -> computed in the constructor of the sparse matrices

end

#----------------------------------
# CONSTRUCTORS
#----------------------------------
"""
    beams = beams(nodes, connectivity, beamprops, Rₑ⁰=nothing)

Constructor of the beams StructArray:
- `nodes`: nodes StructArray (created with nodes);
- `connectivity`: connectivity of the mesh (Vec2{Int});
- `material`: struct containing the material properties of the mesh (Material{Float64});
- `geometry`: struct containing the geomtrical properties of the mesh (Geometry{Float64});
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
function Beam(ind, node1::Node, node2::Node, E, ν, ρ, radius, damping, Rₑ⁰)
    
    i1 = node1.i   
    i2 = node2.i  
    l₀ = norm(node1.X₀ - node2.X₀)   
    beamprops = BeamProperties(l₀, E, ν, ρ, radius, damping)
    
    return Beam{typeof(beamprops)}(ind, i1, i2, l₀, Rₑ⁰, beamprops, zeros(Int, 144))
    
end 

Beam(ind, node1::Node, node2::Node, E, ν, ρ, radius, damping, Rₑ⁰::Nothing) = Beam(ind, node1, node2, E, ν, ρ, radius, damping, get_Rₑ⁰(node1.X₀, node2.X₀))

function get_Rₑ⁰(X1, X2, T=Float64) 
    
    l0 = norm(X2-X1)
    E1 = Vec3(1, 0, 0)
    E1_0 = (X2-X1)/l0
    v = cross(E1, E1_0)
    
    s = norm(v)
    c = dot(E1, E1_0)
    
    if (c<-(1- 2*eps(T)))
        return Mat33{T}(-1, 0, 0, 0, 1, 0, 0, 0, -1)
    elseif (c > (1- 2*eps(T)))
       return  ID3
    else 
        Sv = skew(v)
        return ID3 + Sv + ((1-c)/s^2)*Sv*Sv
    end 
        
end 