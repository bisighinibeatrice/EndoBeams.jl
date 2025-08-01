
#----------------------------------
# BEAM STRUCTURE DEFINITION
#----------------------------------

struct Beam{Tp}
    
    ind::Int # index of this beam 
    node1::Int # index of node 1   
    node2::Int # index of node 2
    l₀::Float64 # initial beam length
    Rₑ⁰::Mat33{Float64} # initial beam rotation matrix
    properties::Tp # beam internal matrix
    local_sparsity_map::SVector{144, Int} # DOF mapping for local matrix assembly                  
    global_sparsity_map::SVector{144, Int} # DOF mapping for global matrix assembly  

end

struct BeamMaterialProperties{TK, TJ}
    radius::Float64          # Beam radius
    E::Float64               # Young's modulus
    K̄ⁱⁿᵗ::TK                 # Internal stiffness matrix
    Jᵨ::TJ                   # Rotational inertia matrix
    Aᵨ::Float64              # Cross-sectional area
    damping::Float64         # Damping coefficient
end

#----------------------------------
# FUNCTIONS TO BUILD BEAMS
#----------------------------------

function Beams(nodes, connectivity::AbstractMatrix, E, ν, ρ, radius, damping, Rₑ⁰=nothing)
    
    beams = StructArray(Beam(i, 
                            nodes[connectivity[i, 1]], 
                            nodes[connectivity[i, 2]], 
                            E isa AbstractVector ? E[i] : E,
                            ν isa AbstractVector ? ν[i] : ν,
                            ρ isa AbstractVector ? ρ[i] : ρ,
                            radius isa AbstractVector ? radius[i] : radius,
                            damping isa AbstractVector ? damping[i] : damping,
                            Rₑ⁰ isa AbstractMatrix ? Rₑ⁰[i,:] : Rₑ⁰) for i in 1:size(connectivity,1))

    return beams
end 

function Beams(nodes, connectivity::AbstractVector, E, ν, ρ, radius, damping, Rₑ⁰=nothing)
    
    beams = StructArray(Beam(i, 
                            nodes[connectivity[i][1]], 
                            nodes[connectivity[i][2]], 
                            E isa AbstractVector ? E[i] : E,
                            ν isa AbstractVector ? ν[i] : ν,
                            ρ isa AbstractVector ? ρ[i] : ρ,
                            radius isa AbstractVector ? radius[i] : radius,
                            damping isa AbstractVector ? damping[i] : damping,
                            Rₑ⁰ isa AbstractVector ? Rₑ⁰[i] : Rₑ⁰) for i in 1:size(connectivity,1))

    return beams  
end 
    
function Beam(ind, node1::NodeBeam, node2::NodeBeam, E, ν, ρ, radius, damping, Rₑ⁰)

    i1 = node1.i   
    i2 = node2.i  
    l₀ = norm(node1.X₀ - node2.X₀)   
    beamprops = compute_material_properties(l₀, E, ν, ρ, radius, damping)
    
    return Beam{typeof(beamprops)}(ind, i1, i2, l₀, Rₑ⁰, beamprops, zeros(Int, 144), zeros(Int, 144)) 
end 

Beam(ind, node1::NodeBeam, node2::NodeBeam, E, ν, ρ, radius, damping, Rₑ⁰::Nothing) = Beam(ind, node1, node2, E, ν, ρ, radius, damping, get_Rₑ⁰(node1.X₀, node2.X₀))

function compute_material_properties(l₀, E, ν, ρ, radius, damping)
    G = E / (2 * (1 + ν))                        # Shear modulus
    A = pi * radius^2                            # Cross-sectional area
    I₂₂ = pi * radius^4 / 4                      # Second moment of inertia (axis 2)
    I₃₃ = I₂₂                                    # Second moment of inertia (axis 3)
    Iₒ = I₂₂ + I₃₃                               # Polar moment of inertia
    Jᵨ = ρ * Diagonal(Vec3(Iₒ, I₂₂, I₃₃))        # Rotational inertia matrix
    Aᵨ = ρ * A                                   # Area mass density
    K̄ⁱⁿᵗ = K̄ⁱⁿᵗ_beam(E, G, Iₒ, A, I₂₂, I₃₃, l₀)  # Beam stiffness matrix
    
    return BeamMaterialProperties{typeof(K̄ⁱⁿᵗ), typeof(Jᵨ)}(radius, E, K̄ⁱⁿᵗ, Jᵨ, Aᵨ, damping)
end

#----------------------------------
# UTILS TO BUILD BEAMS
#----------------------------------

# Function to compute the rotation matrix Rₑ⁰ between two vectors X1 and X2.
# It returns the rotation matrix from X1 to X2, given the direction cosine between them.
function get_Rₑ⁰(X1, X2, Float64=Float64) 
    
    l0 = norm(X2 - X1)  # Length of the vector difference between X2 and X1
    E1 = Vec3(1, 0, 0)  # Unit vector along the x-axis
    E1_0 = (X2 - X1) / l0  # Normalize the difference vector between X1 and X2
    
    v = cross(E1, E1_0)  # Cross product of E1 and the normalized vector
    
    s = norm(v)  # Magnitude of the cross product (sin(theta))
    c = dot(E1, E1_0)  # Dot product between E1 and the normalized vector (cos(theta))
    
    # Handle the special cases where the cosine of the angle is very close to -1 or 1
    if c < -(1 - 2 * eps(Float64))
        # If the cosine is approximately -1, return a 180-degree rotation matrix
        return Mat33{Float64}(-1, 0, 0, 0, 1, 0, 0, 0, -1)
    elseif c > (1 - 2 * eps(Float64))
        # If the cosine is approximately 1, return the identity matrix (no rotation)
        return ID3
    else 
        # Otherwise, compute the Rodrigues' rotation formula for rotation matrix
        Sv = skew(v)  # Skew-symmetric matrix of v
        return ID3 + Sv + ((1 - c) / s^2) * Sv * Sv
    end      
end

# Computes the stiffness matrix components for a beam element.
# This function calculates the axial stiffness, rotational stiffness, and additional rotational stiffness 
# components for a beam element using its material properties and geometric properties.
#
# Inputs:
# - `E`   : Young's Modulus (Pa) of the beam material
# - `G`   : Shear Modulus (Pa) of the beam material
# - `Iₒ`  : Polar moment of inertia (m⁴) of the beam cross-section
# - `A`   : Cross-sectional area (m²) of the beam
# - `I₂₂` : Second moment of area (m⁴) about the second axis (bending stiffness in the y-axis)
# - `I₃₃` : Second moment of area (m⁴) about the third axis (bending stiffness in the z-axis)
# - `l₀`  : Length of the beam element (m)
#
# Returns:
# - `K̄ⁱⁿᵗū`  : Axial stiffness matrix component (a scalar)
# - `K̄ⁱⁿᵗΘ̅`  : Rotational stiffness matrix (a diagonal matrix with rotational stiffness components)
# - `K̄ⁱⁿᵗΘ̅Θ̅`: Additional rotational stiffness matrix (a diagonal matrix with additional rotational stiffness components)
@inline function K̄ⁱⁿᵗ_beam(E, G, Iₒ, A, I₂₂, I₃₃, l₀)
    
    # Compute the axial stiffness matrix component (K̄ⁱⁿᵗū)
    K̄ⁱⁿᵗū = A * E / l₀  # Formula for axial stiffness (A: cross-sectional area, E: Young's Modulus, l₀: beam length)
    
    # Compute the rotational stiffness matrix components (K̄ⁱⁿᵗΘ̅ and K̄ⁱⁿᵗΘ̅Θ̅)
    # K̄ⁱⁿᵗΘ̅ represents the torsion and bending stiffness
    K̄ⁱⁿᵗΘ̅ = Diagonal(@SVector [G * Iₒ / l₀, 4 * E * I₃₃ / l₀, 4 * E * I₂₂ / l₀])
    
    # K̄ⁱⁿᵗΘ̅Θ̅ represents additional rotational stiffness components
    K̄ⁱⁿᵗΘ̅Θ̅ = Diagonal(@SVector [-G * Iₒ / l₀, 2 * E * I₃₃ / l₀, 2 * E * I₂₂ / l₀])
    
    # Return the computed axial and rotational stiffness matrices
    return K̄ⁱⁿᵗū, K̄ⁱⁿᵗΘ̅, K̄ⁱⁿᵗΘ̅Θ̅
end