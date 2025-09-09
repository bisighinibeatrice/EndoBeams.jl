#----------------------------------
# INTERACTION PROPERTIES DEFINITION
#----------------------------------

# Defines the properties for contact interactions (e.g., friction, penalty method).
struct InteractionProperties
    kₙ::Float64   # Normal penalty parameter (stiffness)
    μ::Float64   # Friction coefficient
    εᵗ::Float64   # Regularized parameter for friction contact
    ηₙ::Float64   # Damping parameter in the normal direction
    kₜ::Float64   # Tangential penalty parameter
    ηₜ::Float64   # Damping parameter in the tangential direction
    u̇ₛ::Float64   # Slip velocity threshold for friction regularization
end

#----------------------------------
# SURFACE DEFINITIONS
#----------------------------------

abstract type RigidBodySurface end
abstract type DeformableSurface end

# Defines an analytical plane used as a master surface in contact problems
struct PlaneSurface <: RigidBodySurface
    position::Float64  # Position along a specified axis
    axis::Char         # Axis identifier (:x, :y, :z)
end

# Defines an analytical sphere used as a master surface in contact problems
struct SphereSurface <: RigidBodySurface
    center::Vec3{Float64} # Sphere center (x, y, z)
    radius::Float64 # Sphere radius
end

# Defines a Signed Distance Function (SDF) surface
struct SDFSurface{T} <: RigidBodySurface
    interpolant::T
    domain::Vector{Float64}
    dx::Float64
    dy::Float64
    dz::Float64
end

# Defines a triangulated surface used as a master surface in contact problems
struct TriangulatedSurface <: RigidBodySurface
    positions::Vector{Vec3{Float64}}
    triangles::Vector{Vec3{Int}}
end

# Defines a beam surface composed of beam segments (for contact with rigid bodies)
struct BeamElementSurface <: DeformableSurface
    contact_beams::Vector{Int}
end

#----------------------------------
# SURFACE CREATION FUNCTIONS
#----------------------------------

# Create SDFSurface from a VTK SDF file
function SDFSurface(filename_sdf::String)

    npx, npy, npz, dx, dy, dz, domain, sdf_values = read_vtk_sdf(filename_sdf)
    field = reshape(sdf_values, (npx, npy, npz))

    x = range(domain[1]; step=dx, stop=domain[2])
    y = range(domain[3]; step=dy, stop=domain[4])
    z = range(domain[5]; step=dz, stop=domain[6])

    # Quadratic B-spline interpolation for smooth gradient
    itp = interpolate(field, BSpline(Quadratic(Reflect(OnCell()))))
    interpolant = scale(itp, x, y, z)

    return SDFSurface{typeof(interpolant)}(interpolant, domain, dx, dy, dz)
end
# Creates a BeamElementSurface from a connectivity matrix.
function BeamElementSurface(beams_connectivity::Array{Int, 2})
    contact_beams = collect(1:size(beams_connectivity, 1)) # Creates a StructArray of beam indices
    return BeamElementSurface(contact_beams)
end

#------------------------------------------------
# BEAM-TRIANGLE CONNECTION DEFINITION 
#------------------------------------------------

# Represents a beam-triangle interaction connection.
struct BeamTriangleConnection
    beam_node::Int  # Beam node index
    triangle_nodes::Tuple{Int, Int, Int}  # Triangle node indices
    global_sparsity_map::SVector{144, Int}
end

#----------------------------------
# INTERACTION STRUCTURE DEFINITIONS
#----------------------------------

abstract type Interaction end

# Defines a rigid interaction between a master surface and a slave surface.
struct RigidInteraction <: Interaction
    master::RigidBodySurface
    slave::BeamElementSurface
    properties::InteractionProperties
end
