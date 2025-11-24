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
    point::Vec3{Float64}        # any point on the plane
    normal::Vec3{Float64}      # unit normal vector
end

# Defines an analytical sphere used as a master surface in contact problems
struct SphereSurface <: RigidBodySurface
    center::Vec3{Float64}  # Sphere center (x, y, z)
    radius::Float64  # Sphere radius
end

# Defines a discrete signed distance field used as a master surface in contact problems
struct DiscreteSignedDistanceField{Tsitp} <: RigidBodySurface
    flag_load_from_file_iterative::Bool
    sitp::Tsitp                  # Scaled interpolation of the SDF field
    dom::NTuple{6, Float64}  # Domain boundary coordinates
    dx::Float64              # Grid spacing in x direction
    dy::Float64              # Grid spacing in y direction
    dz::Float64              # Grid spacing in z direction
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

# Creates a BeamElementSurface from a connectivity matrix.
function BeamElementSurface(beams_connectivity)
    contact_beams = collect(1:size(beams_connectivity, 1)) # Creates a StructArray of beam indices
    return BeamElementSurface(contact_beams)
end

# Creates a DiscreteSignedDistanceField reading the SDF data from a VTK file and setting up interpolation.
function DiscreteSignedDistanceField(filename::String, inside::Bool, flag_load_from_file_iterative=false)

    # Read SDF data from VTK file
    npx, npy, npz, dx, dy, dz, dom, sdf = read_vtk_sdf(filename)
    field = inside ? reshape(sdf, (npx, npy, npz)) : reshape(-sdf, (npx, npy, npz))

    # Define coordinate ranges for interpolation
    x, y, z = range(dom[1]; step=dx, stop=dom[2]), range(dom[3]; step=dy, stop=dom[4]), range(dom[5]; step=dz, stop=dom[6])

    # Create a quadratic interpolation for smooth gradients
    itp = interpolate(field, BSpline(Quadratic(Reflect(OnCell()))))
    sitp = scale(itp, x, y, z)  # Scaled interpolation over the coordinate grid

    return DiscreteSignedDistanceField{typeof(sitp)}(flag_load_from_file_iterative, sitp, dom, dx, dy, dz)  
    
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
