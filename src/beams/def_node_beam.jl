#------------------------------------
# STRUCTURE DEFINITION OF BEAM NODES
#------------------------------------

# Defines the properties and state variables for a beam node in the system
struct NodeBeam
   
    i::Int                     # Index of the node
    X₀::Vec3{Float64}          # Initial position (3D vector)

    # Local dofs (DoFs)
    local_dofs::Vec6{Int}            # Local indices for all 6 DoFs (3 displacement + 3 rotation)
    local_dofs_disp::Vec3{Int}       # Local indices for displacement DoFs
    local_dofs_rot::Vec3{Int}        # Local indices for rotational DoFs
    
    # Global dofs (DoFs)
    global_dofs::Vec6{Int}            # Global indices for all 6 DoFs (3 displacement + 3 rotation)
    global_dofs_disp::Vec3{Int}       # Global indices for displacement DoFs
    global_dofs_rot::Vec3{Int}        # Global indices for rotational DoFs

    # Current configuration of the node (@n+1)
    u::Vec3{Float64}           # Displacement (3D vector)
    u̇::Vec3{Float64}          # Velocity (3D vector)
    ü::Vec3{Float64}          # Acceleration (3D vector)
    w::Vec3{Float64}           # Angular displacement (spin vector)
    ẇ::Vec3{Float64}          # Angular velocity (3D vector)
    ẅ::Vec3{Float64}          # Angular acceleration (3D vector)
    R::Mat33{Float64}          # Local rotation matrix (3x3)
    ΔR::Mat33{Float64}         # Variation in rotation matrix (3x3)

    # Last configuration of the node (@n)
    uⁿ::Vec3{Float64}          # Displacement (3D vector)
    u̇ⁿ::Vec3{Float64}         # Velocity (3D vector)
    üⁿ::Vec3{Float64}         # Acceleration (3D vector)
    wⁿ::Vec3{Float64}          # Angular displacement (spin vector)
    ẇⁿ::Vec3{Float64}         # Angular velocity (3D vector)
    ẅⁿ::Vec3{Float64}         # Angular acceleration (3D vector)
    Rⁿ::Mat33{Float64}         # Local rotation matrix (3x3)
    ΔRⁿ::Mat33{Float64}        # Variation in rotation matrix (3x3)

    # Transformation matrices
    R_carthesian_to_cylindrical::Mat33{Float64} # Transformation matrix from R_carthesian to cylindrical coordinates
        
    # Variables to check contact
    incontact::Int        # Variable to check if in contact
    contactforce::Vec3{Float64}         # Contact force
    contactdistance::Float64
end

#----------------------------------
# FUNCTIONS TO BUILD BEAM NODES
#----------------------------------

# Function to construct a collection of beam nodes from given initial configurations.
function NodesBeams(X::AbstractMatrix, u⁰::AbstractMatrix, u̇⁰::AbstractMatrix, ü⁰::AbstractMatrix, w⁰::AbstractMatrix, ẇ⁰::AbstractMatrix, ẅ⁰::AbstractMatrix, plane=nothing, R=nothing) 

    nnodes = size(X,1)

    nodes = StructArray(NodeBeam(
            i, 
            X[i, :],                          # Initial position
            Vec6(6*(i-1).+(1,2,3,4,5,6)),     # Local indices for 6 DoFs
            Vec3(6*(i-1).+(1,2,3)),           # Local indices for displacement DoFs
            Vec3(6*(i-1).+(4,5,6)),           # Local indices for rotational DoFs
            Vec6(6*(i-1).+(1,2,3,4,5,6)),     # Global indices for 6 DoFs
            Vec3(6*(i-1).+(1,2,3)),           # Global indices for displacement DoFs
            Vec3(6*(i-1).+(4,5,6)),           # Global indices for rotational DoFs
            Vec3(u⁰[i,:]),                    # Initial displacement
            Vec3(u̇⁰[i,:]),                    # Initial velocity
            Vec3(ü⁰[i,:]),                    # Initial acceleration
            Vec3(w⁰[i,:]),                    # Initial angular displacement
            Vec3(ẇ⁰[i,:]),                    # Initial angular velocity
            Vec3(ẅ⁰[i,:]),                    # Initial angular acceleration
            R isa AbstractMatrix ? R[i,:] : ID3,   # Initial rotation matrix (default: identity)
            ID3,                                   # Initial rotation matrix variation (default: identity)
            Vec3(u⁰[i,:]),                         # Initial displacement (@n)
            Vec3(u̇⁰[i,:]),                         # Initial velocity (@n)
            Vec3(ü⁰[i,:]),                         # Initial acceleration (@n)
            Vec3(w⁰[i,:]),                         # Initial angular displacement (@n)
            Vec3(ẇ⁰[i,:]),                         # Initial angular velocity (@n)
            Vec3(ẅ⁰[i,:]),                         # Initial angular acceleration (@n)
            R isa AbstractMatrix ? R[i,:] : ID3,   # Last rotation matrix (@n, default: identity)
            ID3,                                   # Last rotation matrix variation (@n, default: identity)
            plane isa String ? compute_cartesian_to_cylindrical_matrix(X[i, :], plane) : ID3, # Carthesian-to-cylindrical transformation matrix
            0, 
            Vec3(0,0,0), 
            0.0,
            ) for i in 1:nnodes)
    return nodes
end

# Similar function for vector input (alternative format for initial configurations).
function NodesBeams(X::AbstractVector, u⁰::AbstractVector, u̇⁰::AbstractVector, ü⁰::AbstractVector, w⁰::AbstractVector, ẇ⁰::AbstractVector, ẅ⁰::AbstractVector, plane=nothing, R=nothing) 

    nnodes = length(X)
    nodes = StructArray{NodeBeam}(( 
            1:nnodes,                         
            convert(Vector{Vec3{Float64}}, X), 
            [Vec6(6*(i-1).+(1,2,3,4,5,6)) for i in 1:nnodes], 
            [Vec3(6*(i-1).+(1,2,3)) for i in 1:nnodes], 
            [Vec3(6*(i-1).+(4,5,6)) for i in 1:nnodes], 
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
            plane isa String ? [Mat33(compute_cartesian_to_cylindrical_matrix(X[i][:], plane)) for i in 1:nnodes] : [Mat33{Float64}(ID3) for i in 1:nnodes],
            [0 for i in 1:nnodes], 
            [Vec3(0,0,0) for i in 1:nnodes],
            [0.0 for i in 1:nnodes],
            )) 

    return nodes
end

#----------------------------------------
# CARTESIAN-TO-CYLIDRICAL TRANSFORMATION
#----------------------------------------

# Computes the transformation matrix from Cartesian to cylindrical coordinates for a given node.
function compute_cartesian_to_cylindrical_matrix(Xi, plane)

    if plane == "xy"
        thetai = atan(Xi[2], Xi[1])           # Angle in the xy-plane (polar angle)
        return Mat33(cos(thetai), -sin(thetai), 0,
                     sin(thetai),  cos(thetai), 0,
                     0,            0,          1)
    elseif plane == "yz"
        thetai = atan(Xi[3], Xi[2])           # Angle in the yz-plane
        return Mat33(1, 0,            0,
                     0, cos(thetai), -sin(thetai),
                     0, sin(thetai),  cos(thetai))
    elseif plane == "xz"
        thetai = atan(Xi[3], Xi[1])           # Angle in the xz-plane
        return Mat33(cos(thetai), 0, -sin(thetai),
                     0,           1,  0,
                     sin(thetai), 0,  cos(thetai))
    end  
end
