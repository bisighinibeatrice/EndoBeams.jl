#----------------------------------------------------
# PACKAGE INITIALIZATION
#----------------------------------------------------

using Pkg
using Pkg; Pkg.activate(".")
using DelimitedFiles, ReadVTK, WriteVTK, LinearAlgebra, Test, Revise, EndoBeams, BenchmarkTools

#----------------------------------------------------
# BEAM MESH DEFINITIONS
#----------------------------------------------------

# Define node positions and connectivity for beams
positions =  readdlm("test/input/pos_net.txt")        
nnodes = size(positions, 1)    # Total number of nodes

nbeams = nnodes - 1                    # Number of beam elements
connectivity = readdlm("test/input/conn_net.txt", Int)   # Connectivity for beam elements

# Initial conditions for displacements, velocities, accelerations, and rotations
initial_displacements = zeros(size(positions))       # Zero initial displacements
initial_velocities = zeros(size(positions))       # Zero initial velocities
initial_velocities[:, 2] .= -0.001
initial_accelerations = zeros(size(positions))       # Zero initial accelerations
initial_rotations = zeros(size(positions))        # Zero initial rotations
initial_angular_velocities = zeros(size(positions))   # Zero initial angular velocities
initial_angular_accelerations = zeros(size(positions))  # Zero initial angular accelerations
plane = "xy"                         # Plane of the problem

# Material and geometric properties for beams
E = 10                              # Young's modulus (Pa)
ν = 0                             # Poisson's ratio
ρ = 1500                                # Density (kg/m³)
radius = 1e-4*2                       # Beam radius (m)
damping = 0                          # Damping coefficient

# Build nodes and beams
nodes = NodesBeams(
positions, 
initial_displacements, 
initial_velocities, 
initial_accelerations, 
initial_rotations, 
initial_angular_velocities, 
initial_angular_accelerations, 
plane
)
beams = Beams(nodes, connectivity, E, ν, ρ, radius, damping)

#----------------------------------
# BEAMS CONFIGURATION DEFINITIONS
#----------------------------------

# Initialize the loads and boundary conditions as nothing
loads = nothing
bcs = nothing

# Mesh configuration initialization
conf = BeamsConfiguration(nodes, beams, loads, bcs)

#----------------------------------
# INTERACTIONS
#----------------------------------

# Create interaction properties (using defaults for some values)
kₙ = 0.01 # Penalty parameter
μ = 0.01 # Friction coefficient
εᵗ = 0.5 # Regularized parameter for friction contact
ηₙ = 0
kₜ = 0
ηₜ = 0
u̇ₛ = 0  
inter_properties = InteractionProperties(kₙ, μ, εᵗ, ηₙ, kₜ, ηₜ, u̇ₛ)

# Create master and slave surfaces
surface_master = DiscreteSignedDistanceField("test/input/sdf_sphere_80.vtk", true, false)
surface_slave = BeamElementSurface(connectivity) 

# Create the interaction instance
inter = RigidInteraction(surface_master, surface_slave, inter_properties)

#----------------------------------------------------
# SOLVER DEFINITIONS
#----------------------------------------------------

# HHT (Houbolt-Hughes-Taylor) time stepping parameters
α = -0.05   # Typically between 0 and -1, used for numerical stability
β = 0.25 * (1 - α)^2  # Damping parameter based on α
γ = 0.5 * (1 - 2 * α)  # Time-stepping parameter

# General time stepping parameters
initial_timestep = 0.5    # Initial time step size
min_timestep = 0.5     # Minimum allowed time step
max_timestep = 0.5    # Maximum allowed time step (could be adjusted based on system behavior)
output_timestep = 0.5    # Time step for output plotting or visualization
simulation_end_time = 100 # End time for the simulation (duration of the analysis)

# Convergence criteria for the solver
tolerance_residual = 1e-5  # Residual tolerance for convergence checks
tolerance_displacement = 1e-5    # Tolerance for changes in displacement (ΔD)
max_iterations = 10      # Maximum number of iterations for the solver

# Store solver parameters in a structured Params object
params = SimulationParams(;
α, β, γ, initial_timestep, min_timestep, max_timestep, output_timestep, simulation_end_time, tolerance_residual, tolerance_displacement, max_iterations, output_dir = "test/output3D", verbose = true)

#----------------------------------------------------
# START SIMULATION
#----------------------------------------------------

# Run the solver with the defined configurations and parameters
run_simulation!(conf, params, inter)

