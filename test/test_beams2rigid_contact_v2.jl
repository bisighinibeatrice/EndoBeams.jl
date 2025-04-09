#----------------------------------------------------
# PACKAGE INITIALIZATION
#----------------------------------------------------

using Pkg
using Pkg; Pkg.activate(".")
using DelimitedFiles, ReadVTK, WriteVTK, LinearAlgebra, Test, Revise, EndoBeams, BenchmarkTools

#----------------------------------------------------# BEAM DEFINITIONS
#----------------------------------------------------

# Define node positions and connectivity for beams
L = 10 / 100  # Total length of the beam
dx = L / 20   # Step size for positions along the z-axis
positions = hcat(zeros(Int(L / dx) + 1), zeros(Int(L / dx) + 1), collect(0:dx:L))     
nnodes = size(positions, 1)    # Total number of nodes

# Connectivity for beam elements
nbeams = nnodes - 1                  # Number of beam elements
connectivity = [(i, i+1) for i in 1:nnodes-1] 

# Initial conditions for displacements, velocities, accelerations, and rotations
initial_displacements = zeros(size(positions))       # Zero initial displacements
initial_velocities = zeros(size(positions))       # Zero initial velocities
initial_accelerations = zeros(size(positions))       # Zero initial accelerations
initial_rotations = zeros(size(positions))        # Zero initial rotations
initial_angular_velocities = zeros(size(positions))   # Zero initial angular velocities
initial_angular_accelerations = zeros(size(positions))  # Zero initial angular accelerations
plane = "xy"                         # Plane of the problem

# Material and geometric properties for beams
E = 5*1e7                              # Young's modulus (Pa)
ν = 0.33                              # Poisson's ratio
ρ = 7850                                # Density (kg/m³)
radius = 3e-4                         # Beam radius (m)
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

# Initialize the loads as nothing
loads = nothing

# Encastre
ndofs = nnodes * 6                     # Total number of DOFs (6 per node for displacement and rotation)
blocked_dofs = 2:6                       # First 6 DOFs are fixed (e.g., at the boundary)
encastre = Encastre(blocked_dofs)

# Imposed diplacement
displaced_dof = 1                     # DOFs that will have the imposed displacement
max_displacement = 10 / 100                     # Maximum displacement
time_threshold = 2.5                            # Time threshold for displacement change
velocity = max_displacement / time_threshold    # Velocity coefficient

displacement_function(t) =
(t <= time_threshold) * velocity * t + 
(t > time_threshold && t <= 2 * time_threshold) * (-velocity * t + 2 * max_displacement) + 
(t > 2 * time_threshold && t <= 3 * time_threshold) * velocity * (t - 2 * time_threshold) + 
(t > 3 * time_threshold && t <= 4 * time_threshold) * (-velocity * t + 4 * max_displacement)

imposed_displacement = ImposedDisplacement(displaced_dof, displacement_function)

# Beam configuration struct initialization
bcs =  BoundaryConditions(encastre, imposed_displacement, ndofs)
conf = BeamsConfiguration(nodes, beams, loads, bcs)  # Store the configuration for beams

#----------------------------------
# INTERACTIONS
#----------------------------------

# Create interaction properties (using defaults for some values)
kₙ = 10 # Penalty parameter
μ = 0.0 # Friction coefficient
εᵗ = 100 # Regularized parameter for friction contact
ηₙ = 1e-3 
kₜ = kₙ
ηₜ = ηₙ
u̇ₛ = 0  
inter_properties = InteractionProperties(kₙ, μ, εᵗ, ηₙ, kₜ, ηₜ, u̇ₛ)

# Create master and slave surfaces
center_sphere = [.05,.005,.08]
radius_sphere = 0.04
surface_master = SphereSurface(center_sphere, radius_sphere)
surface_slave = BeamElementSurface(1:nbeams) 

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
initial_timestep = 1e-2    # Initial time step size
min_timestep = 1e-6    # Minimum allowed time step
max_timestep = 1e-2    # Maximum allowed time step (could be adjusted based on system behavior)
output_timestep = 1e-2    # Time step for output plotting or visualization
simulation_end_time = 2.5 # End time for the simulation (duration of the analysis)

# Convergence criteria for the solver
tolerance_residual = 1e-5   # Residual tolerance for convergence checks
tolerance_displacement = 1e-5    # Tolerance for changes in displacement (ΔD)
max_iterations = 5      # Maximum number of iterations for the solver

# Store solver parameters in a structured Params object
params = SimulationParams(;
α, β, γ, initial_timestep, min_timestep, max_timestep, output_timestep, simulation_end_time, tolerance_residual, tolerance_displacement, max_iterations, output_dir = "test/output3D"
, verbose = true)

#----------------------------------------------------
# START SIMULATION
#----------------------------------------------------

# Run the solver with the defined configurations and parameters
run_simulation!(conf, params, inter)
