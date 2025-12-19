
#----------------------------------------------------
# PACKAGE INITIALIZATION
#----------------------------------------------------

using Pkg; Pkg.activate(".")
using DelimitedFiles, ReadVTK, WriteVTK, LinearAlgebra, Test, Revise, EndoBeams, BenchmarkTools

#----------------------------------------------------# BEAM DEFINITIONS
#----------------------------------------------------

# Define node positions and connectivity for beams
positions = [0.0  0.0  0.0;
0.0  2.5  0.0;
0.0  5.0  0.0;
0.0  7.5  0.0;
0.0 10.0  0.0;
-2.5 10.0  0.0;
-5.0 10.0  0.0;
-7.5 10.0  0.0;
-10.0 10.0  0.0]             
nnodes = size(positions, 1)    # Total number of nodes

nbeams = nnodes - 1                    # Number of beam elements
connectivity = [1 2;
2 3;
3 4;
4 5;
5 6;
6 7;
7 8;
8 9]  # Connectivity for beam elements # Connectivity for beam elements

# Initial conditions for displacements, velocities, accelerations, and rotations
initial_displacements = zeros(size(positions))       # Zero initial displacements
initial_velocities = zeros(size(positions))       # Zero initial velocities
initial_accelerations = zeros(size(positions))       # Zero initial accelerations
initial_rotations = zeros(size(positions))        # Zero initial rotations
initial_angular_velocities = zeros(size(positions))   # Zero initial angular velocities
initial_angular_accelerations = zeros(size(positions))  # Zero initial angular accelerations
plane = "xy"                         # Plane of the problem

# Material and geometric properties for beams
E = 1e6                              # Young's modulus (Pa)
ν = 0.3                              # Poisson's ratio
ρ = 1                                # Density (kg/m³)
radius = 0.3                         # Beam radius (m)
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

# External force
loaded_dofs = [27]   # Degree of freedom where the external force is applied
force_function(t,i) = 100 * t # Time-dependent external force function
concentrated_force = ConcentratedForce(force_function, loaded_dofs)  

# Degrees of freedom (DOFs) definition
ndofs = nnodes * 6                     # Total number of DOFs (6 per node for displacement and rotation)
blocked_dofs = 1:6                      # First 6 DOFs are fixed (e.g., at the boundary)
encastre = Encastre(blocked_dofs)

# Beam configuration struct initialization
conf = BeamsConfiguration(nodes, beams, Loads(concentrated_force), BoundaryConditions(encastre, ndofs))  # Store the configuration for beams

#----------------------------------------------------
# SOLVER DEFINITIONS
#----------------------------------------------------

# HHT (Houbolt-Hughes-Taylor) time stepping parameters
α = -0.2   # Typically between 0 and -1, used for numerical stability
β = 0.25 * (1 - α)^2  # Damping parameter based on α
γ = 0.5 * (1 - 2 * α)  # Time-stepping parameter

# General time stepping parameters
initial_timestep = 1e-4    # Initial time step size
min_timestep = 1e-4    # Minimum allowed time step
max_timestep = 1e-4    # Maximum allowed time step (could be adjusted based on system behavior)
output_timestep = 1e-2    # Time step for output plotting or visualization
simulation_end_time = 0.15  # End time for the simulation (duration of the analysis)

# Convergence criteria for the solver
tolerance_residual = 1e-3   # Residual tolerance for convergence checks
tolerance_displacement = 1e-3    # Tolerance for changes in displacement (ΔD)
max_iterations = 20      # Maximum number of iterations for the solver

# Store solver parameters in a structured Params object
params = SimulationParams(;
α, β, γ, initial_timestep, min_timestep, max_timestep, output_timestep, simulation_end_time, tolerance_residual, tolerance_displacement, max_iterations, output_dir = "test/output3D"
, verbose = true)

#----------------------------------------------------
# START SIMULATION
#----------------------------------------------------

# Run the solver with the defined configurations and parameters
run_simulation!(conf, params)

if simulation_end_time == 0.15
    @test conf.nodes[end].uⁿ[3] ≈ -0.0100 atol=1e-2
elseif simulation_end_time == 1.5
    @test conf.nodes[end].uⁿ[3] ≈ 5.8332 atol=1e-2
elseif simulation_end_time == 15
    @test conf.nodes[end].uⁿ[3] ≈ 10.3849 atol=1e-2
end 
