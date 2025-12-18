using Pkg; Pkg.activate(".")
using DelimitedFiles, ReadVTK, WriteVTK, LinearAlgebra, Test, Revise, EndoBeams, BenchmarkTools, StaticArrays, Dierckx, StructArrays, Parameters

# cd("stent")
include("utils/utils_IO.jl")
include("utils/utils_build_braided_stent_geometry.jl")
include("utils/utils_morping_cl.jl")
include("utils/utils_positioning.jl")

# --------------------------------------------
# Create stent
# --------------------------------------------

free_positions = readdlm("stent_3PB/input/pos_wallstent_3PB.txt")
connectivity = readdlm("stent_3PB/input/conn_wallstent_3PB.txt", Int)
constraints_connectivity = readdlm("stent_3PB/input/constr_wallstent_3PB.txt", Int)

# --------------------------------------------
# Crimping
# --------------------------------------------

output_dir = "stent_3PB/output3D/"
if !isdir(output_dir) mkpath(output_dir) end

num_nodes = size(free_positions, 1)

# Create nodes structure
zeros_6d = zeros(size(free_positions))  
nodes = NodesBeams(
free_positions,
zeros_6d, zeros_6d, zeros_6d,  # zero initial displacements, velocities, accelerations
zeros_6d, zeros_6d, zeros_6d,  # zero initial rotations, angular velocities, accelerations
"xy"                           # working plane for the model (for cylindrical bcs)
)

#-------------------------------
# Beams and constraints
#-------------------------------
E = 206*1e3
ν = 0.33
ρ = 8*1e-9
radius = 0.051                 # Beam radius (mm)
damping = 0                 # Rayleigh damping coefficient

# Create beams structure
beams = Beams(nodes, connectivity, E, ν, ρ, radius, damping)
num_beams = length(beams)

# Create penalty constraints
k_penalty = 1e3
damping_penalty = 0.1
constraints = Constraints(constraints_connectivity, k_penalty, damping_penalty)

#-------------------------------
# Loads and boundary conditions
#-------------------------------

# Loads
loads = nothing

# Boundary conditions container
bcs = nothing

#-------------------------------
# Assembled configuration
#-------------------------------
conf = BeamsConfiguration(nodes, beams, constraints, loads, bcs)

#----------------------------------
# INTERACTIONS
#----------------------------------

kₙ = 100
μ = 0.0
εᵗ = 1 #regularized parameter for friction contact
ηₙ = 0.0
kₜ = 0
ηₜ = 0
u̇ₛ = 0

inter_properties = InteractionProperties(kₙ, μ, εᵗ, ηₙ, kₜ, ηₜ, u̇ₛ)

# Create master and slave surfaces
surface_master = DiscreteSignedDistanceField("stent_3PB/inputSDF/sdf_0.vtk", false, true)
surface_slave = BeamElementSurface(connectivity) 

# Create the interaction instance
inter = RigidInteraction(surface_master, surface_slave, inter_properties)

#-------------------------------
# Simulation parameters
#-------------------------------
# HHT (Hilber–Hughes–Taylor) α-method parameters
α = -0.05                               
β = 0.25 * (1 - α)^2                    
γ = 0.5 * (1 - 2 * α)                   

params = SimulationParams(
α = α, β = β, γ = γ,
initial_timestep = 1e-1,
min_timestep = 1e-12,
max_timestep = 1e-1,
output_timestep = 1,
simulation_end_time = 100.0,
tolerance_residual = 1e-3,
tolerance_displacement = 1e-3,
max_iterations = 5,
output_dir = output_dir,
verbose = true
)

#-------------------------------
# Run simulation and export results
#-------------------------------
run_simulation!(conf, params, inter)

# export_stent_solution_to_txt(
# nodes, beams, num_nodes, num_beams,
# output_dir_crimping,
# true
# )
