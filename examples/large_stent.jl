using EndoBeams
using DelimitedFiles

const T = Float64



# -------------------------------------------------------------------------------------------
# Read stent information
# -------------------------------------------------------------------------------------------

# read expanded configuration
initial_positions =  readdlm("examples/input_large_stent/positions.txt")
connectivity = readdlm("examples/input_large_stent/connectivity.txt", Int)

# read crimped configuration
final_positions = initial_positions + readdlm("examples/input_large_stent/u.txt")

# -------------------------------------------------------------------------------------------
# Building the nodes
# -------------------------------------------------------------------------------------------

# total number of nodes
nnodes = size(final_positions, 1)

# initial conditions
u⁰ = final_positions - initial_positions

u̇⁰ = zeros(nnodes, 3)
ü⁰ = zeros(nnodes, 3)
w⁰ = zeros(nnodes, 3)
ẇ⁰ = zeros(nnodes, 3)
ẅ⁰ = zeros(nnodes, 3)

R = readdlm("examples/input_large_stent/R.txt")

# nodes StructArray
nodes = build_nodes(initial_positions, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, R, T)



# -------------------------------------------------------------------------------------------
# Building the beams
# -------------------------------------------------------------------------------------------


# geometric and material properties
E = 225*1e3
ν = 0.33
ρ = 9.13*1e-6
radius = 0.065
damping = 1e3

# read initial rotations for the beams
Re₀ = readdlm("examples/input_large_stent/Re0.txt")

# beams vector
beams = build_beams(nodes, connectivity, E, ν, ρ, radius, damping, Re₀)



# contact parameters
kₙ = 10 #penalty parameter
μ = 0.3
εᵗ = 0.1 #regularized parameter for friction contact
ηₙ = 0.1

contact = ContactParameters(kₙ, μ, εᵗ, ηₙ, T)

# -------------------------------------------------------------------------------------------
# External forces
# -------------------------------------------------------------------------------------------

# external force and applied dof
loaded_dofs = T[]
force(t, node_idx) = 0

ext_forces = ExternalForces(force, loaded_dofs)

# -------------------------------------------------------------------------------------------
# Boundary conditions
# -------------------------------------------------------------------------------------------

# number of dof (6 per node)
ndofs = nnodes*6

# penalty constraints
kᶜᵒⁿ = 1e3
ηᶜᵒⁿ = 1
nodespairs = readdlm("examples/input_large_stent/constraints.txt")
constraints = build_constraints(nodespairs, kᶜᵒⁿ, ηᶜᵒⁿ)


# Dirichlet boundary conditions: blocked positions
fixed_dofs = T[]
free_dofs = setdiff(1:ndofs, fixed_dofs)

# Dirichlet dof (x6)
disp_dofs = Int[]
disp_vals = T[]
disp(t, node_idx) = 0


# boundary conditions strucutre
bcs = BoundaryConditions(fixed_dofs, free_dofs, disp, disp_vals, disp_dofs, T)

# -------------------------------------------------------------------------------------------
# SDF
# -------------------------------------------------------------------------------------------

sdf = Discrete_SDF("examples/input_stent/arcStretchObj.vtk", radius, false)

# -------------------------------------------------------------------------------------------
# Final configuration
# -------------------------------------------------------------------------------------------

# configuration: mesh, external forces and boundary conditions
conf = Configuration(nodes, beams, constraints, ext_forces, bcs, contact, sdf)

# -------------------------------------------------------------------------------------------
# Time stepping parameters
# -------------------------------------------------------------------------------------------

# initial time step and total time
ini_Δt = 1e-5
max_Δt = 1.
Δt_plot =  1e-5
tᵉⁿᵈ = 1

params = Params{T}(;ini_Δt, Δt_plot, max_Δt, tᵉⁿᵈ, output_dir = "examples/output3D", stop_on_energy_threshold=true, energy_threshold=1e-10)


# -------------------------------------------------------------------------------------------
# Start simulation
# -------------------------------------------------------------------------------------------


solver!(conf, params, T)