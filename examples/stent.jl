using EndoBeams
using DelimitedFiles

const T = Float64

# -------------------------------------------------------------------------------------------
# Read stent information
# -------------------------------------------------------------------------------------------

# read expanded configuration
intial_positions =  readdlm("examples/input_stent/pos_stent.txt")
connectivity = readdlm("examples/input_stent/conn_stent.txt", Int)

# read crimped configuration
crimped_positions = intial_positions + readdlm("examples/input_stent/u_crimping.txt")

# read positioned configuration
final_positions =  crimped_positions + readdlm("examples/input_stent/u_positioning.txt")

# -------------------------------------------------------------------------------------------
# Building the nodes
# -------------------------------------------------------------------------------------------

# total number of nodes
nnodes = size(final_positions, 1)

# initial conditions
u⁰ = final_positions - intial_positions

u̇⁰ = zeros(nnodes, 3)
ü⁰ = zeros(nnodes, 3)
w⁰ = zeros(nnodes, 3)
ẇ⁰ = zeros(nnodes, 3)
ẅ⁰ = zeros(nnodes, 3)
R₀ = readdlm("examples/input_stent/R_positioning.txt")


# nodes StructArray
nodes = constructor_nodes(intial_positions, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, R₀, T)

# -------------------------------------------------------------------------------------------
# Building the beams
# -------------------------------------------------------------------------------------------


# geometric and material properties
E = 225*1e3
ν = 0.33
ρ = 9.13*1e-6
radius = 0.065
geom =  Geometry(radius, T)
mat = Material(E, ν, ρ, radius, T)

# read initial rotations for the beams
Re₀ = readdlm("examples/input_stent/Re0_positioning.txt")

# beams vector
beams = constructor_beams(nodes, connectivity, mat, geom, Re₀)

# -------------------------------------------------------------------------------------------
# Simulation parameters
# -------------------------------------------------------------------------------------------

# integration parameters
α = -0.05
β = 0.25*(1-α)^2
γ = 0.5*(1-2*α)
damping = 1e3

# time step and total time
Δt = 0.01
Δt_plot =  0.01
tᵉⁿᵈ = 1

# tolerance and maximum number of iterations
tol_res = 1e-5
tol_ΔD = 1e-5
max_it = 10

# Gauss points
nG = 3
ωG = Vec3(5/9, 8/9, 5/9)
zG = Vec3(-sqrt(3/5), 0, sqrt(3/5))

# contact parameters
kₙ = 10 #penalty parameter
μ = 0.3
εᵗ = 0.1 #regularized parameter for friction contact
ηₙ = 0.1

comp = SimulationParameters(α, β, γ, damping, Δt, Δt_plot, tᵉⁿᵈ, tol_res, tol_ΔD, max_it, nG, ωG, zG, kₙ, μ, εᵗ, ηₙ, T)

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
nodespairs = readdlm("examples/input_stent/constr_stent.txt")
cons = constructor_constraints(nodespairs, 1E3, 1)

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
conf = Configuration(mat, geom, ndofs, ext_forces, bcs, T)

# -------------------------------------------------------------------------------------------
# Start simulation
# -------------------------------------------------------------------------------------------

params = Params(output_dir = "examples/output3D", ENERGY_STOP = true, scale=2, SHOW_TIME_SECTIONS=false)
solver!(nodes, beams, conf, comp, sdf, cons, params, T)