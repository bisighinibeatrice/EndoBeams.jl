using EndoBeams
using Parameters

include("utils_stent.jl")

# -------------------------------------------------------------------------------------------
# Building the stent
# -------------------------------------------------------------------------------------------

@unpack nbWires, rStent, rCrimpedStent, rWireSection, wireGap, lengthStent, nbTotalCells, braidingPattern = BraidedStentB2B()
X₀, connectivity = compute_bs_geom_given_nbTotalCells(nbWires, rStent, rWireSection, wireGap, lengthStent, nbTotalCells, braidingPattern)

write_vtk_configuration("initial_position_stent.vtk", X₀, connectivity)

# total number of nodes
nnodes = length(X₀)

# initial conditions
u⁰ = zeros(Vec3, nnodes)
u̇⁰ = zeros(Vec3, nnodes)
ü⁰ = zeros(Vec3, nnodes)
w⁰ = zeros(Vec3, nnodes)
ẇ⁰ = zeros(Vec3, nnodes)
ẅ⁰ = zeros(Vec3, nnodes)

# nodes StructArray
nodes = build_nodes(X₀, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, nothing)

# total number of beams
nbeams = nnodes-1

# geometric and material properties
E = 5*1e9
ν = 0.33
ρ = 7850
damping = 0

beams = build_beams(nodes, connectivity, E, ν, ρ, rWireSection, damping)

# -------------------------------------------------------------------------------------------
# Contact parameters
# -------------------------------------------------------------------------------------------

# contact parameters
kₙ = 0.1 #penalty parameter
μ = 0.0
εᵗ = 0.1 #regularized parameter for friction contact
ηₙ = 0.1
kₜ = kₙ
ηₜ = ηₙ
u̇ₛ = 0.005
beam2beam = true

contact = ContactParameters(kₙ, μ, εᵗ, ηₙ, kₜ, ηₜ, u̇ₛ, beam2beam)

# -------------------------------------------------------------------------------------------
# External forces
# -------------------------------------------------------------------------------------------

# external force and applied dof
loaded_dofs = []
force(t, node_idx) = 0 
   
ext_forces = ExternalForces(force, loaded_dofs)

# -------------------------------------------------------------------------------------------
# Boundary conditions
# -------------------------------------------------------------------------------------------

# multifreedom constraints
constraints = nothing

# Dirichlet boundary conditions: fixed positions
ndofs = nnodes*6
fixed_dofs = Int[]
free_dofs = setdiff(1:ndofs, fixed_dofs)

# Dirichlet dof (x6)
disp_dofs = []
disp_vals = []
disp(t, node_idx) = 0

# boundary conditions strucutre 
bcs = BoundaryConditions(fixed_dofs, free_dofs, disp, disp_vals, disp_dofs)

# -------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------

conf = Configuration(nodes, beams, nothing, ext_forces, bcs, contact, nothing)

# -------------------------------------------------------------------------------------------
# Time stepping parameters
# -------------------------------------------------------------------------------------------

# initial time step and total time
ini_Δt = 0.0001
max_Δt = 0.01
Δt_plot = 0.01
tᵉⁿᵈ = 1

params = Params(;ini_Δt, max_Δt, Δt_plot, tᵉⁿᵈ, output_dir = "examples/output3D", verbose=false , record_timings=true)

# -------------------------------------------------------------------------------------------
# Solve
# -------------------------------------------------------------------------------------------

solver!(conf, params);

