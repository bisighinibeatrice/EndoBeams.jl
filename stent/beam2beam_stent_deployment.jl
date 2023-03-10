using EndoBeams
using Parameters
using DelimitedFiles

include("utils_stent.jl")

# -------------------------------------------------------------------------------------------
# Building the stent
# -------------------------------------------------------------------------------------------

@unpack nbWires, rStent, rCrimpedStent, rWireSection, wireGap, lengthStent, nbTotalCells, braidingPattern = WallStent()
X₀mat, connectivity = compute_bs_geom_given_nbTotalCells(nbWires, rStent, rWireSection, wireGap, lengthStent, nbTotalCells, braidingPattern)
X₀ = reshape(reinterpret(Float64, X₀mat), (3, length(X₀mat)))'

# total number of nodes
nnodes = size(X₀, 1)

# initial conditions
u⁰ = readdlm("stent/output3D/outputCrimping3D/" * "u.txt")
u̇⁰ = zeros(nnodes, 3)
ü⁰ = zeros(nnodes, 3)
w⁰ = zeros(nnodes, 3)
ẇ⁰ = zeros(nnodes, 3)
ẅ⁰ = zeros(nnodes, 3)

R⁰ = readdlm("stent/output3D/outputCrimping3D/" * "R.txt") 

# nodes StructArray
nodes = build_nodes(X₀, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, nothing, R⁰)

# total number of beams
nbeams = nnodes-1

# geometric and material properties
mass_scaling = 1e4
E = 225*1e3
ν = 0.33
ρ = 9.13*1e-6 * mass_scaling
damping = 1E4

Re₀ = read_ics_mat(readdlm("stent/output3D/outputCrimping3D/" * "Re0.txt"))

beams = build_beams(nodes, connectivity, E, ν, ρ, rWireSection, damping, Re₀)

# -------------------------------------------------------------------------------------------
# Contact parameters
# -------------------------------------------------------------------------------------------

# contact parameters
kₙ = 4/3 * 5/(1-0.5^2)*sqrt(rWireSection) # Approximate Hertz contact with 5 MPa wall stiffness
μ = 0.01
εᵗ = 0.001 #regularized parameter for friction contact
ηₙ = 0.1
kₜ = kₙ
ηₜ = ηₙ
u̇ₛ = 0.005
beam2beam = true

contact = ContactParameters(kₙ, μ, εᵗ, ηₙ, kₜ, ηₜ, u̇ₛ, beam2beam)
# sdf = Cylinder_SDF(rWireSection, 2)
sdf = nothing

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

# Dirichlet boundary conditions: fixed positions
ndofs = nnodes*6
fixed_dofs = Int[]
free_dofs = setdiff(1:ndofs, fixed_dofs)

# Dirichlet dof (x6)
disp_dofs = Int[]
disp_vals = Float64[]
disp(t, node_idx) = 0

# boundary conditions strucutre 
bcs = BoundaryConditions(fixed_dofs, free_dofs, disp, disp_vals, disp_dofs)

# -------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------

conf = Configuration(nodes, beams, nothing, ext_forces, bcs, contact, sdf)

# -------------------------------------------------------------------------------------------
# Time stepping parameters
# -------------------------------------------------------------------------------------------

# initial time step and total time
ini_Δt = 1e-6
max_Δt = 100
Δt_plot = 0
tᵉⁿᵈ = 10000000

params = Params(;ini_Δt, max_Δt, Δt_plot, tᵉⁿᵈ, output_dir = "stent/output3D/outputDeployment3D", verbose=true , record_timings=true)

# -------------------------------------------------------------------------------------------
# Solve
# -------------------------------------------------------------------------------------------

solver!(conf, params);
s
