using EndoBeams
using DelimitedFiles
using LinearAlgebra

#-------------------------------------------------------------------------------------------
# Read the mesh
# -------------------------------------------------------------------------------------------

# positions
X₀ =  readdlm("examples/input_net/pos_net.txt")
nnodes = size(X₀, 1)

# connectivity
connectivity = readdlm("examples/input_net/conn_net.txt", Int)

#-------------------------------------------------------------------------------------------
# Building the nodes
# -------------------------------------------------------------------------------------------

# initial conditions
u⁰ = zeros(nnodes, 3)
u̇⁰ = zeros(nnodes, 3)
for i in 1:nnodes
    u̇⁰[i, 2] = -0.001
end
ü⁰ = zeros(nnodes, 3)
w⁰ = zeros(nnodes, 3)
ẇ⁰ = zeros(nnodes, 3)
ẅ⁰ = zeros(nnodes, 3)


# nodes StructArray
nodes = build_nodes(X₀, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, nothing)

# -------------------------------------------------------------------------------------------
# Building the beams
# -------------------------------------------------------------------------------------------


# geometric and material properties
E = 10
ν = 0
ρ = 1500
radius = 0.2e-3
damping = 1e-3


beams = build_beams(nodes, connectivity, E, ν, ρ, radius, damping)



# contact parameters
kₙ = 0.1 #penalty parameter
μ = 0.1
εᵗ = 0.1 #regularized parameter for friction contact
ηₙ = 0.1

contact = ContactParameters(kₙ, μ, εᵗ, ηₙ)
  
# -------------------------------------------------------------------------------------------
# External forces
# -------------------------------------------------------------------------------------------

# external force and applied dof
loaded_dofs = Float64[]
force(t, node_idx) = 0

ext_forces = ExternalForces(force, loaded_dofs)

# -------------------------------------------------------------------------------------------
# Boundary conditions
# -------------------------------------------------------------------------------------------

# number of dof (6 per node)
ndofs = nnodes*6

# Dirichlet boundary conditions: blocked positions
fixed_dofs = Float64[]
free_dofs = setdiff(1:ndofs, fixed_dofs)

# Dirichlet dof (x6)
disp_dofs = Int[]
disp_vals = Float64[]
disp(t, node_idx) = 0


# boundary conditions strucutre
bcs = BoundaryConditions(fixed_dofs, free_dofs, disp, disp_vals, disp_dofs)
# -------------------------------------------------------------------------------------------
# SDF
# -------------------------------------------------------------------------------------------

# analytical SDF
sdf = Sphere_SDF(radius, 0.05, 0, 0, 0)

# -------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------

# configuration: mesh, external forces and boundary conditions
conf = Configuration(nodes, beams, nothing, ext_forces, bcs, contact, sdf)

# -------------------------------------------------------------------------------------------
# Time stepping parameters
# -------------------------------------------------------------------------------------------

# initial time step and total time
ini_Δt = 0.5
max_Δt = 0.5
Δt_plot =  0.5
tᵉⁿᵈ = 100

params = Params(;ini_Δt, Δt_plot, max_Δt, tᵉⁿᵈ, output_dir = "examples/output3D")



# -------------------------------------------------------------------------------------------
# Solve
# -------------------------------------------------------------------------------------------


solver!(conf, params);