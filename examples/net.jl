using EndoBeams
using DelimitedFiles
using LinearAlgebra

const T = Float64

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
nodes = constructor_nodes(X₀, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, nothing, T)

# -------------------------------------------------------------------------------------------
# Building the beams
# -------------------------------------------------------------------------------------------


# geometric and material properties
r = 0.2/1000
E = 10
vu = 0.0
G = E/(2*(1+vu))
ρ = 1500
A = pi*r^2
I₂₂ = 0.25*pi*r^4
I₃₃ = 0.25*pi*r^4
Iₒ = I₂₂+I₃₃
Iᵣᵣ = I₂₂+I₃₃
J = I₂₂+I₃₃
Jᵨ = ρ*Diagonal(Vec3(J, I₂₂, I₃₃))
Aᵨ = ρ*A

geom = Geometry{T}(A, I₂₂, I₃₃, Iₒ, Iᵣᵣ, J)
mat = Material{T, typeof(Jᵨ)}(E, G, Aᵨ, Jᵨ)

# beams vector
beams = constructor_beams(nodes, connectivity, mat, geom, nothing)
    
#-----------------------------------------------------------------------------------
# Simulation parameters
# -------------------------------------------------------------------------------------------

# integration parameters
α = -0.05
β = 0.25*(1-α)^2
γ = 0.5*(1-2*α)
damping = 1e-3

# time step and total time
Δt = 0.5
Δt_plot = 0.5
tᵉⁿᵈ = 100

# tolerance and maximum number of iterations
tol_res = 1e-5
tol_ΔD = 1e-5
max_it = 10

# Gauss points
nG = 3
ωG = Vec3(5/9, 8/9, 5/9)
zG = Vec3(-sqrt(3/5), 0, sqrt(3/5))

# penalty parameters
kₙ = 0.01
μ = 0.01
εᵗ = 0.5
ηₙ = 0.01

comp = SimulationParameters(α, β, γ, damping,  Δt, Δt_plot, tᵉⁿᵈ, tol_res, tol_ΔD, max_it, nG, ωG, zG, kₙ, μ, εᵗ, ηₙ, T)

# -------------------------------------------------------------------------------------------
# External forces
# -------------------------------------------------------------------------------------------

# external force and applied dof
loaded_dofs = Int[]
force(t, node_idx) = 0

ext_forces = ExternalForces(force, loaded_dofs)

# -------------------------------------------------------------------------------------------
# Boundary conditions
# -------------------------------------------------------------------------------------------

# multifreedom constraints
cons = nothing


# Dirichlet boundary conditions: fixed positions
ndofs = nnodes*6
fixed_dofs = Int[]
free_dofs = 1:ndofs


# Dirichlet dof (x6)
disp_dofs = Int[]
disp_vals = T[]
disp(t, node_idx) = 0

# boundary conditions strucutre 
bcs = BoundaryConditions(fixed_dofs, free_dofs, disp, disp_vals, disp_dofs, T)

# -------------------------------------------------------------------------------------------
# SDF
# -------------------------------------------------------------------------------------------

# analytical SDF
sdf = Sphere_SDF{T}(r, 0.05, 0, 0, 0)

# discrete SDF
# inside = true
# sdf = Discrete_SDF("examples/input_net/sdf_sphere_20.vtk", r, inside)

# -------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------

# configuration: mesh, external forces and boundary conditions
conf = Configuration(mat, geom, ndofs, ext_forces, bcs, T)

# -------------------------------------------------------------------------------------------
# Solve
# -------------------------------------------------------------------------------------------

params = Params(scale =2,output_dir = "examples/output3D", SHOW_TIME_SECTIONS=true)
solver!(nodes, beams, conf, comp, sdf, cons, params, T)