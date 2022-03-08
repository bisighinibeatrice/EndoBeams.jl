using EndoBeams

const T = Float64

#-------------------------------------------------------------------------------------------
# Building the nodes
# -------------------------------------------------------------------------------------------

# positions
dx = 2.5; L = 10

X₀ = [Vec3(0,x,0) for x in 0:dx:L]
append!(X₀, [Vec3(-x,L,0) for x in dx:dx:L] )

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
nodes = nodes(X₀, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, nothing, T)

# -------------------------------------------------------------------------------------------
# Building the beams
# -------------------------------------------------------------------------------------------

# total number of beams
nbeams = nnodes-1

connectivity = [Vec2(i, i+1) for i in 1:nbeams]


# geometric and material properties
E = 1e6
G = 1e6
Jᵨ = EndoBeams.Diagonal(Vec3(20, 10, 10))
Aᵨ = 1
A = 1
I₂₂ = 1e-3
I₃₃ = 1e-3
Iₒ = 0
Iᵣᵣ = 0
J = 1e-3

geometry = Geometry{T}(A, I₂₂, I₃₃, Iₒ, Iᵣᵣ, J)
material = Material{T, typeof(Jᵨ)}(E, G, Aᵨ, Jᵨ)

# beams vector
beams = beams(nodes, connectivity, material, geometry, nothing)

#-----------------------------------------------------------------------------------
# Simulation parameters
# -------------------------------------------------------------------------------------------

# integration parameters
α = -0.05
β = 0.25*(1-α)^2
γ = 0.5*(1-2*α)
damping = 0

# time step and total time
Δt = 0.25
Δt_plot = 0.25
tᵉⁿᵈ = 30

# tolerance and maximum number of iterations
tol_res = 1e-5
tol_ΔD = 1e-5
max_it = 10

# Gauss points
nᴳ = 3
ωᴳ = Vec3(5/9, 8/9, 5/9)
zᴳ = Vec3(-sqrt(3/5), 0, sqrt(3/5)) 

# penalty parameters
kₙ = 5000
μ = 0
εᵗ = 0.1

comp = SimulationParameters(α, β, γ, damping,  Δt, Δt_plot, tᵉⁿᵈ, tol_res, tol_ΔD, max_it, nᴳ, ωᴳ, zᴳ, kₙ, μ, εᵗ, T)

# -------------------------------------------------------------------------------------------
# External forces
# -------------------------------------------------------------------------------------------

# external force and applied dof
loaded_dofs = [6*length(0:dx:L)-3]
force(t, node_idx) = 1*(t.<=1).*(50*t) .+1*((t.>1) .& (t.<=2)).*(-50*t.+100).+(t.>2).*0

ext_forces = ExternalForces(force, loaded_dofs)

# -------------------------------------------------------------------------------------------
# Boundary conditions
# -------------------------------------------------------------------------------------------

# multifreedom constraints
constraints = nothing

# Dirichlet boundary conditions: fixed positions
ndofs = nnodes*6
fixed_dofs = 1:6
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

#there is no contact
sdf = nothing

# -------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------

# configuration: mesh, external forces and boundary conditions
conf = Configuration(material, geometry, ndofs, ext_forces, bcs, T)

# -------------------------------------------------------------------------------------------
# Solve
# -------------------------------------------------------------------------------------------

params = Params(output_dir = "examples/output3D", record_timings=false)
solver!(nodes, beams, conf, comp, sdf, constraints, params, T)       