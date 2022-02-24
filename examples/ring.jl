using EndoBeams

const T = Float64

# -------------------------------------------------------------------------------------------
# Building the nodes
# -------------------------------------------------------------------------------------------

# external and internal raidus of the ring
Rmid = 10
nelem = 32
alpha_div = 2*pi/(nelem)

# positions
X₀ =  Vec3{T}[]

push!(X₀, Vec3(Rmid, 0, 0))

for idiv = 1:nelem-1
    alphai = idiv*alpha_div
    Xi=Rmid*cos(alphai)
    Zi=Rmid*sin(alphai)
    push!(X₀, Vec3(Xi,  0, Zi))
end 

# total number of nodes
nnodes = size(X₀,1)

# initial conditions
u⁰ = zeros(size(X₀,1)*3)

vx_ini = 2*sqrt(2)/2
vz_ini = -2*sqrt(2)/2
u̇⁰ = zeros(size(X₀,1)*3,1)
u̇⁰[1:3:end-2] .= vx_ini
u̇⁰[3:3:end] .= vz_ini

ü⁰ = zeros(size(X₀,1)*3)
w⁰ = zeros(size(X₀,1)*3)
ẇ⁰ = zeros(size(X₀,1)*3)
ẅ⁰ = zeros(size(X₀,1)*3)

# plane for cylindrical coordinates
plane = fill("xy", length(X₀))

# nodes StructArray
nodes = constructor_nodes(X₀, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, nothing, T)

# -------------------------------------------------------------------------------------------
# Building the beams
# -------------------------------------------------------------------------------------------

# total number of beams
nbeams = nnodes # the ring is closed

# conn
conn = Vec2{Int}[]
aux1 = 1:nbeams-1
aux2 = 2:nbeams

for i in 1:nbeams-1
    push!(conn, (aux1[i], aux2[i])) 
end

push!(conn, (conn[1][1], conn[end][2]))

# interpolation points per beam
nbInterpolationPoints = 30

# geometric dimension of the cylinder
r = 0.5
E = 100
nu = 0.0001
ρ = 0.01
G = E/(2*(1+nu))
A = pi*r^2
I₂₂ = 0.25*pi*r^4
I₃₃ = 0.25*pi*r^4
Iₒ = I₂₂+I₃₃
Iᵣᵣ = Iₒ
J = Iₒ
Jᵨ = Mat33(ρ*Iₒ, 0, 0, 0, ρ*I₂₂, 0, 0, 0, ρ*I₃₃)
Aᵨ = ρ*A

geom = Geometry{T}(A, I₂₂, I₃₃, Iₒ, Iᵣᵣ, J)
mat = Material{T}(E, G, Aᵨ, Jᵨ)

# beams vector
beams = constructor_beams(nodes, conn, mat, geom, nbInterpolationPoints, nothing)

#-----------------------------------------------------------------------------------
# Simulation parameters
# -------------------------------------------------------------------------------------------

# integration parameters
α = -0.05
β = 0.25*(1-α)^2
γ = 0.5*(1-2*α)
damping = 0

# time step and total time
Δt = 0.2/16
Δt_plot = 0.2/16
tᵉⁿᵈ = 20

# tolerance and maximum number of iterations
tol_res = 1e-5
tol_ΔD = 1e-5
max_it = 10

# Gauss points
nG = 3
ωG = Vec3(5/9, 8/9, 5/9)
zG = Vec3(-sqrt(3/5), 0, sqrt(3/5))

# penalty parameters
εᶜ = 1E6
μ = 0.3
εᵗ = 0.5

comp = SimulationParameters(α, β, γ, damping,  Δt, Δt_plot, tᵉⁿᵈ, tol_res, tol_ΔD, max_it, nG, ωG, zG, εᶜ, μ, εᵗ, T)

# -------------------------------------------------------------------------------------------
# External forces
# -------------------------------------------------------------------------------------------

# external force and applied dof
flag_crimping = false
F(t) = 0
dofs_load = T[]

ext_forces = ExternalForces(flag_crimping, F, dofs_load, T)

# -------------------------------------------------------------------------------------------
# Boundary conditions
# -------------------------------------------------------------------------------------------

# multifreedom constraints
cons = Int[]

# number of dof (6 per node)
ndofs = nnodes*6

# Dirichlet boundary conditions: fixed positions
fixed_dofs = Int[]
free_dofs = setdiff(1:ndofs, fixed_dofs) 

# Dirichlet boundary conditions: moving positions
flag_cylindrical = false
disp_dofs = Int[]
Fdisp(t) = 0
flag_disp_vector = false
disp_vals = T[]

# boundary conditions strucutre 
bcs = BoundaryConditions(fixed_dofs, free_dofs, flag_cylindrical, flag_disp_vector, Fdisp, disp_vals, disp_dofs, T)

# -------------------------------------------------------------------------------------------
# SDF
# -------------------------------------------------------------------------------------------

r = 0.5
z0 = -12

sdf = SDF_Plane_z{T}(r, z0)

# -------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------


# configuration: mesh, external forces and boundary conditions
conf = Configuration(mat, geom, nnodes, ndofs, ext_forces, bcs, T)

# -------------------------------------------------------------------------------------------
# Solve
# -------------------------------------------------------------------------------------------

params = Params(scale = 2, output_dir = "examples/output3D")
solver!(nodes, beams, conf, comp, sdf, cons, params, T)                                        
