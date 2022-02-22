using EndoBeams

const T = Float64

#-------------------------------------------------------------------------------------------
# Read the mesh
# -------------------------------------------------------------------------------------------

# positions
X₀ =  read_TXT_file_pos("examples/input_net/pos_net.txt")
nnodes = length(X₀)

# connectivity
conn = read_TXT_file_conn("examples/input_net/conn_net.txt")

#-------------------------------------------------------------------------------------------
# Building the nodes
# -------------------------------------------------------------------------------------------

# initial conditions
u⁰ = zeros(nnodes*3)
u̇⁰ = zeros(nnodes*3)
for i in 2:3:nnodes*3
    u̇⁰[i] = -0.001
end
ü⁰ = zeros(nnodes*3)
w⁰ = zeros(nnodes*3)
ẇ⁰ = zeros(nnodes*3)
ẅ⁰ = zeros(nnodes*3)

# plane for cylindrical coordinates
plane = fill("xy", length(X₀))

# nodes StructArray
nodes = constructor_nodes(X₀, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, plane, nothing, T)

# -------------------------------------------------------------------------------------------
# Building the beams
# -------------------------------------------------------------------------------------------

# interpolation points per beam
nbInterpolationPoints = 30

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
Jᵨ = Mat33(ρ*J, 0, 0, 0, ρ*I₂₂, 0, 0, 0, ρ*I₃₃)
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
damping = 1e-3

# time step and total time
Δt = 0.5
Δt_plot = 0.5
tᵉⁿᵈ = 100

# tolerance and maximum number of iterations
res_tol = 1e-5
tol_ddk = 1e-5
max_it = 10

# Gauss points
nG = 3
ωG = Vec3(5/9, 8/9, 5/9)
zG = Vec3(-sqrt(3/5), 0, sqrt(3/5))

# penalty parameters
εᶜ = 0.01
μ = 0.01
εᵗ = 0.5

comp = constructor_simulation_parameters(α, β, γ, damping,  Δt, Δt_plot, tᵉⁿᵈ, res_tol, tol_ddk, max_it, nG, ωG, zG, εᶜ, μ, εᵗ, T)

# -------------------------------------------------------------------------------------------
# External forces
# -------------------------------------------------------------------------------------------

# external force and applied dof
flag_crimping = false
Fext(t) = 0*t
dofs_load = Int[]

ext_forces = constructor_ext_force(flag_crimping, Fext, dofs_load)

# -------------------------------------------------------------------------------------------
# Boundary conditions
# -------------------------------------------------------------------------------------------

# multifreedom constrains
cons = nothing

# Dirichlet boundary conditions: fixed positions
ndofs = nnodes*6
fixed_dofs = T[]
free_dofs = setdiff(1:ndofs, fixed_dofs)

# Dirichlet boundary conditions: moving positions
flag_cylindrical = false
Fdisp(t) = 0
dofs_disp = Int[]
flag_disp_vector = false
udisp = T[]

# boundary conditions strucutre
bcs = constructor_boundary_conditions(fixed_dofs, free_dofs, flag_cylindrical, flag_disp_vector, Fdisp, udisp, dofs_disp, T)

# -------------------------------------------------------------------------------------------
# SDF
# -------------------------------------------------------------------------------------------

# analytical SDF
sdf = SDF_Sphere{T}(r, 0.05, 0, 0, 0)

# discrete SDF
# inside = true
# sdf = constructor_discrete_sdf("examples/input_net/sdf_sphere_20.vtk", r, inside)

# -------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------

# configuration: mesh, external forces and boundary conditions
conf = constructor_configuration(mat, geom, nnodes, ndofs, ext_forces, bcs, T)

# -------------------------------------------------------------------------------------------
# Solve
# -------------------------------------------------------------------------------------------

params = Params(scale =2,thisDirOutputPath = "examples/output3D", SHOW_TIME_SECTIONS=true)
solver!(nodes, beams, conf, comp, sdf, cons, params, T)