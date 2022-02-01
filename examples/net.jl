using EndoBeams

const T = Float64

#-------------------------------------------------------------------------------------------
# Read the mesh
# -------------------------------------------------------------------------------------------

# positions
pos =  read_TXT_file_pos("examples/input_net/pos_net.txt")
nnodes = length(pos)

# connectivity
conn = read_TXT_file_conn("examples/input_net/conn_net.txt")

#-------------------------------------------------------------------------------------------
# Building the nodes
# -------------------------------------------------------------------------------------------

# initial conditions
u_0 = zeros(nnodes*3)
udt_0 = zeros(nnodes*3)
for i in 2:3:nnodes*3
    udt_0[i] = -0.001
end
udtdt_0 = zeros(nnodes*3)
w_0 = zeros(nnodes*3)
wdt_0 = zeros(nnodes*3)
wdtdt_0 = zeros(nnodes*3)

# plane for cylindrical coordinates
plane = fill("xy", length(pos))

# nodes StructArray
allnodes = constructor_nodes(pos, u_0, udt_0, udtdt_0, w_0, wdt_0, wdtdt_0, plane, nothing, T)

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
rho = 1500
A = pi*r^2
I22 = 0.25*pi*r^4
I33 = 0.25*pi*r^4
Io = I22+I33
Irr = I22+I33
J = I22+I33
Jᵨ = Mat33(rho*J, 0, 0, 0, rho*I22, 0, 0, 0, rho*I33)
Arho = rho*A

geom = Geometry{T}(A, I22, I33, Io, Irr, J)
mat = Material{T}(E, G, Arho, Jᵨ)

# beams vector
allbeams = constructor_beams(allnodes, conn, mat, geom, nbInterpolationPoints, nothing, T)
    
#-----------------------------------------------------------------------------------
# Simulation parameters
# -------------------------------------------------------------------------------------------

# integration parameters
alpha = -0.05
beta = 0.25*(1-alpha)^2
gamma = 0.5*(1-2*alpha)
damping = 0

# time step and total time
dt = 0.5
dt_plot = 0.5
tend = 100

# tolerance and maximum number of iterations
tol_res = 1e-5
tol_ddk = 1e-5
max_it = 10

# Gauss points
nG = 3
ωG = Vec3(5/9, 8/9, 5/9)
zG = Vec3(-sqrt(3/5), 0, sqrt(3/5))

# penalty parameters
eps_C = 0.01
μ = 0.01
εₜ = 0.5

comp = constructor_simulation_parameters(alpha, beta, gamma, damping,  dt, dt_plot, tend, tol_res, tol_ddk, max_it, nG, ωG, zG, eps_C, μ, εₜ, T)

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
cons = T[]

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

params = Params(scale =2,thisDirOutputPath = "examples/output3D")
solver!(allnodes, allbeams, conf, comp, sdf, cons, params, T)