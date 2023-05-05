using EndoBeams

const T = Float64

#-------------------------------------------------------------------------------------------
# Building the nodes
# -------------------------------------------------------------------------------------------

# positions
L = 10
nelem = 20

pos = Vec3{T}[]
for x in range(0.0, L, length=nelem + 1)
    push!(pos, Vec3(x, 0, 0))
end

# total number of nodes
nnodes = size(pos, 1)

# initial conditions
u_0 = zeros(nnodes * 3)
udt_0 = zeros(nnodes * 3)
udtdt_0 = zeros(nnodes * 3)
w_0 = zeros(nnodes * 3)
wdt_0 = zeros(nnodes * 3)
wdtdt_0 = zeros(nnodes * 3)

# plane for cylindrical coordinates
plane = fill("xy", length(pos))

# nodes StructArray
allnodes = constructor_nodes(pos, u_0, udt_0, udtdt_0, w_0, wdt_0, wdtdt_0, plane, nothing, T)

# -------------------------------------------------------------------------------------------
# Building the beams
# -------------------------------------------------------------------------------------------

# total number of beams
nbeams = nnodes - 1

# conn
conn = Vec2{Int}[]
aux1 = 1:nnodes-1
aux2 = 2:nnodes

for i in 1:nbeams
    push!(conn, (aux1[i], aux2[i]))
end

# interpolation points per beam
nbInterpolationPoints = 30

E = 210 * 1e9               # [Pa]
Iy = Iz = 4.7619 * 1e-7     # [m^4]
A = 4.7619 * 1e-4           # [m^2]
G = 78 * 1e9                # [Pa]
It = 2 * Iy                 # [m^4]
# geometric and material properties
rho = 7750
Io = It
Irr = It
J = It
Arho = A*rho
Jrho = Mat33(rho*Io, 0, 0, 0, rho*Iy, 0, 0, 0, rho*Iz)

geom = Geometry{T}(A, Iz, Iy, Io, Irr, J)
mat = Material{T}(E, G, Arho, Jrho)

# beams vector
allbeams = constructor_beams(allnodes, conn, mat, geom, nbInterpolationPoints, nothing, T)

#-----------------------------------------------------------------------------------
# Simulation parameters
# -------------------------------------------------------------------------------------------

# integration parameters
alpha = -0.05
beta = 0.25 * (1 - alpha)^2
gamma = 0.5 * (1 - 2 * alpha)
damping = 0
dynamics = false

# time step and total time
dt = 1
dt_plot = 10
tend = 100

# tolerance and maximum number of iterations
tol_res = 1e-5
tol_ddk = 1e-5
max_it = 10

# Gauss points
nG = 3
wG = Vec3(5 / 9, 8 / 9, 5 / 9)
zG = Vec3(-sqrt(3 / 5), 0, sqrt(3 / 5))

# penalty parameters
eps_C = 5000
mu_T = 0
eps_tol_fric = 0.1

comp = constructor_simulation_parameters(alpha, beta, gamma, damping, dynamics, dt, dt_plot, tend, tol_res, tol_ddk, max_it, nG, wG, zG, eps_C, mu_T, eps_tol_fric, T)

# -------------------------------------------------------------------------------------------
# External forces
# -------------------------------------------------------------------------------------------

# external force and applied dof
flag_crimping = false
Fext(t) = 2pi * E * Iy/L / tend * t 
dofs_load = Int[]
push!(dofs_load, 6*(nelem+1))

ext_forces = constructor_ext_force(flag_crimping, Fext, dofs_load, T)

# -------------------------------------------------------------------------------------------
# Boundary conditions
# -------------------------------------------------------------------------------------------

# multifreedom constrains
cons = T[]

# Dirichlet boundary conditions: fixed positions
ndofs = nnodes * 6
fixed_dofs = 1:6
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

#there is no contact
sdf = nothing

# -------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------

# configuration: mesh, external forces and boundary conditions
conf = constructor_configuration(mat, geom, nnodes, ndofs, ext_forces, bcs, T)

# -------------------------------------------------------------------------------------------
# Solve
# -------------------------------------------------------------------------------------------

params = Params(thisDirOutputPath="examples/outputTipMoment")
solver!(allnodes, allbeams, conf, comp, sdf, cons, params, T)