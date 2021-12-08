using EndoBeams

T = Float64

#-------------------------------------------------------------------------------------------
# Building the nodes
# -------------------------------------------------------------------------------------------

# positions
dx = 2.5; L = 10

pos =  Vec3{T}[]
for x in 0:dx:L
    push!(pos, Vec3(0, x, 0))
end
for x in dx:dx:L
    push!(pos, Vec3(-x,  L, 0))
end

# total number of nodes
nnodes = size(pos,1)

# initial conditions
u_0 = zeros(nnodes*3)
udt_0 = zeros(nnodes*3)
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

# total number of beams
nbeams = nnodes-1

# conn
conn = Vec2{Int}[]
aux1 =  1:nnodes-1
aux2 =  2:nnodes

for i in 1:nbeams
    push!(conn, (aux1[i], aux2[i])) 
end

# interpolation points per beam
nbInterpolationPoints = 30

# geometric and material properties
E = 1e6
G = 1e6
Jrho = Mat33(20, 0, 0, 0, 10, 0, 0, 0, 10)
Arho = 1
A = 1
I22 = 1e-3
I33 = 1e-3
Io = 0
Irr = 0
J = 1e-3

geom = Geometry{T}(A, I22, I33, Io, Irr, J)
mat = Material{T}(E, G, Arho, Jrho)

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
dt = 0.25
dt_plot = 0.25
tend = 30

# tolerance and maximum number of iterations
tol_res = 1e-5
tol_ddk = 1e-5
max_it = 10

# Gauss points
nG = 3
wG = Vec3(5/9, 8/9, 5/9)
zG = Vec3(-sqrt(3/5), 0, sqrt(3/5)) 

# penalty parameters
eps_C = 5000
mu_T = 0
eps_tol_fric = 0.1

comp = constructor_simulation_parameters(alpha, beta, gamma, damping,  dt, dt_plot, tend, tol_res, tol_ddk, max_it, nG, wG, zG, eps_C, mu_T, eps_tol_fric, T)

# -------------------------------------------------------------------------------------------
# External forces
# -------------------------------------------------------------------------------------------

# external force and applied dof
flag_crimping = false
Fext(t) = 1*(t.<=1).*(50*t) .+1*((t.>1) .& (t.<=2)).*(-50*t.+100).+(t.>2).*0
dofs_load = Int[]
push!(dofs_load, 6*size((0:dx:L),1)-3)

ext_forces = constructor_ext_force(flag_crimping, Fext, dofs_load, T)

# -------------------------------------------------------------------------------------------
# Boundary conditions
# -------------------------------------------------------------------------------------------

# multifreedom constrains
cons = T[]

# Dirichlet boundary conditions: fixed positions
ndofs = nnodes*6
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

params = Params(thisDirOutputPath = "examples/output3D")
solver!(allnodes, allbeams, conf, comp, sdf, cons, params, T)       