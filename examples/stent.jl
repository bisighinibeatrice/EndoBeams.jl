using EndoBeams

# -------------------------------------------------------------------------------------------
# Read stent information
# -------------------------------------------------------------------------------------------

# read expanded configuration
posStentExp =  read_TXT_file_pos("examples/input_stent/pos_stent.txt")
connStent = read_TXT_file_conn("examples/input_stent/conn_stent.txt")

# read crimped configuration
posStentCrimped = copy(posStentExp + read_TXT_file_pos("examples/input_stent/u_crimping.txt"))

# read positioned configuration
posStentCrimpedPositioned =  posStentCrimped + read_TXT_file_pos("examples/input_stent/u_positioning.txt")

# -------------------------------------------------------------------------------------------
# Building the nodes
# -------------------------------------------------------------------------------------------

# total number of nodes
nnodes = size(posStentExp,1)

# initial conditions
u_0 = zeros(nnodes*3)
for (j,i) in enumerate(1:3:nnodes*3)
    u_0[i+0] = posStentCrimpedPositioned[j][1] - posStentExp[j][1]
    u_0[i+1] = posStentCrimpedPositioned[j][2] - posStentExp[j][2]
    u_0[i+2] = posStentCrimpedPositioned[j][3] - posStentExp[j][3]
end  

udt_0 = zeros(nnodes*3)
udtdt_0 = zeros(nnodes*3)
w_0 =  zeros(nnodes*3)
wdt_0 = zeros(nnodes*3)
wdtdt_0 = zeros(nnodes*3)
R_0 = read_TXT_file_ICs_matrix("examples/input_stent/R_positioning.txt")

# plane for cylindrical coordinates
plane = fill("xy", length(posStentExp))

# nodes StructArray
allnodes = constructor_nodes(posStentExp, u_0, udt_0, udtdt_0, w_0, wdt_0, wdtdt_0, plane, R_0)

# -------------------------------------------------------------------------------------------
# Building the beams
# -------------------------------------------------------------------------------------------

# total number of beams
nbeams = nnodes-1

# interpolation points per beam
nbInterpolationPoints = 30

# geometric and material properties
E = 225*1e3
ν = 0.33
ρ = 9.13*1e-6
rWireSection = 0.065
geom =  constructor_geometry_properties(rWireSection)
mat = constructor_material_properties(E, ν, ρ, rWireSection)

# read initial rotations for the beams
Re_0 = read_TXT_file_ICs_matrix("examples/input_stent/Re0_positioning.txt")

# beams vector
allbeams = constructor_beams(allnodes, connStent, mat, geom, nbInterpolationPoints, Re_0)

# -------------------------------------------------------------------------------------------
# Simulation parameters
# -------------------------------------------------------------------------------------------

# integration parameters
alpha = -0.05
beta = 0.25*(1-alpha)^2
gamma = 0.5*(1-2*alpha)
damping = 1E4*5

# time step and total time
dt = 0.1
dt_plot =  0.1
tend = 1E10

# tolerance and maximum number of iterations
tol_res = 1e-5
tol_ddk = 1e-5
max_it = 10

# Gaussian points
nG = 3
wG = Vec3(5/9, 8/9, 5/9)
zG = Vec3(-sqrt(3/5), 0, sqrt(3/5))

# contact parameters
eps_C = 500 #penalty parameter
mu_T = 0.01
eps_tol_fric = 0.1 #regularized parameter for friction cont act

comp = constructor_simulation_parameters(alpha, beta, gamma, damping, dt, dt_plot, tend, tol_res, tol_ddk, max_it, nG, wG, zG, eps_C, mu_T, eps_tol_fric)

# -------------------------------------------------------------------------------------------
# External forces
# -------------------------------------------------------------------------------------------

# external force and applied dof
flag_crimping = false
Fext(t) = 0*t
dofs_load = T[]

ext_forces = constructor_ext_force(flag_crimping, Fext, dofs_load)

# -------------------------------------------------------------------------------------------
# Boundary conditions
# -------------------------------------------------------------------------------------------

# number of dof (6 per node)
ndofs = nnodes*6

# lagrange constrains
nodespairs = read_TXT_file_conn("examples/input_stent/constr_stent.txt")
cons = constructor_constraints(nodespairs, 1E6, 1E3)

# Dirichlet boundary conditions: blocked positions
fixed_dofs = T[]
free_dofs = setdiff(1:ndofs, fixed_dofs)

# Dirichlet dof (x6)
flag_cylindrical = false
dofs_disp = T[]
Fdisp(t) = 0
flag_disp_vector = false
udisp = T[]

# boundary conditions strucutre
bc = constructor_boundary_conditions(fixed_dofs, free_dofs, flag_cylindrical, flag_disp_vector, Fdisp, udisp, dofs_disp)

# -------------------------------------------------------------------------------------------
# SDF
# -------------------------------------------------------------------------------------------

sdf = constructor_discrete_sdf("examples/input_stent/arcStretchObj.vtk", rWireSection, false)

# -------------------------------------------------------------------------------------------
# Final configuration
# -------------------------------------------------------------------------------------------

# configuration: mesh, external forces and boundary conditions
conf = constructor_configuration(mat, geom, nnodes, ndofs, ext_forces, bc)

# -------------------------------------------------------------------------------------------
# Start simulation
# -------------------------------------------------------------------------------------------

params = Params(thisDirOutputPath = "examples/output3D", ENERGY_STOP = true, SAVE_ENERGY = true)
solver!(allnodes, allbeams, conf, comp, sdf, cons, params)