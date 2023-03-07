using EndoBeams



#-------------------------------------------------------------------------------------------
# Building the nodes
# -------------------------------------------------------------------------------------------

# positions
dx = 1; L = 10
aux = 0:dx:L

X₀ = [Vec3(0, x-4.31, 3.23) for x in 0:dx:L]
append!(X₀, [Vec3(0.06, 0, x) for x in 0:dx:L])
X₀  = X₀./100

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
nodes = build_nodes(X₀, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, nothing)

# -------------------------------------------------------------------------------------------
# Building the beams
# -------------------------------------------------------------------------------------------

# total number of beams
nbeams = nnodes-1

connectivity = [Vec2(i, i+1) for i in 1:Int(nnodes*0.5)-1]
append!(connectivity, [Vec2(i, i+1) for i in Int(nnodes*0.5)+1:nnodes-1])

# geometric and material properties
E = 5*1e9
ν = 0.33
ρ = 7850
radius = 0.3/1000
damping = 0

beams = build_beams(nodes, connectivity, E, ν, ρ, radius, damping)

# -------------------------------------------------------------------------------------------
# Contact parameters
# -------------------------------------------------------------------------------------------

# contact parameters
kₙ = 0.1 #penalty parameter
μ = 0.0
εᵗ = 0.1 #regularized parameter for friction contact
ηₙ = 0.1
kₜ = kₙ
ηₜ = ηₙ
u̇ₛ = 0.005
beam2beam = true

contact = ContactParameters(kₙ, μ, εᵗ, ηₙ, kₜ, ηₜ, u̇ₛ, beam2beam)

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

# multifreedom constraints
constraints = nothing

# Dirichlet boundary conditions: fixed positions
ndofs = nnodes*6
fixed_dofs = Int[1:6; 6*(nnodes/2+1-1).+(1:6)]
free_dofs = setdiff(1:ndofs, fixed_dofs)

# Dirichlet dof (x6)
disp_dofs = [1]
disp_vals = zeros(size(disp_dofs))
dispA = 10/100
tA = 2.5
k = dispA/tA
disp(t, node_idx) = (k*t).*(t<=tA) + (-k*t + 2*dispA).*(t>tA && t<=2*tA) + (k*(t-2*tA)).*(t>2*tA && t<=3*tA) +(-k*t + 4*dispA).*(t>3*tA && t<=4*tA)

# boundary conditions strucutre 
bcs = BoundaryConditions(fixed_dofs, free_dofs, disp, disp_vals, disp_dofs)

# -------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------

conf = Configuration(nodes, beams, nothing, ext_forces, bcs, contact, nothing)

# -------------------------------------------------------------------------------------------
# Time stepping parameters
# -------------------------------------------------------------------------------------------

# initial time step and total time
ini_Δt = 0.0001
max_Δt = 0.01
Δt_plot = 0.01
tᵉⁿᵈ = 1

params = Params(;ini_Δt, max_Δt, Δt_plot, tᵉⁿᵈ, output_dir = "examples/output3D", verbose=false , record_timings=true)

# -------------------------------------------------------------------------------------------
# Solve
# -------------------------------------------------------------------------------------------

solver!(conf, params);
