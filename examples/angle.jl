using EndoBeams



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
nodes = build_nodes(X₀, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, nothing)

# -------------------------------------------------------------------------------------------
# Building the beams
# -------------------------------------------------------------------------------------------

# total number of beams
nbeams = nnodes-1

connectivity = [Vec2(i, i+1) for i in 1:nbeams]


# geometric and material properties
E = 1e6
ν = 0.3
ρ = 1
radius = 0.3
damping = 0

beams = build_beams(nodes, connectivity, E, ν, ρ, radius, damping)



# -------------------------------------------------------------------------------------------
# External forces
# -------------------------------------------------------------------------------------------

# external force and applied dof
loaded_dofs = [6*length(0:dx:L)-3]
function force(t, node_idx)
    if t<=1
        return 50. *t
    elseif t>1 && t<=2
        return -50. *t + 100.
    else
        return 0.
    end
end

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
disp_vals = Float64[]
disp(t, node_idx) = 0

# boundary conditions strucutre 
bcs = BoundaryConditions(fixed_dofs, free_dofs, disp, disp_vals, disp_dofs)

# -------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------

conf = Configuration(nodes, beams, nothing, ext_forces, bcs, nothing, nothing)

# -------------------------------------------------------------------------------------------
# Time stepping parameters
# -------------------------------------------------------------------------------------------

# initial time step and total time
ini_Δt = 0.25
max_Δt = 0.25
Δt_plot =  0.25
tᵉⁿᵈ = 30

params = Params(;ini_Δt, max_Δt, Δt_plot, tᵉⁿᵈ, output_dir = "examples/output3D")


# -------------------------------------------------------------------------------------------
# Solve
# -------------------------------------------------------------------------------------------


solver!(conf, params);