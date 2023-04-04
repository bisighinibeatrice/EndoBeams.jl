using EndoBeams
using Parameters
using DelimitedFiles
using LinearAlgebra 

include("utils_stent.jl")

# -------------------------------------------------------------------------------------------
# Time stepping parameters
# -------------------------------------------------------------------------------------------

# initial time step and total time
ini_Δt = 1e-6
max_Δt = 1e-2
Δt_plot = 1e-2
tᵉⁿᵈ = 1

params = Params(;ini_Δt, max_Δt, Δt_plot, tᵉⁿᵈ, output_dir = "stent/output3D/outputCrimping3D", verbose=true , record_timings=true)

# -------------------------------------------------------------------------------------------
# Building the stent
# -------------------------------------------------------------------------------------------

@unpack nbWires, rStent, rCrimpedStent, rWireSection, wireGap, lengthStent, nbTotalCells, braidingPattern = WallStent()
X₀, connectivity = compute_bs_geom_given_nbTotalCells(nbWires, rStent, rWireSection, wireGap, lengthStent, nbTotalCells, braidingPattern)
connectivity = reshape(reinterpret(Int64, connectivity), (2, length(connectivity)))'

# ringInit = get_nodespairs_stent_first_ring(rWireSection, X₀) 

# for i in eachindex(ringInit)
#     replace!(connectivity, ringInit[i][2] => ringInit[i][1])
# end 

# for i in 1:length(X₀)-length(ringInit)

#     if !(i in connectivity)
#         connectivity[findall(x->x>=i, connectivity)] = connectivity[findall(x->x>=i, connectivity)]  .- 1
#         deleteat!(X₀, i)
#     end 

# end

# ringEnd = get_nodespairs_stent_last_ring(rWireSection, X₀) 
# ringEnd = reshape(reinterpret(Int64, ringEnd), (2, length(ringEnd)))'
# ringEnd = ringEnd[sortperm(ringEnd[:, 2]), :] 

# for i in 1:size(ringEnd,1)
#     replace!(connectivity, ringEnd[i,2] => ringEnd[i,1])
# end 

# for i in 1:length(X₀)-size(ringEnd,1)+1

#     if !(i in connectivity)
#         connectivity[findall(x->x>=i, connectivity)] = connectivity[findall(x->x>=i, connectivity)]  .- 1
#         deleteat!(X₀, i)
#     end 

# end

# total number of nodes
nnodes = length(X₀)

# initial conditions
u⁰ = zeros(Vec3, nnodes)
u̇⁰ = zeros(Vec3, nnodes)
ü⁰ = zeros(Vec3, nnodes)
w⁰ = zeros(Vec3, nnodes)
ẇ⁰ = zeros(Vec3, nnodes)
ẅ⁰ = zeros(Vec3, nnodes)

plane = "xy"

# nodes StructArray
nodes = build_nodes(X₀, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, plane)

# geometric and material properties
# E = 2.6 * 1E5 * 1E4
# ν = 0.33
# ρ = 8 * 1E-9 * 1E12
# damping = 15 
E = 2.6 * 1E5  * 1E3
ν = 0.33
ρ = 8 * 1E-9 * 1E6
damping = 1E2

beams = build_beams(nodes, connectivity, E, ν, ρ, rWireSection, damping)

# -------------------------------------------------------------------------------------------
# Contact parameters
# -------------------------------------------------------------------------------------------

# contact parameters
kₙ = 1E1 #penalty parameter
μ = 0.0
εᵗ = 0.1 #regularized parameter for friction contact
ηₙ = 0.1
kₜ = kₙ
ηₜ = ηₙ
u̇ₛ = 0.0
beam2beam = true

contact = ContactParameters(kₙ, μ, εᵗ, ηₙ, kₜ, ηₜ, u̇ₛ, beam2beam)

sdf = Cylinder_SDF(rWireSection, rStent*0.5)

# -------------------------------------------------------------------------------------------
# External forces
# -------------------------------------------------------------------------------------------
ndofs = nnodes*6

# external force and applied dof

loaded_dofs = Int[]
force(t, node_idx) = 0
   
ext_forces = ExternalForces(force, loaded_dofs, true)

# -------------------------------------------------------------------------------------------
# Boundary conditions
# -------------------------------------------------------------------------------------------

constraints = nothing

# Dirichlet boundary conditions: fixed positions
fixed_dofs = Int[]
free_dofs = setdiff(1:ndofs, fixed_dofs)

# Dirichlet dof (x6)
flag_cylindrical = false
disp_dofs = Int[]
disp_vals = zeros(ndofs)
dispTot = -(rStent-rCrimpedStent)
disp(t, node_idx) = 0

# boundary conditions strucutre 
bcs = BoundaryConditions(fixed_dofs, free_dofs, disp, disp_vals, disp_dofs, flag_cylindrical)

# -------------------------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------------------------

conf = Configuration(nodes, beams, constraints, ext_forces, bcs, contact, sdf)

# -------------------------------------------------------------------------------------------
# Solve
# -------------------------------------------------------------------------------------------

solver!(conf, params);
write_txt_solution(nodes, beams, "stent/output3D/outputCrimping3D")
