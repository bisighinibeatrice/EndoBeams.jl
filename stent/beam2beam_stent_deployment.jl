using EndoBeams
using Parameters
using DelimitedFiles
using LinearAlgebra 

include("utils_stent.jl")

# -------------------------------------------------------------------------------------------
# Building the stent
# -------------------------------------------------------------------------------------------

@unpack nbWires, rStent, rCrimpedStent, rWireSection, wireGap, lengthStent, nbTotalCells, braidingPattern = WallStent()
X₀, connectivity = compute_bs_geom_given_nbTotalCells(nbWires, rStent, rWireSection, wireGap, lengthStent, nbTotalCells, braidingPattern)
connectivity = reshape(reinterpret(Int64, connectivity), (2, length(connectivity)))'

ringInit = get_nodespairs_stent_first_ring(rWireSection, X₀) 

for i in eachindex(ringInit)
    replace!(connectivity, ringInit[i][2] => ringInit[i][1])
end 

for i in 1:length(X₀)-length(ringInit)

    if !(i in connectivity)
        connectivity[findall(x->x>=i, connectivity)] = connectivity[findall(x->x>=i, connectivity)]  .- 1
        deleteat!(X₀, i)
    end 

end

ringEnd = get_nodespairs_stent_last_ring(rWireSection, X₀) 
ringEnd = reshape(reinterpret(Int64, ringEnd), (2, length(ringEnd)))'
ringEnd = ringEnd[sortperm(ringEnd[:, 2]), :] 

for i in 1:size(ringEnd,1)
    replace!(connectivity, ringEnd[i,2] => ringEnd[i,1])
end 

for i in 1:length(X₀)-size(ringEnd,1)+1

    if !(i in connectivity)
        connectivity[findall(x->x>=i, connectivity)] = connectivity[findall(x->x>=i, connectivity)]  .- 1
        deleteat!(X₀, i)
    end 

end

X₀mat = deepcopy(X₀)
X₀ = reshape(reinterpret(Float64, X₀), (3, length(X₀)))'

# total number of nodes
nnodes = size(X₀, 1)

# initial conditions
u⁰ = readdlm("stent/output3D/outputCrimping3D/" * "u.txt")
u̇⁰ = zeros(nnodes, 3)
ü⁰ = zeros(nnodes, 3)
w⁰ = zeros(nnodes, 3)
ẇ⁰ = zeros(nnodes, 3)
ẅ⁰ = zeros(nnodes, 3)

R⁰ = readdlm("stent/output3D/outputCrimping3D/" * "R.txt") 

# nodes StructArray
nodes = build_nodes(X₀, u⁰, u̇⁰, ü⁰, w⁰, ẇ⁰, ẅ⁰, nothing, R⁰)

# geometric and material properties
E = 2.6 * 1E5 * 1E4
ν = 0.33
ρ = 8 * 1E-9 * 1E12
damping = 50

Re₀ = read_ics_mat(readdlm("stent/output3D/outputCrimping3D/" * "Re0.txt"))
# Re₀ = reshape(reinterpret(Int64, Re₀), (9, length(Re₀)))'

beams = build_beams(nodes, connectivity, E, ν, ρ, rWireSection, damping, Re₀)

# -------------------------------------------------------------------------------------------
# Contact parameters
# -------------------------------------------------------------------------------------------

# contact parameters
kₙ = 4/3 * 5/(1-0.5^2)*sqrt(rWireSection) # Approximate Hertz contact with 5 MPa wall stiffness
μ = 0.01
εᵗ = 0.001 #regularized parameter for friction contact
ηₙ = 0.1
kₜ = kₙ
ηₜ = ηₙ
u̇ₛ = 0.005
beam2beam = true

contact = ContactParameters(kₙ, μ, εᵗ, ηₙ, kₜ, ηₜ, u̇ₛ, beam2beam)
# sdf = Cylinder_SDF(rWireSection, 2)
sdf = nothing

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

# Dirichlet boundary conditions: fixed positions
ndofs = nnodes*6
fixed_dofs = Int[]
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

conf = Configuration(nodes, beams, nothing, ext_forces, bcs, contact, sdf)

# -------------------------------------------------------------------------------------------
# Time stepping parameters
# -------------------------------------------------------------------------------------------

# initial time step and total time
ini_Δt = 1e-2
max_Δt = 1e2
Δt_plot = 1e-2
tᵉⁿᵈ = 10000000

params = Params(;ini_Δt, max_Δt, Δt_plot, tᵉⁿᵈ, output_dir = "stent/output3D/outputDeployment3D", verbose=true , record_timings=true)

# -------------------------------------------------------------------------------------------
# Solve
# -------------------------------------------------------------------------------------------

solver!(conf, params);

