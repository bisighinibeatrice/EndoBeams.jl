#----------------------------------
# STRUCTURES
#----------------------------------

struct ContactParameters
    kₙ::Float64
    μ::Float64
    εᵗ::Float64
    ηₙ::Float64
    kₜ::Float64
    ηₜ::Float64
    u̇ₛ::Float64
    beam2beam::Bool
end


function ContactParameters(kₙ, μ, εᵗ, ηₙ, kₜ, ηₜ, u̇ₛ, beam2beam=false)

    return ContactParameters(kₙ, μ, εᵗ, ηₙ, kₜ, ηₜ, u̇ₛ, beam2beam)

end

# Dirichlet boundary conditions
struct BoundaryConditions{TF}
    
    # blocked dofs 
    fixed_dofs::Vector{Int}
    
    # free dofs 
    free_dofs::Vector{Int}
    
    # displacement function
    u::TF
    
    # displacement vector
    disp_vals::Vector{Float64}
    
    # dofs where the displacement function is applied
    disp_dofs::Vector{Int}   
    
    # cylindrical coordinates for the bcs
    flag_cylindrical::Bool
    
    # load bcs from files containing the displacements value
    flag_load_from_file::Bool
    
    # folder with the files containing the displacements value
    dir_folder_load::String
    
end


"""
bcs = BoundaryConditions(fixed_dofs, free_dofs, flag_cylindrical, flag_disp_vector, Fdisp, disp_vals, disp_dofs)

Constructor of the structure containing the information about the Dirichlet boundary conditions (BCs), if present:
- `fixed_dofs`: fixed DOFs;
- `free_dofs`: free DOFs;
- `flag_cylindrical`: true if the BCs are defined in cylindral coordinates;
- `flag_disp_vector`: true if the BCs are defined as a vector;
- `Fdisp`: imposed dispacement amplitude;
- `disp_vals`: vector; containing the imposed dispacements;
- `disp_dofs`: DOFs interested by the Dirichlet BCs.

Returns a BoundaryConditions structure.
"""
function BoundaryConditions(fixed_dofs, free_dofs, disp_fun::TF, disp_vals, disp_dofs, flag_cylindrical=false, flag_load_from_file=false, dir_folder_load="") where TF
    
    return BoundaryConditions{TF}(fixed_dofs, free_dofs, disp_fun, disp_vals, disp_dofs, flag_cylindrical, flag_load_from_file, dir_folder_load)
    
end 



# External forces
struct ExternalForces{TF}
    
    # concentrated force function 
    f::TF
    
    # dofs where the concentrated force is applied
    loaded_dofs::Vector{Int}
    
end 


"""
fᵉˣᵗ = ExternalForces(F, loaded_dofs)

Constructor of the structure containing the information about the external load, if present:
- `F`: external load amplitude;
- `loaded_dofs`: DOFs interested by the external load.

Returns a ExternalForces structure.
"""
function ExternalForces(force_fun::TF, loaded_dofs) where TF
    
    return ExternalForces{TF}(force_fun, loaded_dofs)
    
end 



# Material properties
struct BeamProperties{TK, TJ}
    
    radius::Float64
    E::Float64
    K̄ⁱⁿᵗ::TK
    Jᵨ::TJ
    Aᵨ::Float64
    damping::Float64
    
end 

"""
material = Material(E, ν, ρ, radius, damping)

Constructor of the structure containing the material properties:
- `E`: Young modulus;
- `ν`: Poisson coefficient;
- `ρ`: density;
- `radius`: beam radius;
- `damping`: viscous damping

Returns a Material structure.
"""
function BeamProperties(l₀, E, ν, ρ, radius, damping)
    
    G = E/(2*(1+ν))
    A = pi*radius^2
    I₂₂ = pi*radius^4/4
    I₃₃ = I₂₂
    Iₒ = I₂₂ + I₃₃
    Jᵨ = ρ * Diagonal(Vec3(Iₒ, I₂₂, I₃₃))
    Aᵨ = ρ*A
    
    K̄ⁱⁿᵗ = K̄ⁱⁿᵗ_beam(E, G, Iₒ, A, I₂₂, I₃₃, l₀)
    
    beamprops = BeamProperties{typeof(K̄ⁱⁿᵗ), typeof(Jᵨ)}(radius, E, K̄ⁱⁿᵗ, Jᵨ, Aᵨ, damping)
    
    return beamprops
    
end


# Configuration of the mesh
struct Configuration{Tn, Tb, Tc, Te, Tbc, Tcon, Tsdf, Tgps}
    
    nodes::Tn
    beams::Tb
    constraints::Tc
    
    # dofs
    ndofs::Int
    disp_dofs::Vector{Int}
    rot_dofs::Vector{Int}
    
    # external forces
    ext_forces::Te
    
    # boundary conditions
    bcs::Tbc
    
    contact::Tcon
    
    sdf::Tsdf
    
    colors::Vector{Vector{Int}}
    
    gausspoints::Tgps
    
end



function first_available(color_set)
    # Return smallest positive integer not in the given list of colors.
    count = 1
    while true
        if count ∉ color_set
            return count
        end
        count += 1
    end
end




function greedy_color(beams)
    
    elcolors = zeros(Int, length(beams))
    colors = Vector{Vector{Int}}()
    for (i,b) in enumerate(LazyRows(beams))
        neighbors = [j for (j, nodes) in enumerate(zip(beams.node1, beams.node2)) if (b.node1 in nodes) || (b.node2 in nodes)]
        used_neighbour_colors = [elcolors[nbr] for nbr in neighbors if elcolors[nbr]>0]
        c = first_available(used_neighbour_colors)
        elcolors[i] = c
        if c>length(colors)
            push!(colors, [i])
        else
            push!(colors[c], i)
        end
        
    end
    
    return colors
    
end

struct GaussPoint
    
    pos::Vec3{Float64}
    δₜ::Float64
    status::Int
    
end 

"""
conf = Configuration(material, geometry, nnodes, ndofs, ext_forces, bcs)

Constructor of the structure collecting the information for the simulation:
- `material`: material properties (Material{Float64});
- `geometry`: geoltrical properties (Geometry{Float64});
- `nnodes`: number of nodes in the system;
- `ndofs`: number of DOFs in the system;
- `bcs`: Dirichlet BCs.

Returns a Configuration structure.
"""
function Configuration(nodes::StructVector, beams::StructVector, constraints::Union{StructVector, Nothing}, ext_forces::ExternalForces, bcs::BoundaryConditions, contact::Union{ContactParameters, Nothing}, sdf::Union{SignedDistanceField, Nothing})
    
    ndofs = length(nodes)*6
    
    disp_dofs = [i for i in 1:ndofs if mod1(i, 6)≤3]
    rot_dofs = [i for i in 1:ndofs if mod1(i, 6)>3]
    
    ngps = length(beams)*3
    gausspoints = StructArray(GaussPoint(
    (0,0,0),
    0, 
    0) for i in 1:ngps)
    
    return Configuration{typeof(nodes), typeof(beams), typeof(constraints), typeof(ext_forces), typeof(bcs), typeof(contact), typeof(sdf), typeof(gausspoints)}(nodes, beams, constraints, ndofs, disp_dofs, rot_dofs, ext_forces, bcs, contact, sdf, greedy_color(beams), gausspoints)
    
end 

# Force vectors needed at the next time step by the solver
struct Solution 
    
    fᵉˣᵗ::Vector{Float64}
    Tⁱⁿᵗ::Vector{Float64}
    Tᵏ::Vector{Float64}
    Tᶜ::Vector{Float64}
    Tᶜᵒⁿ::Vector{Float64}
    
end 

# Constructor of the structure where the last force vectors are saved in order to be used in the next step by the solver
function Solution(conf::Configuration)
    
    ndofs = conf.ndofs
    
    # external force @t=0
    fᵉˣᵗ  = zeros(ndofs)
    for i in conf.ext_forces.loaded_dofs
        fᵉˣᵗ[i] = conf.ext_forces.f(0, i)
    end
    
    Tᵢₙₜ = zeros(ndofs)
    Tₖ = zeros(ndofs)
    Tₜ = zeros(ndofs)
    Tₓ = zeros(ndofs)
    
    return Solution(fᵉˣᵗ, Tᵢₙₜ, Tₖ, Tₜ, Tₓ)
    
end 

# Current nodal solutions (preallocation)
struct NodalSolution 
    
    D::Vector{Float64}
    Ḋ::Vector{Float64}
    D̈::Vector{Float64}
    
    r::Vector{Float64}
    Ktan::Matrix{Float64}
    ΔD::Vector{Float64}    
    temp::Vector{Float64}
    
    r_free::Vector{Float64}
    Ktan_free::Matrix{Float64}
    ΔD_free::Vector{Float64}
    
end 


# Constructor of the structure containing the preallocated variables used in the solver
function NodalSolution(ndofs, nfreedofs)

    return NodalSolution(
    
    zeros(ndofs), #D
    zeros(ndofs), #Ḋ 
    zeros(ndofs), #D̈ 
    
    zeros(ndofs), #r
    zeros(ndofs, ndofs), #Ktan
    zeros(ndofs), #ΔD    
    zeros(ndofs), #temp vector
    
    zeros(nfreedofs), #r_free
    zeros(nfreedofs, nfreedofs), #Ktan_free
    zeros(nfreedofs)) #ΔD_free
    
end


# Global matrices structure 
struct Matrices
    
    K ::Matrix{Float64}
    C::Matrix{Float64}
    M::Matrix{Float64}
    
    Tⁱⁿᵗ::Vector{Float64}
    Tᵏ::Vector{Float64}
    Tᵛ::Vector{Float64}
    Tᶜ::Vector{Float64}
    Tᶜᵒⁿ::Vector{Float64}
        
end 

# Constructor of the structure containing the global matrices 
function Matrices(ndofs)
    
    C = zeros(ndofs, ndofs)
    M = zeros(ndofs, ndofs)
    K =  zeros(ndofs, ndofs)
    
    Tⁱⁿᵗ = zeros(ndofs)
    Tᵏ =  zeros(ndofs)
    Tᵛ =  zeros(ndofs)
    Tᶜ = zeros(ndofs)
    Tᶜᵒⁿ = zeros(ndofs)
    
    return Matrices(K, C, M, Tⁱⁿᵗ, Tᵏ, Tᵛ, Tᶜ, Tᶜᵒⁿ)
    
end 

# Energy contributions structure
mutable struct Energy
    
    strain_energy::Float64
    kinetic_energy::Float64
    contact_energy::Float64
    
end

# Constructor of the structure containing the energy contributions
function Energy()
    
    return Energy(0, 0, 0)
    
end

struct VTKData
    
    VTKcollection::WriteVTK.CollectionFile
    output_dir::String
    
    intermediate_points::Int
    
    interpolated_points::Vector{Vec3{Float64}}
    interpolated_lines::Vector{MeshCell{WriteVTK.PolyData.Lines, UnitRange{Int}}}
    
    stress::Vector{Float64}
    strain::Vector{Float64}
    displacement::Vector{Vec3{Float64}}
    velocity::Vector{Vec3{Float64}}
    
    contact_distance::Vector{Float64}
    normal_contact_force::Vector{Vec3{Float64}}
    tangential_contact_force::Vector{Vec3{Float64}}
    tangential_velocity::Vector{Float64}
    incontact::Vector{Int}
    status::Vector{Int}
    
end

function VTKData(nbeams, output_dir, sdf, intermediate_points = 5)
    
    interpolated_points = zeros(Vec3{Float64}, nbeams*intermediate_points)
    
    interpolated_lines = [MeshCell(PolyData.Lines(), (i-1)*intermediate_points+1:i*intermediate_points) for i in 1:nbeams]
    
    stress = zeros(length(interpolated_points))
    strain = zeros(length(interpolated_points))
    displacement = zeros(Vec3{Float64}, length(interpolated_points))
    velocity = zeros(Vec3{Float64}, length(interpolated_points))
    
    if !isnothing(sdf)
        contact_distance = zeros(length(interpolated_points))
        normal_contact_force = zeros(Vec3{Float64}, length(interpolated_points))
        tangential_contact_force = zeros(Vec3{Float64}, length(interpolated_points))
        tangential_velocity = zeros(length(interpolated_points))
        incontact = zeros(Int, length(interpolated_points))
        status = zeros(Int, length(interpolated_points))
    else
        contact_distance = Float64[]
        normal_contact_force = Vec3{Float64}[]
        tangential_contact_force = Vec3{Float64}[]
        tangential_velocity = Float64[]
        incontact = Int[]
        status = Int[]
    end
    
    
    clean_folders(output_dir)
    
    collection = paraview_collection("$output_dir/simulation")
    
    return VTKData(collection, output_dir, intermediate_points, interpolated_points, interpolated_lines, stress, strain, displacement, velocity, contact_distance, normal_contact_force, tangential_contact_force, tangential_velocity, incontact, status)
    
end


