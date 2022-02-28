#----------------------------------
# STRUCTURES
#----------------------------------

# Variables for the time and space integration and contact linearization
struct SimulationParameters{T, NG}
    
    # integration parameters
    α::T 
    β::T 
    γ::T 
    
    # numerical damping
    damping::T
    
    # time step and total time
    Δt::Float64
    Δt_plot::Float64
    tᵉⁿᵈ::Float64
    
    # tolerance and maximum number of iterations
    tol_res::Float64
    tol_ΔD::Float64
    max_it::Int
    
    # Gauss points
    nᴳ::Int
    ωᴳ::SVector{NG, T}
    zᴳ::SVector{NG, T}
    
    # penalty parameters
    εᶜ::T
    μ::T
    εᵗ::T
    γᵈᵃᵐᵖ::T
    
end


"""
    comp = SimulationParameters(α, β, γ, damping, Δt, Δt_plot, tᵉⁿᵈ, tol_res, tol_ΔD, max_it, nG, ωG, zG, εᶜ, μ, εᵗ, T=Float64)    

Constructor of the structure containing the simulation parameters:
    - `α`: integration parameter;
    - `β`: integration parameter;
    - `γ`: integration parameter;
    - `damping`: damping coefficient;
    - `Δt`: time step;
    - `Δt_plot`: time step for the saving of output files;
    - `tᵉⁿᵈ`: total time of the simulation;
    - `tol_res`: residual tolerance for the Newton-Raphson algorithm;
    - `tol_ΔD`: solution vector tolerance for the Newton-Raphson algorithm;
    - `max_it`: maximum number of iterations for the Newton-Raphson algorithm;
    - `nG`: number of Gauss points;
    - `ωG`: Gauss points weights;
    - `zG`: Gauss points positions along centreline;
    - `εᶜ`: penalty coefficient for the contact normal contributions;
    - `μ`: friction coefficient;
    - `εᵗ`: regularisation coefficient for the contact tangential contributions.
    - `γᵈᵃᵐᵖ`: contact normal damping parameter

Returns a SimulationParameters structure.
"""
function SimulationParameters(α, β, γ, damping, Δt, Δt_plot, tᵉⁿᵈ, tol_res, tol_ΔD, max_it, nG, ωG, zG, εᶜ, μ, εᵗ, γᵈᵃᵐᵖ, T=Float64)
    
    return SimulationParameters{T, nG}(α, β, γ, damping, Δt, Δt_plot, tᵉⁿᵈ, tol_res, tol_ΔD, max_it, nG, ωG, zG, εᶜ, μ, εᵗ, γᵈᵃᵐᵖ)
    
end 

# Dirichlet boundary conditions
struct BoundaryConditions{T, TF}
    
    # blocked dofs 
    fixed_dofs::Vector{Int}
    
    # free dofs 
    free_dofs::Vector{Int}
    
    # displacement function
    u::TF

    # displacement vector
    disp_vals::Vector{T}

    # dofs where the displacement function is applied
    disp_dofs::Vector{Int}   

end


"""
    bcs = BoundaryConditions(fixed_dofs, free_dofs, flag_cylindrical, flag_disp_vector, Fdisp, disp_vals, disp_dofs, T=Float64)

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
function BoundaryConditions(fixed_dofs, free_dofs, disp_fun::TF, disp_vals, disp_dofs, T=Float64) where TF

    return  BoundaryConditions{T, TF}(fixed_dofs, free_dofs, disp_fun, disp_vals, disp_dofs)
    
end 



# External forces
struct ExternalForces{TF}
    
    # concentrated force function 
    f::TF
    
    # dofs where the concentrated force is applied
    loaded_dofs::Vector{Int}
    
end 


"""
    fᵉˣᵗ = ExternalForces(F, loaded_dofs, T=Float64)

Constructor of the structure containing the information about the external load, if present:
- `F`: external load amplitude;
- `loaded_dofs`: DOFs interested by the external load.

Returns a ExternalForces structure.
"""
function ExternalForces(force_fun::T, loaded_dofs) where T
        
    return ExternalForces{T}(force_fun, loaded_dofs)
    
end 



# Material properties
struct Material{T, TJ}
    
    E::T
    G::T
    Aᵨ::T
    Jᵨ::TJ
    
end 

"""
    mat = Material(E, nu, ρ, radius)

Constructor of the structure containing the material properties:
- `E`: Young modulus;
- `nu`: Poisson coefficient;
- `ρ`: density;
- `radius`: beam radius.

Returns a Material structure.
"""
function Material(E, nu, ρ, radius, T=Float64)

    G = E/(2*(1+nu))
    A = pi*radius^2
    I₂₂ = pi*radius^4/4
    I₃₃ = pi*radius^4/4
    Iₒ = I₂₂ + I₃₃
    Jᵨ = ρ * Diagonal(Vec3(Iₒ, I₂₂, I₃₃))
    Aᵨ = ρ*A

    mat = Material{T, typeof(Jᵨ)}(E, G, Aᵨ, Jᵨ)

    return mat

end

# Geometrical properties
struct Geometry{T}
    
    A::T
    I₂₂::T
    I₃₃::T
    Iₒ::T
    Iᵣᵣ::T
    J::T   
    
end 


"""
    geom = Geometry(radius)

Constructor of the structure containing the geometrical properties:
- `radius`: beam radius.

Returns a Geometry structure.
"""
function Geometry(radius, T=Float64)

    A = pi*radius^2
    I₂₂ = pi*radius^4/4
    I₃₃ = pi*radius^4/4
    Iₒ = I₂₂+I₃₃
    Iᵣᵣ = Iₒ
    J = Iₒ
    
    geom = Geometry{T}(A, I₂₂, I₃₃, Iₒ, Iᵣᵣ, J)

    return geom

end 

# Configuration of the mesh
struct Configuration{T, TJ, TF, TU}
    
    # material
    mat::Material{T, TJ}
    
    # geometry
    geom::Geometry{T}
    
    # dof
    ndofs::Int
    disp_dofs::Vector{Int}
    rot_dofs::Vector{Int}
    
    # external forces
    ext_forces::ExternalForces{TF}
    
    # boundary conditions
    bcs::BoundaryConditions{T, TU}
    
end




"""
    conf = Configuration(mat, geom, nnodes, ndofs, ext_forces, bcs, T=Float64)

Constructor of the structure collecting the information for the simulation:
- `mat`: material properties (Material{T});
- `geom`: geoltrical properties (Geometry{T});
- `nnodes`: number of nodes in the system;
- `ndofs`: number of DOFs in the system;
- `bcs`: Dirichlet BCs.

Returns a Configuration structure.
"""
function Configuration(mat::Material{TT, TJ}, geom::Geometry, ndofs, ext_forces::ExternalForces{TF}, bcs::BoundaryConditions{TT, TU}, T=Float64) where {TT, TF, TU, TJ}
    
    # displacement dofs
    disp_dofs = [i for i in 1:ndofs if mod1(i, 6)≤3]
    rot_dofs = [i for i in 1:ndofs if mod1(i, 6)>3]

    return Configuration{T, TJ, TF, TU}(mat, geom, ndofs, disp_dofs, rot_dofs, ext_forces, bcs)
    
end 

# Force vectors needed at the next time step by the solver
struct Solution{T}  
    
    fᵉˣᵗ::Vector{T}
    Tⁱⁿᵗ::Vector{T}
    Tᵏ::Vector{T}
    Tᶜ::Vector{T}
    Tᶜᵒⁿ::Vector{T}
    
end 

# Constructor of the structure where the last force vectors are saved in order to be used in the next step by the solver
function Solution(conf, T=Float64)
    
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
    
    return Solution{T}(fᵉˣᵗ, Tᵢₙₜ, Tₖ, Tₜ, Tₓ)
    
end 

# Current nodal solutions (preallocation)
struct NodalSolution{T}  
    
    D::Vector{T}
    Ḋ::Vector{T}
    D̈::Vector{T}
    
    r::Vector{T}
    Ktan_mat::SparseMatrixCSC{T,Int}
    Ktan::Vector{Float64}
    ΔD::Vector{T}    
    temp::Vector{T}

    r_free::Vector{Float64}
    Ktan_free::SparseMatrixCSC{Float64,Int}
    ΔD_free::Vector{Float64}
    
end 


# Constructor of the structure containing the preallocated variables used in the solver
function NodalSolution(Ktan, Ktan_free, ndofs, nfreedofs, T)

    return NodalSolution{T}(

    zeros(T, ndofs), #D
    zeros(T, ndofs), #Ḋ 
    zeros(T, ndofs), #D̈ 
    
    zeros(T, ndofs), #r
    Ktan, #Ktan_mat
    nonzeros(Ktan), #Ktan
    zeros(T, ndofs), #ΔD    
    zeros(T, ndofs), #temp vector

    zeros(T, nfreedofs), #r_free
    Ktan_free, #Ktan_free
    zeros(T, nfreedofs)) #ΔD_free

end


# Global matrices structure (sparse arrays)
struct Matrices{T}
    
    K_mat ::SparseMatrixCSC{T,Int}
    K::Vector{T}
    C_mat::SparseMatrixCSC{T,Int}
    C::Vector{T}
    M_mat::SparseMatrixCSC{T,Int}
    M::Vector{T}
    
    Tⁱⁿᵗ::Vector{T}
    Tᵏ::Vector{T}
    Tᶜ::Vector{T}
    Tᶜᵒⁿ::Vector{T}

    sparsity_free::Vector{Int}
    
end 

# Constructor of the structure containing the global matrices 
function Matrices(I, J, sparsity_free, ndofs, T=Float64)
    
    C_mat = sparse(I, J, zero(T))
    C = nonzeros(C_mat)
    
    M_mat = sparse(I, J, zero(T))
    M = nonzeros(M_mat)

    K_mat = sparse(I, J, zero(T))
    K =  nonzeros(K_mat)

    Tⁱⁿᵗ = zeros(T, ndofs)
    Tᵏ =  zeros(T, ndofs)
    Tᶜ = zeros(T, ndofs)
    Tᶜᵒⁿ = zeros(T, ndofs)

    return Matrices{T}(K_mat, K, C_mat, C, M_mat, M, Tⁱⁿᵗ, Tᵏ, Tᶜ, Tᶜᵒⁿ, sparsity_free)
    
end 

# Energy contributions structure
mutable struct Energy{T}
    
    strain_energy::T
    kinetic_energy::T
    contact_energy::T
    
end



# Constructor of the structure containing the energy contributions
function Energy(T=Float64)
    
    return Energy{T}(0, 0, 0)
    
end






struct VTKData{T}

    VTKcollection::WriteVTK.CollectionFile
    output_dir::String

    intermediate_points::Int

    interpolated_points::Vector{Vec3{T}}
    interpolated_lines::Vector{MeshCell{VTKCellType, Tuple{Int, Int}} }

    stress::Vector{T}
    strain::Vector{T}
    displacement::Vector{Vec3{T}}
    velocity::Vector{Vec3{T}}

    contact_distance::Vector{T}
    normal_contact_force::Vector{Vec3{T}}
    tangential_contact_force::Vector{Vec3{T}}
    incontact::Vector{Int}

end

function VTKData(nbeams, output_dir, sdf, T)

    intermediate_points = 5

    interpolated_points = zeros(Vec3{T}, nbeams*intermediate_points)
    interpolated_lines = [MeshCell(VTKCellTypes.VTK_LINE, ((i-1)*intermediate_points+j, (i-1)*intermediate_points+(j+1))) for i in 1:nbeams for j in 1:intermediate_points-1]


    stress = zeros(T, length(interpolated_points))
    strain = zeros(T, length(interpolated_points))
    displacement = zeros(Vec3{T}, length(interpolated_points))
    velocity = zeros(Vec3{T}, length(interpolated_points))

    if !isnothing(sdf)
        contact_distance = zeros(T, length(interpolated_points))
        normal_contact_force = zeros(Vec3{T}, length(interpolated_points))
        tangential_contact_force = zeros(Vec3{T}, length(interpolated_points))
        incontact = zeros(Int, length(interpolated_points))
    else
        contact_distance = T[]
        normal_contact_force = Vec3{T}[]
        tangential_contact_force = Vec3{T}[]
        incontact = Int[]
    end


    clean_folders(output_dir)

    collection = paraview_collection("$output_dir/simulation")

    return VTKData{T}(collection, output_dir, intermediate_points, interpolated_points, interpolated_lines, stress, strain, displacement, velocity, contact_distance, normal_contact_force, tangential_contact_force, incontact)

end







# Constructor of the sparse matrices
function constructor_sparse_matrices!(beams, nodes, constraints, conf::Configuration{T}) where T

    ndofs =  conf.ndofs
    free_dofs = conf.bcs.free_dofs
    fixed_dofs = conf.bcs.fixed_dofs
    nfreedofs = length(free_dofs)

    I, J, sparsity_free = sparsity(nodes, beams, constraints, fixed_dofs)

    Ktan = sparse(I, J, 0.)
    
    Ktan_free = Ktan[free_dofs, free_dofs]

    matrices = Matrices(I, J, sparsity_free, ndofs, T)
    nodes_sol = NodalSolution(Ktan, Ktan_free, ndofs, nfreedofs, T)

    return matrices, nodes_sol

end 

